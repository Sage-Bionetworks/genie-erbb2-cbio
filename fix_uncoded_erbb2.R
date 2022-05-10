# Description: Custom fixes to the ERBB2 uncoded REDCap data.
# - 1. Ensure all SAMPLE IDs do not contain spaces.
# - 2. Ensure all SAMPLE IDs start with "GENIE-"
# Author: Haley Hunter-Zinck
# Date: 2022-05-09

# setup ----------------------------

tic = as.double(Sys.time())

library(glue)
library(dplyr)
library(synapser)
library(optparse)

# constants
sample_nos <- c(1:10)

waitifnot <- function(cond, msg) {
  if (!cond) {
    
    for (str in msg) {
      message(str)
    }
    message("Press control-C to exit and try again.")
    
    while(T) {}
  }
}

# user input ----------------------------

option_list <- list( 
  make_option(c("-i", "--synid_file_input"), type = "character",
              help="Synapse ID of input file"),
  make_option(c("-o", "--synid_folder_output"), type = "character", default = NA,
              help="Synapse ID of output folder"),
  make_option(c("-f", "--file_name"), type = "character", default = "fixed.csv", 
              help="Name of the output file (default: 'fixed.csv')"),
  make_option(c("-v", "--verbose"), action="store_true", default = FALSE, 
              help="Output script messages to the user."),
  make_option(c("-a", "--auth"), 
              type = "character",
              default = NA,
              help="Synapse personal access token or path to .synapseConfig (default: normal synapse login behavior)")
)
opt <- parse_args(OptionParser(option_list=option_list))
waitifnot(!is.null(opt$synid_file_input) && !is.null(opt$synid_folder_output),
          msg = "Rscript fix_uncoded_erbb2.R -h")

synid_file_input <- opt$synid_file_input
synid_folder_output <- opt$synid_folder_output
file_output <- opt$file_name
verbose <- opt$verbose
auth <- opt$auth

if (verbose) {
  print(glue("Parameters:"))
  print(glue("----------"))
  print(glue("synid_file_input:\t{synid_file_input}"))
  print(glue("synid_folder_output:\t{synid_folder_output}"))
  print(glue("file_output:\t'{file_output}'"))
  print(glue("verbose:\t\t{verbose}"))
  print(glue("----------"))
}

# functions ----------------------------

#' Return current time as a string.
#' 
#' @param timeOnly If TRUE, return only time; otherwise return date and time
#' @param tz Time Zone
#' @return Time stamp as string
#' @example 
#' now(timeOnly = T)
now <- function(timeOnly = F, tz = "US/Pacific") {
  
  Sys.setenv(TZ=tz)
  
  if(timeOnly) {
    return(format(Sys.time(), "%H:%M:%S"))
  }
  
  return(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}

#' Extract personal access token from .synapseConfig
#' located at a custom path. 
#' 
#' @param path Path to .synapseConfig
#' @return personal acccess token
get_auth_token <- function(path) {
  
  lines <- scan(path, what = "character", sep = "\t", quiet = T)
  line <- grep(pattern = "^authtoken = ", x = lines, value = T)
  
  token <- strsplit(line, split = ' ')[[1]][3]
  return(token)
}

#' Override of synapser::synLogin() function to accept 
#' custom path to .synapseConfig file or personal authentication
#' token.  If no arguments are supplied, performs standard synLogin().
#' 
#' @param auth full path to .synapseConfig file or authentication token
#' @param silent verbosity control on login
#' @return TRUE for successful login; F otherwise
synLogin <- function(auth = NA, silent = T) {
  
  secret <- Sys.getenv("SCHEDULED_JOB_SECRETS")
  if (secret != "") {
    # Synapse token stored as secret in json string
    syn = synapser::synLogin(silent = T, authToken = fromJSON(secret)$SYNAPSE_AUTH_TOKEN)
  } else if (auth == "~/.synapseConfig" || is.na(auth)) {
    # default Synapse behavior
    syn <- synapser::synLogin(silent = silent)
  } else {
    
    # in case pat passed directly
    token <- auth
    
    # extract token from custom path to .synapseConfig
    if (grepl(x = auth, pattern = "\\.synapseConfig$")) {
      token = get_auth_token(auth)
      
      if (is.na(token)) {
        return(F)
      }
    }
    
    # login with token
    syn <- tryCatch({
      synapser::synLogin(authToken = token, silent = silent)
    }, error = function(cond) {
      return(F)
    })
  }
  
  # NULL returned indicates successful login
  if (is.null(syn)) {
    return(T)
  }
  return(F)
}

#' Download and load data stored in csv or other delimited format on Synapse
#' into an R data frame.
#' 
#' @param synapse_id Synapse ID
#' @version Version of the Synapse entity to download.  NA will load current
#' version
#' @param set Delimiter for file
#' @param na.strings Vector of strings to be read in as NA values
#' @param header TRUE if the file contains a header row; FALSE otherwise.
#' @param check_names TRUE if column names should be modified for compatibility 
#' with R upon reading; FALSE otherwise.
#' @param comment.char character designating comment lines to ignore
#' @return data frame
get_synapse_entity_data_in_csv <- function(synapse_id, 
                                           version = NA,
                                           sep = ",", 
                                           na.strings = c("NA"), 
                                           header = T,
                                           check_names = F,
                                           comment.char = "#",
                                           colClasses = "character") {
  
  if (is.na(version)) {
    entity <- synGet(synapse_id)
  } else {
    entity <- synGet(synapse_id, version = version)
  }
  
  data <- read.csv(entity$path, stringsAsFactors = F, 
                   na.strings = na.strings, sep = sep, check.names = check_names,
                   header = header, comment.char = comment.char, colClasses = colClasses)
  return(data)
}

#' Store a file on Synapse with options to define provenance.
#' 
#' @param path Path to the file on the local machine.
#' @param parent_id Synapse ID of the folder or project to which to load the file.
#' @param file_name Name of the Synapse entity once loaded
#' @param prov_name Provenance short description title
#' @param prov_desc Provenance long description
#' @param prov_used Vector of Synapse IDs of data used to create the current
#' file to be loaded.
#' @param prov_exec String representing URL to script used to create the file.
#' @return Synapse ID of entity representing file
save_to_synapse <- function(path, 
                            parent_id, 
                            file_name = NA, 
                            prov_name = NA, 
                            prov_desc = NA, 
                            prov_used = NA, 
                            prov_exec = NA) {
  
  if (is.na(file_name)) {
    file_name = path
  } 
  file <- File(path = path, parentId = parent_id, name = file_name)
  
  if (!is.na(prov_name) || !is.na(prov_desc) || !is.na(prov_used) || !is.na(prov_exec)) {
    act <- Activity(name = prov_name,
                    description = prov_desc,
                    used = prov_used,
                    executed = prov_exec)
    file <- synStore(file, activity = act)
  } else {
    file <- synStore(file)
  }
  
  return(file$properties$id)
}

#' If a string does not begin with the prefix "GENIE-", add it to
#' the beginning of the string.
#' 
#' @param x string
#' @return string with prefix "GENIE-"
add_genie_prefix <- function(x) {
  if (!(is.na(x) || grepl(pattern = "^GENIE-", x = x))) {
    return(glue("GENIE-{x}"))
  }
  return(x)
}

# synapse login --------------------

status <- synLogin(auth = auth)

# read ----------------------------

if (verbose) { print(glue("{now()}: reading raw uncoded data...")) }

df_raw <- get_synapse_entity_data_in_csv(synid_file_input, na.strings = c("NA",""))

# main ----------------------------

df_fix <- df_raw

for (sample_no in sample_nos) {
  varname <- glue("sample_id_{sample_no}")
  mod <- df_fix[[varname]]
  mod <- gsub(pattern = "[[:space:]]+", replacement = "", x = mod)
  mod <- as.character(unlist(sapply(mod, add_genie_prefix)))
  df_fix[varname] <- mod
}

# write --------------------

write.csv(x = df_fix, file = file_output, row.names = F, na = "", quote = F)

if (!is.na(synid_folder_output)) {
  
  if (verbose) {
    print(glue("{now(timeOnly = T)}: Saving uncoded data to Synapse..."))
  }
  
  synid_file_output <- save_to_synapse(path = file_output,
                                       parent_id = synid_folder_output,
                                       prov_name = "Fixed REDCap export",
                                       prov_desc = "REDCap export data with variable values uncoded according to the data dictionary and custom fixes applied",
                                       prov_used = c(synid_file_input),
                                       prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/main/fix_uncoded_erbb2.R")
  
  # clean up locally
  file.remove(file_output)
  
  if (verbose) {
    print(glue("{now(timeOnly = T)}: Saved file to Synapse at {synid_file_output} and removed local file. "))
  }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
