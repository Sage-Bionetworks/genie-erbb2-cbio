# Description: Transform uncoded ERBB2 REDCap export to cBioPortal format
#               to create patients, sample, and timeline files. 
# Author: Haley Hunter-Zinck
# Date: 2022-04-29

# setup ----------------------------

tic = as.double(Sys.time())

library(glue)
library(dplyr)
library(synapser)
library(optparse)

waitifnot <- function(cond, msg) {
  if (!cond) {
    
    for (str in msg) {
      message(str)
    }
    message("Press control-C to exit and try again.")
    
    while(T) {}
  }
}

# parameters
outfiles <- c()

# user input ----------------------------

option_list <- list( 
  make_option(c("-i", "--synid_file_input"), type = "character",
              help="Synapse ID of uncoded REDCap export"),
  make_option(c("-m", "--synid_table_map"), type = "character",
              help="Synapse ID of table of REDCap to cBio mapping"),
  make_option(c("-o", "--synid_folder_output"), type = "character",
              help="Synapse ID of output folder"),
  make_option(c("-v", "--verbose"), action="store_true", default = FALSE, 
              help="Output script messages to the user."),
  make_option(c("-a", "--auth"), 
              type = "character",
              default = NA,
              help="Synapse personal access token or path to .synapseConfig (default: normal synapse login behavior)")
)
opt <- parse_args(OptionParser(option_list=option_list))
waitifnot(!is.null(opt$synid_file_input) && !is.null(opt$synid_folder_output),
          msg = "Rscript template.R -h")

synid_file_input <- opt$synid_file_input
synid_table_map <- opt$synid_table_map
synid_folder_output <- opt$synid_folder_output
verbose <- opt$verbose
auth <- opt$auth

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

get_patient_os_months <- function(df_data) {
  os_months <- df_data$date_death_int - df_data$date_first_met_int
  return(os_months)
}

get_patient_os_status <- function(df_data) {
  
  os_status <- df_data$vital_status
  os_status[which(os_status =="Dead")] = "DECEASED"
  os_status[which(os_status =="Alive")] = "LIVING"
  return(os_status)
}

get_sample_seq_assay_id <- function(df_data, sample_no, synid_table_assay = "syn7517674") {
  var_sample_id <- glue("sample_id_{sample_no}")
  sample_ids <- paste0(df_data[var_sample_id], collapse = "','")
  query <- "SELECT SEQ_ASSAY_ID FROM {synid_table_assay} WHERE SAMPLE_ID IN '({sample_ids})'"
  seq_assay_id <- synTableQuery(query)
  return(seq_assay_id)
}

get_timeline_specimen_start_date <- function(df_data, sample_no) {
  start_date <- df_data$primary_dx_date_int - df_data[glue("seq_report_date_int_{sample_no}")]
  return(start_date)
}

get_fxn <- function(sampleType, code) {
  
  map_fxn_str <- c()
  map_fxn_str["get_patient_os_months"] <- get_patient_os_months()
  map_fxn_str["get_patient_os_status"] <- get_patient_os_status()
  map_fxn_str["get_sample_seq_assay_id"] <- get_sample_seq_assay_id()
  map_fxn_str["get_timeline_specimen_start_date"] <- get_timeline_specimen_start_date()
  
  label <- glue("get_{tolower(gsub('-', '_', sampleType))}_{tolower(code)}")
  
  return(map_fxn_str[label])
}

create_patient <- function(df_data, cbio, code, instructions) {
  
  mat <- matrix(NA, nrow = nrow(df_data), ncol = length(cbio), dimnames = list(c(), cbio))
  
  for (i in 1:length(cbio)) {
    if (code == "ETL_OPERATION") {
      fxn <- get_fxn(sampleType = "PATIENT", code = code[i])
      mat[,cbio[i]] <- fxn(df_data)
    } else {
      mat[,cbio[i]] <- df_data[code[i]]
    }
  }
  
  return(data.frame(mat))
}

create_sample <- function(df_data, cbio, code, instructions, sample_nos = c(1:10)) {
  
  mat <- matrix(nrow = 0, ncol = length(cbio), dimnames = list(c(), cbio))
  
  for (sample_no in sample_nos) {
    
    code_no <- sapply(glue, code)
    mat_no <- matrix(nrow = nrow(df_data), ncol = length(cbio), dimnames = list(c(), cbio))
    
    for (i in 1:length(cbio)) {
      if (code == "ETL_OPERATION") {
        fxn <- get_fxn(sampleType = "SAMPLE", code = code[i])
        mat[,cbio[i]] <- fxn(df_data)
      } else {
        mat_no[,cbio[i]] <- df_data[code_no[i]]
      }
    }
    
    mat <- rbind(mat, mat_no)
  }
  
  return(data.frame(mat))
}

create_timeline_treatment <- function(df_data, cbio, code, instructions, therapy_nos = c(1:25), drug_nos = c(1:5)) {
  
  mat <- matrix(nrow = 0, ncol = length(cbio), dimnames = list(c(), cbio))
  
  for (therapy_no in therapy_nos) {
    for (drug_no in drug_nos) {
      
      code_no <- sapply(glue, code)
      mat_no <- matrix(nrow = nrow(df_data), ncol = length(cbio), dimnames = list(c(), cbio))
      
      for (i in 1:length(cbio)) {
        if (code == "ETL_OPERATION") {
          fxn <- get_fxn(sampleType = "SAMPLE", code = code[i])
          mat[,cbio[i]] <- fxn(df_data)
        } else {
          mat_no[,cbio[i]] <- df_data[code_no[i]]
        }
      }
      
      mat <- rbind(mat, mat_no)
    }
  }
  
  return(data.frame(mat))
}

create_timeline_specimen <- function(df_data, cbio, code, instructions) {
  df_cbio <- NULL
  return(df)
}

create_timeline_status <- function(df_data, cbio, code, instructions) {
  df_cbio <- NULL
  return(df)
}

create_df <- function(df_data, df_map, sampleType) {
  
  df_data_type <- NULL
  
  df_map_type <- mapping %>% filter(sampleType == sampleType)
  cbio = df_map_type$cbio
  code = df_map_type$code
  instructions = df_map_type$instructions
  
  if (sampleType == "PATIENT") {
    df_data_type <- create_patient(df_data, cbio, code, instructions)
  } else if (sampleType == "SAMPLE") {
    df_data_type <- create_sample(df_data, cbio, code, instructions)
  } else if (sampleType == "TIMELINE-TREATMENT") {
    df_data_type <- create_timeline_treatment(df_data, cbio, code, instructions)
  } else if (sampleType == "TIMELINE-SPECIMEN") {
    df_data_type <- create_timeline_specimen(df_data, cbio, code, instructions)
  } else if (sampleType == "TIMELINE-STATUS") {
    df_data_type <- create_timeline_status(df_data, cbio, code, instructions)
  } else {
    message(glue("ERROR: {sampleType} not implemented.  Quitting..."))
    stop()
  }
  
  return(df_data_type)
}

write_patient <- function(df, description, colType, filename = "data_clinical_patient.txt") {
  return(filename)
}

write_sample <- function(df, description, colType, filename = "data_clinical_sample.txt") {
  return(filename)
}

write_timeline <- function(dfs, filename = "data_timeline.txt") {
  return(filename)
}

# synapse login --------------------

synLogin(auth = auth)

# read ----------------------------

if (verbose) {
  print(glue("{now(timeOnly=T)}: reading data and mapping..."))
}

df_data <- get_synapse_entity_data_in_csv(synid_file_data)

query <- "SELECT * FROM {synid_table_map}"
df_map <- as.data.frame(synTableQuery(query, includeRowIdAndVersion = F))

# main ----------------------------

if (verbose) {
  print(glue("{now(timeOnly=T)}: formatting data for cBioPortal files..."))
}

# create each data frame
for (sampleType in mapping %>% sampleType %>% distinct()) {
  dfs[[sampleType]] <- create_df(df_data = df_data, df_map = df_map, sampleType = sampleType)
}

if (verbose) {
  print(glue("{now(timeOnly=T)}: writing files locally..."))
}

# write data frames to file locally
outfile_pat <- write_patient(dfs[["PATIENT"]])
outfile_sam <- write_sample(dfs[["SAMPLE"]])
outfile_tl <- write_timeline(dfs[[grepl(names(dfs), value = T)]])

# store to synpase
if (!is.na(synid_folder_output)) {
  
  if (verbose) {
    print(glue("{now(timeOnly=T)}: storing files to Synapse in folder {synid_folder_output}..."))
  }
  
  synid_file_pat <- save_to_synapse(path = outfile_pat, 
                                    parent_id = synid_folder_output, 
                                    prov_name = "ERBB2 cbio patient file", 
                                    prov_desc = "GENIE ERBB2 cBioPortal patient clinical file", 
                                    prov_used = c(synid_file_data, synid_table_map), 
                                    prov_exec = "")
  
  synid_file_sam <- save_to_synapse(path = outfile_sam, 
                                    parent_id = synid_folder_output, 
                                    prov_name = "ERBB2 cbio sample file", 
                                    prov_desc = "GENIE ERBB2 cBioPortal sample clinical file", 
                                    prov_used = c(synid_file_data, synid_table_map), 
                                    prov_exec = "")
  
  synid_file_tl <- save_to_synapse(path = outfile_tl, 
                                    parent_id = synid_folder_output, 
                                    prov_name = "ERBB2 cbio timeline file", 
                                    prov_desc = "GENIE ERBB2 cBioPortal timeline file", 
                                    prov_used = c(synid_file_data, synid_table_map), 
                                    prov_exec = "")
}

# close out ----------------------------


toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
