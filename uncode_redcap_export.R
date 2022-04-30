# Description: Uncode a REDCap uploads for GENIE Sponsored 
#               Project given the export file and data dictionary.
# Author: Haley Hunter-Zinck
# Date: 2022-04-29

# setup ----------------------------

tic <- as.double(Sys.time())

library(optparse)
library(glue)
library(synapser)
library(dplyr)
library(lubridate)

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
  make_option(c("-e", "--synid_file_export"), type = "character",
              help=glue("Synapse ID of export file as csv or xlsx")),
  make_option(c("-d", "--synid_file_dictionary"), type = "character",
              help=glue("Synapse ID of data dictionary file as csv")),
  make_option(c("-o", "--synid_folder_output"), type = "character", default = NA, 
              help="Syanpse ID of folder in which to upload output file (default: write locally only)"),
  make_option(c("-f", "--file_name"), type = "character", default = "uncoded_redcap_export.csv", 
              help="Name of the output file (default: 'uncoded_redcap_export.csv')"),
  make_option(c("-a", "--synapse_auth"), type = "character", default = "~/.synapseConfig", 
              help="Path to .synapseConfig file or Synapse PAT (default: '~/.synapseConfig')"),
  make_option(c("-v", "--verbose"), action="store_true", default = FALSE, 
              help="Print out verbose output on script progress")
)
opt <- parse_args(OptionParser(option_list=option_list))
waitifnot(!is.null(opt$synid_file_export) && !is.null(opt$synid_file_dictionary),
          msg = "Usage: Rscript uncode_redcap_export.R -h")

synid_file_data <- opt$synid_file_export
synid_file_dd <- opt$synid_file_dictionary
synid_folder_output <- opt$synid_folder_output
file_output <- opt$file_name
auth <- opt$synapse_auth
verbose <- opt$verbose

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
                   header = header, comment.char = comment.char)
  return(data)
}

#' Read contents of an Excel Spreadsheet stored on Synapse.
#' 
#' @param synapse_id Synapse ID of the spreadsheet
#' @param version Version of the file
#' @param sheet Number of the sheet of the spreadsheet
#' @param check.names Whether R should modify names if non-conforming
#' @return Matrix of data
#' @example 
#' get_synapse_entity_data_in_xlsx(synapse_id = "syn123345", sheet = 2)
get_synapse_entity_data_in_xlsx <- function(synapse_id, 
                                            version = NA,
                                            sheet = 1,
                                            check.names = F,
                                            na.strings = c("NA")) {
  library(openxlsx)
  
  if (is.na(version)) {
    entity <- synGet(synapse_id)
  } else {
    entity <- synGet(synapse_id, version = version)
  }
  
  data <- read.xlsx(entity$path, check.names = check.names, 
                    sheet = sheet, na.strings = na.strings)
  
  return(data)
}

#' Read the REDCap data export file stored on Synapse as a '.csv' or '.xlsx' file.
#' 
#' @param synid_file_data Synapse ID of file entity with REDCap data export
#' @return Data frame containing data in the file
read_data_export <- function(synid_file_data) {
  
  data_upload <- c()
  name_upload <- synGet(synid_file_data, downloadFile = T)$path
  
  if (grepl(pattern = "\\.xlsx$", x = name_upload)) {
    data_upload <- get_synapse_entity_data_in_xlsx(synid_file_data, na.strings = c(""))
  } else if(grepl(pattern = "\\.csv$", x = name_upload)) {
    data_upload <- get_synapse_entity_data_in_xlsx(synid_file_data, na.strings = c(""))
  } else {
    print("Error: data export file must be xlxs or csv. File name must end with the appropriate suffix.  Quitting...")
    stop()
  }
  
  return(data_upload)
}

get_root_variable_name <- function(str) {
  splt <- strsplit(x = str, split = "___")
  root <- unlist(lapply(splt, FUN = function(x) {return(x[1])}))
  return(root)
}

trim_string <- function(str) {
  front <- gsub(pattern = "^[[:space:]]+", replacement = "", x = str)
  back <- gsub(pattern = "[[:space:]]+$", replacement = "", x = front)
  
  return(back)
}

merge_last_elements <- function(x, delim) {
  
  y <- c()
  y[1] = x[1]
  y[2] <- paste0(x[2:length(x)], collapse = delim)
  return(y)
}

#' Perform string split operation but only on the first
#' occurrence of the split character.
strsplit_first <- function(x, split) {
  
  unmerged <- strsplit(x = x, split = split)
  remerge <- lapply(unmerged, merge_last_elements, delim = split)
  
  return(remerge)
}

parse_mapping <- function(str) {
  
  clean <- trim_string(gsub(pattern = "\"", replacement = "", x = str))
  splt <- strsplit_first(strsplit(x = clean, split = "|", fixed = T)[[1]], split = ",")
  
  codes <- unlist(lapply(splt, FUN = function(x) {return(trim_string(x[1]))}))
  values <- unlist(lapply(splt, FUN = function(x) {return(trim_string(x[2]))}))
  mapping <- data.frame(cbind(codes, values), stringsAsFactors = F)
  
  return(mapping)
}

parse_mappings <- function(strs, labels) {
  
  mappings <- list()
  
  for (i in 1:length(strs)) {
    
    if (!is.na(strs[[i]])) {
      mappings[[labels[i]]] <- parse_mapping(strs[i])
    } else {
      mappings[[labels[i]]] <- NULL
    }
  }
  
  return(mappings)
}

uncode_data_column <- function(col_coded, mapping) {
  
  # map coded to uncoded
  col_uncoded <- data.frame(codes = col_coded) %>%
    mutate(codes_chr = as.character(codes)) %>%
    left_join(mapping, by = c("codes_chr" = "codes")) %>%
    select("values")
  
  return(col_uncoded)
}

#' Map any coded data to actual values as mapped in the 
#' REDCap Data Dictionary (DD).
#' 
#' @param data Data frame of coded data
#' @param mappings Matrix with two columns, first containing a label and
#' second columns a mapping string.  
#' @param secondary_mappings Another mapping matrix that is used secondarily
#' if the label is not found in the primary mapping matrix.
#' @return Data frame of uncoded data.
#' @example
#' map_code_to_value(data = my_data, dd = dd)
uncode_data <- function(df_coded, dd) {
  
  df_uncoded <- df_coded
  
  # gather mappings
  mappings <- parse_mappings(strs = dd[["Choices, Calculations, OR Slider Labels"]], 
                                labels = dd[["Variable / Field Name"]])

  for (i in 1:ncol(df_coded)) {
    
    var_name <- get_root_variable_name(names(df_coded)[i])
    field_type <- if (is.element(var_name, dd$`Variable / Field Name`)) dd$`Field Type`[which(dd$`Variable / Field Name`  == var_name)] else 'unknown'
    
   if (length(which(names(mappings) == var_name))) {
      
      idx_mapping <- which(names(mappings) == var_name)
      
      if (var_name == names(df_coded)[i]) {
        
        # uncode non-expanded variable according to mapping 
        df_uncoded[,i] <- uncode_data_column(col_coded = df_coded[,i], 
                                             mapping = mappings[[idx_mapping]])
      } else if (var_name != names(df_coded)[i]) {
        
        # replace 0,1 coding for expanded variables according to representative value 
        var_code <- strsplit(names(df_coded)[i], split = "___")[[1]][2]
        code_value <- mappings[[idx_mapping]]$values[which(mappings[[idx_mapping]]$codes == var_code)]
        mapping <- data.frame(codes = c("0","1"), values = c(0, code_value))
        df_uncoded[,i] <- uncode_data_column(col_coded = df_coded[,i], 
                                             mapping = mapping)
        
        # make check box coding for unchecked
        if (field_type == "checkbox") {
          
          idx_zero <- which(df_uncoded[,i] == "0")
          df_uncoded[idx_zero, i] <- NA
        }
      }
    }
  }
  
  return(df_uncoded)
}

convert_string_to_timestamp <- function(x, format = "%Y-%m-%d %H:%M") {
  
  result <- x
  
  idx_slash <- which(grepl(x = x, pattern = "/"))
  if (length(idx_slash)) {
    result[idx_slash] <- format(mdy_hm(x[idx_slash]), format)
  }
  
  idx_19 <- setdiff(which(nchar(x) == 19), idx_slash)
  if (length(idx_19)) {
    result[idx_19] <- format(ymd_hms(x[idx_19]), format)
  }
    
  idx_16 <- setdiff(which(nchar(x) == 16 | nchar(x) == 15 & grepl(x = x, pattern = "-")), idx_slash)
  if (length(idx_16)) {
    result[idx_16] <- format(ymd_hm(x[idx_16]), format)
  }
  
  return(result)
}

convert_string_to_date <- function(x, format = "%Y-%m-%d") {
  
  result <- x
  
  idx_slash <- which(grepl(x = x, pattern = "/"))
  if (length(idx_slash)) {
    result[idx_slash] <- format(mdy(x[idx_slash]), format)
  }
  
  return(result)
}

format_rca <- function(x, dd, mapping_yesno = setNames(c("No", "Yes"), c(0, 1))) {
  
  # modify date time format
  col_ts <- unlist(dd %>% filter(grepl(pattern = "datetime", x = `Text Validation Type OR Show Slider Number`)) %>% select(`Variable / Field Name`))
  if (length(col_ts)) {
    x[,col_ts] <- lapply(x[,col_ts], convert_string_to_timestamp, format = "%Y-%m-%d %H:%M")
  }
  
  # modify date format
  col_dt <- unlist(dd %>% filter(grepl(pattern = "date_mdy", x = `Text Validation Type OR Show Slider Number`)) %>% select(`Variable / Field Name`))
  if (length(col_dt)) {
    x[,col_dt] <- lapply(x[,col_dt], convert_string_to_date, format = "%Y-%m-%d")
  }
  
  # modify boolean values values
  mapping_yesno <- data.frame(codes = names(mapping_yesno),
                              values = as.character(mapping_yesno),
                              stringsAsFactors = F)
  idx_ind <- dd$`Variable / Field Name`[which(dd$`Field Type` == "yesno")]
  for (idx in idx_ind) {
    x[,idx] <- as.character(unlist(uncode_data_column(x[,idx], mapping = mapping_yesno)))
  }
  
  return(x)
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
save_to_synapse <- function(path, parent_id, file_name = NA, prov_name = NA, prov_desc = NA, prov_used = NA, prov_exec = NA) {
  
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
#' Override of synapser::synLogin() function to accept 
#' custom path to .synapseConfig file or personal authentication
#' token.  If no arguments are supplied, performs standard synLogin().
#' 
#' @param auth full path to .synapseConfig file or authentication token
#' @param silent verbosity control on login
#' @return TRUE for successful login; F otherwise
synLogin <- function(auth = NA, silent = T) {
  
  # default synLogin behavior
  if (is.na(auth)) {
    syn <- synapser::synLogin(silent = silent)
    return(T)
  }
  
  token = auth
  
  # extract token from .synapseConfig
  if (grepl(x = auth, pattern = "\\.synapseConfig$")) {
    token = get_auth_token(auth)
    
    if (is.na(token)) {
      return(F)
    }
  }
  
  # login
  syn <- tryCatch({
    synapser::synLogin(authToken = token, silent = silent)
  }, error = function(cond) {
    return(F)
  })
  
  if (is.null(syn)) {
    return(T)
  }
  return(syn)
}

# synapse login -------------------

status <- synLogin(auth = auth)

# main ----------------------------
  
if (verbose) {
  print(glue("{now(timeOnly = T)}: Reading data dictionary..."))
}

dd <- get_synapse_entity_data_in_csv(synid_file_dd, na.strings = "")

if (verbose) {
  print(glue("{now(timeOnly = T)}: Reading data uploads..."))
}

coded <- read_data_export(synid_file_data)

if (verbose) {
  print(glue("{now(timeOnly = T)}: Uncoding data uploads..."))
}

# uncode
uncoded <- uncode_data(df_coded = coded, dd = dd)

if (verbose) {
  print(glue("{now(timeOnly = T)}: Formatting uncoded data..."))
}

# format data
uncoded_formatted <- format_rca(uncoded, dd = dd)

if (verbose) {
  print(glue("{now(timeOnly = T)}: Writing uncoded data to file locally to '{getwd()}/{file_output}'..."))
}

write.csv(x = uncoded_formatted, file = file_output, row.names = F, na = "")

if (!is.na(synid_folder_output)) {
  
  if (verbose) {
    print(glue("{now(timeOnly = T)}: Saving uncoded data to Synapse..."))
  }
  
  synid_file_output <- save_to_synapse(path = file_output,
                  parent_id = synid_folder_output,
                  prov_name = "Uncoded REDCap export",
                  prov_desc = "REDCap export data with variable values uncoded according to the data dictionary",
                  prov_used = c(synid_file_data, synid_file_dd),
                  prov_exec = "https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/blob/develop/scripts/uncode_redcap_export.R")
  
  # clean up locally
  file.remove(file_output)
  
  if (verbose) {
    print(glue("{now(timeOnly = T)}: Saved file to Synapse at {synid_file_output} and removed local file. "))
  }
}

# close out ----------------------------

toc <- as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
