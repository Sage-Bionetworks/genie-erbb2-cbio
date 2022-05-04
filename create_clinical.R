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

# user input ----------------------------

option_list <- list( 
  make_option(c("-i", "--synid_file_input"), type = "character",
              help="Synapse ID of uncoded REDCap export"),
  make_option(c("-m", "--synid_table_map"), type = "character",
              help="Synapse ID of table of REDCap to cBio mapping"),
  make_option(c("-o", "--synid_folder_output"), type = "character", default = NA,
              help="Synapse ID of output folder (default: write locally only)"),
  make_option(c("-v", "--verbose"), action="store_true", default = F, 
              help="Output script messages to the user. (default: FALSE)"),
  make_option(c("-a", "--auth"), type = "character", default = NA,
              help="Synapse personal access token or path to .synapseConfig (default: normal synapse login behavior)")
)
opt <- parse_args(OptionParser(option_list=option_list))
waitifnot(!is.null(opt$synid_file_input) && !is.null(opt$synid_table_map),
          msg = "Rscript create_clinical.R -h")

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
  
  os_months <- rep(NA, nrow(df_data))
  idx_not_na <- which(!is.na(df_data$date_death_int) & !is.na(df_data$date_first_met_int))
  os_months[idx_not_na] <- round((as.double(df_data$date_death_int[idx_not_na]) - as.double(df_data$date_first_met_int[idx_not_na])) / 30.4)
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
  sample_ids <- paste0(df_data[,var_sample_id], collapse = "','")
  query <- glue("SELECT SAMPLE_ID, SEQ_ASSAY_ID FROM {synid_table_assay} WHERE SAMPLE_ID IN ('{sample_ids}')")
  df_assay <- as.data.frame(synTableQuery(query, includeRowIdAndRowVersion = F))
  
  df_data_assay <- df_data %>% 
    left_join(df_assay, by = setNames("SAMPLE_ID", var_sample_id)) %>%
    select(SEQ_ASSAY_ID)
  
  return(unlist(df_data_assay))
}

get_timeline_specimen_start_date <- function(df_data, sample_no) {
  
  start_date <- rep(NA, nrow(df_data))
  var_sample_id <- glue("seq_report_date_int_{sample_no}")
  
  if (length(which(colnames(df_data) == var_sample_id))) {
    idx_not_na <- which(!is.na(df_data$primary_dx_date_int) & !is.na(unlist(df_data[var_sample_id])))
    start_date[idx_not_na] <- as.double(df_data$primary_dx_date_int[idx_not_na]) - as.double(unlist(df_data[var_sample_id])[idx_not_na])
  }
  
  return(start_date)
}

get_fxn <- function(sampleType, cbio) {
  
  map_fxn_str <- list()
  labels <- c("get_patient_os_months", "get_patient_os_status", 
              "get_sample_seq_assay_id", "get_timeline_specimen_start_date")
  
  for (label in labels) {
    map_fxn_str[[label]] <- get(label)
  }
  
  return(map_fxn_str[[glue("get_{tolower(gsub('-', '_', sampleType))}_{tolower(cbio)}")]])
}

filter_empty_rows <- function(df_cbio, primary_key) {
  idx_rm <- which(is.na(df_cbio[primary_key]))
  if (length(idx_rm)) {
    return(df_cbio[-idx_rm,])
  }
  return(df_cbio)
}

filter_null_start_date <- function(df_cbio, col_start = "START_DATE") {
  idx_rm <- which(is.na(df_cbio[col_start]))
  if (length(idx_rm)) {
    return(df_cbio[-idx_rm,])
  }
  return(df_cbio)
}

create_patient <- function(df_data, cbio, code, instructions = NULL) {
  
  mat <- matrix(NA, nrow = nrow(df_data), ncol = length(cbio), dimnames = list(c(), cbio))
  
  for (i in 1:length(cbio)) {
    if (code[i] == "ETL_OPERATION") {
      fxn <- get_fxn(sampleType = "PATIENT", cbio = cbio[i])
      mat[,cbio[i]] <- fxn(df_data)
    } else {
      mat[,cbio[i]] <- df_data[,code[i]]
    }
  }
  
  return(data.frame(mat))
}

create_sample <- function(df_data, cbio, code, instructions = NULL, sample_nos = c(1:10)) {
  
  mat <- matrix(nrow = 0, ncol = length(cbio), dimnames = list(c(), cbio))
  
  for (sample_no in sample_nos) {
    
    mat_no <- matrix(NA, nrow = nrow(df_data), ncol = length(cbio), dimnames = list(c(), cbio))
    
    for (i in 1:length(cbio)) {
      code_no <- glue(code[i])
      if (code_no == "ETL_OPERATION") {
        fxn <- get_fxn(sampleType = "SAMPLE", cbio = cbio[i])
        mat_no[,cbio[i]] <- fxn(df_data, sample_no)
      } else {
        if (length(which(colnames(df_data) == code_no))) {
          mat_no[,cbio[i]] <- df_data[,code_no]
        }
      }
    }
    
    mat <- rbind(mat, mat_no)
  }
  
  df_final <- filter_empty_rows(df_cbio = data.frame(mat), primary_key = "SAMPLE_ID") 
  
  return(df_final)
}

create_timeline_treatment <- function(df_data, cbio, code, instructions, 
                                      therapy_nos = c(1:25), drug_nos = c(1:5)) {
  
  mat <- matrix(nrow = 0, ncol = length(cbio), dimnames = list(c(), cbio))
  
  for (therapy_no in therapy_nos) {
    for (drug_no in drug_nos) {
      
      mat_no <- matrix(nrow = nrow(df_data), ncol = length(cbio), dimnames = list(c(), cbio))
      
      for (i in 1:length(cbio)) {
        
        code_no <- glue(code[i])
        instructions_no <- glue(instructions[i])
        
        if (code_no == "ETL_CONSTANT") {
          mat_no[,cbio[i]] <- instructions_no
        } else {
          if (length(which(colnames(df_data) == code_no))) {
            mat_no[,cbio[i]] <- df_data[,code_no]
          }
        }
      }
      
      mat <- rbind(mat, mat_no)
    }
  }
  
  df_final <- filter_empty_rows(df_cbio = data.frame(mat), primary_key = "AGENT") 
  df_final <- filter_null_start_date(df_cbio = df_final)
  
  return(df_final)
}

create_timeline_specimen <- function(df_data, cbio, code, instructions, sample_nos = c(1:10)) {
  mat <- matrix(nrow = 0, ncol = length(cbio), dimnames = list(c(), cbio))
  
  for (sample_no in sample_nos) {
    
    mat_no <- matrix(NA, nrow = nrow(df_data), ncol = length(cbio), dimnames = list(c(), cbio))
    
    for (i in 1:length(cbio)) {
      code_no <- glue(code[i])
      if (code_no == "ETL_OPERATION") {
        fxn <- get_fxn(sampleType = "TIMELINE-SPECIMEN", cbio = cbio[i])
        mat_no[,cbio[i]] <- fxn(df_data, sample_no)
      } else if (code_no == "ETL_CONSTANT") {
        mat_no[,cbio[i]] <- instructions[i]
      } else {
        if (length(which(colnames(df_data) == code_no))) {
          mat_no[,cbio[i]] <- df_data[,code_no]
        }
      }
    }
    
    mat <- rbind(mat, mat_no)
  }
  
  df_final <- filter_null_start_date(df_cbio = data.frame(mat))
  
  return(df_final)
}

create_timeline_status <- function(df_data, cbio, code, instructions) {
  mat <- matrix(NA, nrow = nrow(df_data), ncol = length(cbio), dimnames = list(c(), cbio))
  
  for (i in 1:length(cbio)) {
    if (code[i] == "ETL_CONSTANT") {
      mat[,cbio[i]] <- instructions[i]
    } else {
      mat[,cbio[i]] <- df_data[,code[i]]
    }
  }
  
  df_final <- filter_null_start_date(df_cbio = data.frame(mat))
  
  return(df_final)
}

create_df <- function(df_data, df_map, sampleType) {
  
  df_cbio_type <- NULL
  
  df_map_type <- df_map %>% filter(sampleType == (!!sampleType))
  cbio = df_map_type$cbio
  code = df_map_type$code
  instructions = df_map_type$instructions
  
  if (sampleType == "PATIENT") {
    df_cbio_type <- create_patient(df_data, cbio, code, instructions)
  } else if (sampleType == "SAMPLE") {
    df_cbio_type <- create_sample(df_data, cbio, code, instructions)
  } else if (sampleType == "TIMELINE-TREATMENT") {
    df_cbio_type <- create_timeline_treatment(df_data, cbio, code, instructions)
  } else if (sampleType == "TIMELINE-SPECIMEN") {
    df_cbio_type <- create_timeline_specimen(df_data, cbio, code, instructions)
  } else if (sampleType == "TIMELINE-STATUS") {
    df_cbio_type <- create_timeline_status(df_data, cbio, code, instructions)
  } else {
    message(glue("ERROR: {sampleType} not implemented.  Quitting..."))
    stop()
  }
  
  return(df_cbio_type)
}

create_cbio_clinical_header <- function(df, label, description, colType) {
  header <- rbind(label, description, colType, rep(1))
  header <- t(apply(header, 1, function(x) {return(c(paste0("#", x[1]), x[2:length(x)]))}))
  header <- rbind(header, colnames(df))
  rownames(header) <- NULL
  colnames(header) <- colnames(df)
  return(header)
}

get_cbio_filename <- function(file_type) {
  mapping <- setNames(c("data_clinical_sample.txt", "data_clinical_patient.txt", "data_timeline.txt"),
                      c("SAMPLE", "PATIENT", "TIMELINE"))
  return(mapping[file_type])
}

write_cbio_clinical <- function(df, label, description, colType, sampleType, delim = "\t") {
  
  filename <- get_cbio_filename(sampleType)
  header <- create_cbio_clinical_header(df, label, description, colType)
  
  df_write <- rbind(header, df)
  write.table(df_write, file = filename, sep = delim, na = "",
              col.names = F, row.names = F, quote = F)
  
  return(filename)
}

write_cbio_timeline <- function(df, delim = "\t") {
  filename <- get_cbio_filename("TIMELINE")
  write.table(df, file = filename, sep = delim, na = "",
              col.names = T, row.names = F, quote = F)
  return(filename)
}

# synapse login --------------------

status <- synLogin(auth = auth)

# read ----------------------------

if (verbose) {
  print(glue("{now(timeOnly=T)}: reading data and mapping..."))
}

df_data <- get_synapse_entity_data_in_csv(synid_file_input, na.strings = c("NA",""))

query <- glue("SELECT * FROM {synid_table_map} ORDER BY sampleType, columnOrder")
df_map <- as.data.frame(synTableQuery(query, includeRowIdAndRowVersion = F))

# main ----------------------------

if (verbose) {
  print(glue("{now(timeOnly=T)}: formatting data for cBioPortal files..."))
}

# create each data frame
dfs <- list()
for (sampleType in unlist(df_map %>% select(sampleType) %>% distinct())) {
  if (sampleType != "ETL") {
    if (verbose) {
      print(glue("{now(timeOnly=T)}: generating {sampleType} data frame..."))
    }
    dfs[[sampleType]] <- create_df(df_data = df_data, df_map = df_map, sampleType = sampleType)
  }
}

# write -----------------

if (verbose) {
  print(glue("{now(timeOnly=T)}: writing files locally..."))
}

# write clinical data frames to file locally
outfiles <- c()
for (sampleType in c("PATIENT", "SAMPLE")) {
  outfiles[sampleType] <- write_cbio_clinical (dfs[[sampleType]], 
                                      label = unlist(df_map %>% filter(sampleType == (!!sampleType)) %>% select(labels)),
                                      description = unlist(df_map %>% filter(sampleType == (!!sampleType)) %>% select(description)),
                                      sampleType = sampleType,
                                      colType = unlist(df_map %>% filter(sampleType == (!!sampleType)) %>% select(colType)))
}

# consolidate write timeline data frame locally
df_timeline <- NULL
for (sampleType in grep(names(dfs), pattern = "^TIMELINE-", value = T)) {
  if (is.null(df_timeline)) {
    df_timeline <- dfs[[sampleType]]
  } else {
    df_timeline <- df_timeline %>% bind_rows(dfs[[sampleType]])
  }
}
outfiles["TIMELINE"] <- write_cbio_timeline(df_timeline)

# store to synpase
if (!is.na(synid_folder_output)) {
  
  if (verbose) {
    print(glue("{now(timeOnly=T)}: storing files to Synapse in folder {synid_folder_output}..."))
  }
  
  for (outfile in outfiles) {
    synid_file_cbio <- save_to_synapse(path = outfile, 
                                      parent_id = synid_folder_output, 
                                      prov_name = "ERBB2 cBioPortal file", 
                                      prov_desc = "GENIE ERBB2 Sponsored Project data in cBioPortal format", 
                                      prov_used = c(synid_file_input, synid_table_map), 
                                      prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/main/create_clinical.R")
  }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
