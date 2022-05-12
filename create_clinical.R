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
source("shared_fxns.R")

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

if (verbose) {
  print(glue("Parameters:"))
  print(glue("----------"))
  print(glue("synid_file_input:\t{synid_file_input}"))
  print(glue("synid_folder_output:\t{synid_folder_output}"))
  print(glue("synid_table_map:\t{synid_table_map}"))
  print(glue("verbose:\t\t{verbose}"))
  print(glue("----------"))
}

# functions ----------------------------

#' ETL_OPERATION function to calculate overall survival in months.
#' sampleType: PATIENT
#' cbio: OS_MONTHS
#' 
#' @param df_data data frame representing the raw dataset.
#' @return vector of integers 
get_patient_os_months <- function(df_data) {
  
  os_months <- rep(NA, nrow(df_data))
  idx_not_na <- which(!is.na(df_data$date_death_int) & !is.na(df_data$date_first_met_int))
  os_months[idx_not_na] <- round((as.double(df_data$date_death_int[idx_not_na]) - as.double(df_data$date_first_met_int[idx_not_na])) / 30.4)
  return(os_months)
}

#' ETL_OPERATION function to generate overall patient survival status.
#' sampleType: PATIENT
#' cbio: OS_STATUS
#' 
#' @param df_data data frame representing the raw dataset.
#' @return vector of strings
get_patient_os_status <- function(df_data) {
  
  os_status <- df_data$vital_status
  os_status[which(os_status =="Dead")] = "DECEASED"
  os_status[which(os_status =="Alive")] = "LIVING"
  return(os_status)
}

#' ETL_OPERATION function to query the panel sequencing assay ID for a sample.
#' sampleType: SAMPLE
#' cbio: SEQ_ASSAY_ID
#' 
#' @param df_data data frame representing the raw dataset.
#' @param sample_no integer representing numeric index of the sample.
#' @param synid_table_assay string representing Synapse ID of reference table with sample sequencing assay ID information.
#' @return vector of strings
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

#' ETL_OPERATION function to calculate sample specimen start date.
#' sampleType: TIMELINE-SPECIMEN
#' cbio: START_DATE
#' 
#' @param df_data data frame representing the raw dataset.
#' @param sample_no integer representing numeric index of the sample.
#' @return vector of integers
get_timeline_specimen_start_date <- function(df_data, sample_no) {
  
  start_date <- rep(NA, nrow(df_data))
  var_sample_id <- glue("seq_report_date_int_{sample_no}")
  
  if (length(which(colnames(df_data) == var_sample_id))) {
    idx_not_na <- which(!is.na(df_data$primary_dx_date_int) & !is.na(unlist(df_data[var_sample_id])))
    start_date[idx_not_na] <- as.double(df_data$primary_dx_date_int[idx_not_na]) - as.double(unlist(df_data[var_sample_id])[idx_not_na])
  }
  
  return(start_date)
}

#' ETL_OPERATION function to format sample IDs for GENIE for sample file.
#' sampleType: SAMPLE
#' cbio: SAMPLE_ID
#' 
#' @param df_data data frame representing the raw dataset.
#' @param sample_no integer representing numeric index of the sample.
#' @return vector of strings
get_sample_sample_id <- function(df_data, sample_no) {
  var_sample_id <- glue("sample_id_{sample_no}")
  res <- as.character(unlist(sapply(df_data[,var_sample_id], add_genie_prefix)))
  return(res)
}

#' ETL_OPERATION function to calculate format SAMPLE IDs for GENIE for timeline file.
#' sampleType: TIMELINE-SPECIMEN
#' cbio: SAMPLE_ID
#' 
#' @param df_data data frame representing the raw dataset.
#' @param sample_no integer representing numeric index of the sample.
#' @return vector of strings
get_timeline_specimen_sample_id <- function(df_data, sample_no) {
  var_sample_id <- glue("sample_id_{sample_no}")
  res <- as.character(unlist(sapply(df_data[,var_sample_id], add_genie_prefix)))
  return(res)
}

#' Map function names to implementations.
#' 
#' @param sampleType string representing cBio sample type
#' @param cbio string representing cbio variable name
#' @return object representing function
get_fxn <- function(sampleType, cbio) {
  
  map_fxn_str <- list()
  labels <- c("get_patient_os_months", 
              "get_patient_os_status", 
              "get_sample_seq_assay_id", 
              "get_timeline_specimen_start_date",
              "get_sample_sample_id",
              "get_timeline_specimen_sample_id")
  
  for (label in labels) {
    map_fxn_str[[label]] <- get(label)
  }
  
  return(map_fxn_str[[glue("get_{tolower(gsub('-', '_', sampleType))}_{tolower(cbio)}")]])
}

#' Remove rows from a data frame with a missing primary key.
#' 
#' @param df_cbio data frame 
#' @param primary_key column name on which to remove if NA
#' @return data frame with rows removed
filter_empty_rows <- function(df_cbio, primary_key) {
  idx_rm <- which(is.na(df_cbio[primary_key]))
  if (length(idx_rm)) {
    return(df_cbio[-idx_rm,])
  }
  return(df_cbio)
}

#' Remove rows from a data frame with a start date.
#' 
#' @param df_cbio data frame 
#' @param col_name column name on which to remove if NA
#' @return data frame with rows removed
filter_null_start_date <- function(df_cbio, col_start = "START_DATE") {
  idx_rm <- which(is.na(df_cbio[col_start]))
  if (length(idx_rm)) {
    return(df_cbio[-idx_rm,])
  }
  return(df_cbio)
}

#' Create data frame with patient cBioPortal data.
#' 
#' @param df_data data frame representing the raw dataset.
#' @param cbio vector of string representing cbio variable name.
#' @param code parallel vector of strings representing variable name in the raw data.
#' @param instructions instructions for cbio variable generation, if applicable.
#' @return data frame representing cBioPortal patient data
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

#' Create data frame with sample cBioPortal data.
#' 
#' @param df_data data frame representing the raw dataset.
#' @param cbio vector of string representing cbio variable name
#' @param code parallel vector of strings representing variable name in the raw data
#' @param instructions instructions for cbio variable generation, if applicable
#' @param sample_nos vector of integers representing sample indexes
#' @return data frame representing cBioPortal sample data
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

#' Create data frame with timeline-treatment cBioPortal data.
#' 
#' @param df_data data frame representing the raw dataset.
#' @param cbio vector of string representing cbio variable name
#' @param code parallel vector of strings representing variable name in the raw data
#' @param instructions instructions for cbio variable generation, if applicable
#' @param therapy_nos vector of integers representing therapy indexes
#' @param drug_nos vector of integers representing drug indexes
#' @return data frame representing cBioPortal timeline-treatment data
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

#' Create data frame with timeline-specimen cBioPortal data.
#' 
#' @param df_data data frame representing the raw dataset.
#' @param cbio vector of string representing cbio variable name
#' @param code parallel vector of strings representing variable name in the raw data
#' @param instructions instructions for cbio variable generation, if applicable
#' @param sample_nos vector of integers representing sample indexes
#' @return data frame representing cBioPortal timeline-specimen data
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

#' Create data frame with timeline-status cBioPortal data.
#' 
#' @param df_data data frame representing the raw dataset.
#' @param cbio vector of string representing cbio variable name
#' @param code parallel vector of strings representing variable name in the raw data
#' @param instructions instructions for cbio variable generation, if applicable
#' @param sample_nos vector of integers representing sample indexes
#' @return data frame representing cBioPortal timeline-status data
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

#' Create a data frame for the requested cBioPortal file dataset.
#' 
#' @param df_data data frame representing the raw dataset.
#' @param df_map data frame representing the mapping table
#' @param sampleType string representing the cBioPortal sample type
#' @return data frame representing the requested dataset
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

#' Create the custom header for cBioPortal clinical files.
#' 
#' @param df data frame representing clinical dataset
#' @param label vector of string representing a short label for each column in the dataset
#' @param description vector of string representing a long descriptions for each column in the dataset
#' @param colType vector of string representing the data type of each column in the dataset
create_cbio_clinical_header <- function(df, label, description, colType) {
  header <- rbind(label, description, colType, rep(1))
  header <- t(apply(header, 1, function(x) {return(c(paste0("#", x[1]), x[2:length(x)]))}))
  header <- rbind(header, colnames(df))
  rownames(header) <- NULL
  colnames(header) <- colnames(df)
  return(header)
}

#' Get cBioPortal clinical file name based on file type.
#' 
#' @param sampleType string representing cBioPortal sample type.  
#' @return string
get_cbio_filename <- function(sampleType) {
  mapping <- setNames(c("data_clinical_sample.txt", "data_clinical_patient.txt", "data_timeline.txt"),
                      c("SAMPLE", "PATIENT", "TIMELINE"))
  return(mapping[sampleType])
}

#' Write cBioPortal clinical file. 
#' 
#' @param df data frame representing clinical dataset.
#' @param label vector of string representing a short label for each column in the dataset.
#' @param description vector of string representing a long descriptions for each column in the dataset.
#' @param colType vector of string representing the data type of each column in the dataset.
#' @param sampleType string representing cBioPortal sample type.
#' @param delim character representing delimiter used in writing to file
#' @return file name of ouptut file
write_cbio_clinical <- function(df, label, description, colType, sampleType, delim = "\t") {
  
  filename <- get_cbio_filename(sampleType)
  header <- create_cbio_clinical_header(df, label, description, colType)
  
  df_write <- rbind(header, df)
  write.table(df_write, file = filename, sep = delim, na = "",
              col.names = F, row.names = F, quote = F)
  
  return(filename)
}

#' Write cBioPortal timeline file. 
#' 
#' @param df data frame representing clinical dataset.
#' @return file name of ouptut file
write_cbio_timeline <- function(df, delim = "\t") {
  filename <- get_cbio_filename("TIMELINE")
  write.table(df, file = filename, sep = delim, na = "",
              col.names = T, row.names = F, quote = F)
  return(filename)
}

# synapse login --------------------

status <- synLogin(auth = auth)

# read ----------------------------

if (verbose) { print(glue("{now()}: reading data and mapping...")) }

df_data <- get_synapse_entity_data_in_csv(synid_file_input, na.strings = c("NA",""))

query <- glue("SELECT * FROM {synid_table_map} ORDER BY sampleType, columnOrder")
df_map <- as.data.frame(synTableQuery(query, includeRowIdAndRowVersion = F))

# main ----------------------------

if (verbose) { print(glue("{now()}: formatting data for cBioPortal files...")) }

# create each data frame
dfs <- list()
for (sampleType in unlist(df_map %>% select(sampleType) %>% distinct())) {
  if (verbose) { print(glue("{now()}: generating {sampleType} data frame...")) }
  dfs[[sampleType]] <- create_df(df_data = df_data, df_map = df_map, sampleType = sampleType)
}

# write -----------------

if (verbose) { print(glue("{now()}: writing files locally...")) }

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
  
  if (verbose) { print(glue("{now()}: storing files to Synapse in folder {synid_folder_output}...")) }
  
  for (outfile in outfiles) {
    synid_file_cbio <- save_to_synapse(path = outfile, 
                                      parent_id = synid_folder_output, 
                                      prov_name = "ERBB2 cBioPortal file", 
                                      prov_desc = "GENIE ERBB2 Sponsored Project data in cBioPortal format", 
                                      prov_used = c(synid_file_input, synid_table_map), 
                                      prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/main/create_clinical.R")
    file.remove(outfile)
  }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
