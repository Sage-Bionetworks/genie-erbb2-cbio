# Description: Create case files for ERBB2 cBioPortal files. 
# Format: 
# - https://docs.cbioportal.org/file-formats/#case-lists
# Author: Haley Hunter-Zinck
# Date: 2022-05-05

# setup ----------------------------

tic = as.double(Sys.time())

library(glue)
library(dplyr)
library(synapser)
library(optparse)
source("shared_fxns.R")

# constants
cancer_study_identifier <- "genie_erbb2"
filenames <- c("case-breast" = "cases_Breast_Cancer.txt", 
               "case-all" = "cases_all.txt", 
               "case-cna" = "cases_cna.txt", 
               "case-cnaseq" = "cases_cnaseq.txt", 
               "case-seq" = "cases_sequenced.txt")

# user input ----------------------------

option_list <- list( 
  make_option(c("-i", "--synid_folder_cbio"), type = "character",
              help="Synapse ID of folder containing cBioPortal clinical and genomic files"),
  make_option(c("-o", "--synid_folder_output"), type = "character", default = NA,
              help="Synapse ID of output folder (default: write locally only)"),
  make_option(c("-v", "--verbose"), action="store_true", default = F, 
              help="Output script messages to the user. (default: FALSE)"),
  make_option(c("-a", "--auth"), type = "character", default = NA,
              help="Synapse personal access token or path to .synapseConfig (default: normal synapse login behavior)")
)
opt <- parse_args(OptionParser(option_list=option_list))
waitifnot(!is.null(opt$synid_folder_cbio),
          msg = "Rscript create_case.R -h")

synid_folder_cbio <- opt$synid_folder_cbio
synid_folder_output <- opt$synid_folder_output
verbose <- opt$verbose
auth <- opt$auth

# functions ----------------------------

#' Create a vector representing each row of a generic cBioPortal case file given the values for all fields.
#' 
#' @param cancer_study_identifier study label
#' @param stable_id cancer_study_identifier followed by underscore and custom suffix
#' @param case_list_description description of the patients included in the case list
#' @param case_list_ids vector of strings representing sample IDs
#' @return vector of strings representing each row of the case file.
create_case_generic <- function(cancer_study_identifier, stable_id, case_list_name, case_list_description, case_list_ids) {
  
  rows <- rep("", 5)
  rows[1] <- c(glue("cancer_study_identifier: {cancer_study_identifier}"))
  rows[2] <- c(glue("stable_id: {stable_id}"))
  rows[3] <- c(glue("case_list_name: {case_list_name}"))
  rows[4] <- c(glue("case_list_description: {case_list_description}"))
  rows[5] <- c(glue("case_list_ids: {paste0(case_list_ids, collapse = '\t')}"))
  
  return(rows)
}

#' Create case list for all samples in the study.
#' 
#' @param synid_files_cbio named vector of strings representing Synapse IDs of cBioPortal files; names represent file labels.
#' @param cancer_study_identifier study label
#' @return vector of strings representing each row of the case file.
create_case_all <- function(synid_files_cbio, cancer_study_identifier) {
  
  synid_file_sample <- as.character(synid_files_cbio["sample"])
  df_sample <- get_synapse_entity_data_in_csv(synid_file_sample, sep = "\t")
  
  df_case <- create_case_generic(cancer_study_identifier=cancer_study_identifier, 
                                 stable_id = glue("{cancer_study_identifier}_all"), 
                                 case_list_name = "All samples", 
                                 case_list_description = "All samples", 
                                 case_list_ids = df_sample$SAMPLE_ID)
  return(df_case)
}

#' Create case list for all samples in the study with breast cancer.
#' 
#' @param synid_files_cbio named vector of strings representing Synapse IDs of cBioPortal files; names represent file labels.
#' @param cancer_study_identifier study label
#' @return vector of strings representing each row of the case file.
create_case_breast <- function(synid_files_cbio, cancer_study_identifier) {
  synid_file_sample <- as.character(synid_files_cbio["sample"])
  df_sample <- get_synapse_entity_data_in_csv(synid_file_sample, sep = "\t")
  
  sample_ids <- unlist(df_sample %>% 
    filter(grepl(pattern = "Breast" , x = ONCOTREE_CODE)) %>%
    select(SAMPLE_ID))
  
  df_case <- create_case_generic(cancer_study_identifier=cancer_study_identifier, 
                                 stable_id = glue("{cancer_study_identifier}_Breast_Cancer"), 
                                 case_list_name = "Tumor Type: Breast Cancer", 
                                 case_list_description = "All tumors with cancer type Breast Cancer", 
                                 case_list_ids = df_sample$SAMPLE_ID)
  return(df_case)
}

#' Create case list for all samples in the study with CNA data.
#' 
#' @param synid_files_cbio named vector of strings representing Synapse IDs of cBioPortal files; names represent file labels.
#' @param cancer_study_identifier study label
#' @return vector of strings representing each row of the case file.
create_case_cna <- function(synid_files_cbio, cancer_study_identifier) {

  synid_file_matrix<- as.character(synid_files_cbio["matrix"])
  df_matrix <- get_synapse_entity_data_in_csv(synid_file_matrix, sep = "\t", na.strings = "NA")
  
  sample_ids <- as.character(unlist(df_matrix %>% 
    filter(!is.na(cna)) %>%
    select(SAMPLE_ID)))
  
  df_case <- create_case_generic(cancer_study_identifier=cancer_study_identifier, 
                                 stable_id = glue("{cancer_study_identifier}_cna"), 
                                 case_list_name = "Samples with CNA", 
                                 case_list_description = "Samples with CNA", 
                                 case_list_ids = sample_ids)
  return(df_case)
}

#' Create case list for all samples in the study with mutation and CNA data.
#' 
#' @param synid_files_cbio named vector of strings representing Synapse IDs of cBioPortal files; names represent file labels.
#' @param cancer_study_identifier study label
#' @return vector of strings representing each row of the case file.
create_case_cnaseq <- function(synid_files_cbio, cancer_study_identifier) {
  
  synid_file_matrix<- as.character(synid_files_cbio["matrix"])
  df_matrix <- get_synapse_entity_data_in_csv(synid_file_matrix, sep = "\t", na.strings = "NA")
  
  sample_ids <- as.character(unlist(df_matrix %>% 
                                      filter(!is.na(cna) & !is.na(mutations)) %>%
                                      select(SAMPLE_ID)))
  
  df_case <- create_case_generic(cancer_study_identifier=cancer_study_identifier, 
                                 stable_id = glue("{cancer_study_identifier}_cnaseq"), 
                                 case_list_name = "Samples with CNA and mutation", 
                                 case_list_description = "Samples with CNA and mutation", 
                                 case_list_ids = sample_ids)
  return(df_case)
}

#' Create case list for all samples in the study with mutation or CNA data.
#' 
#' @param synid_files_cbio named vector of strings representing Synapse IDs of cBioPortal files; names represent file labels.
#' @param cancer_study_identifier study label
#' @return vector of strings representing each row of the case file.
create_case_seq <- function(synid_files_cbio, cancer_study_identifier) {
  synid_file_matrix<- as.character(synid_files_cbio["matrix"])
  df_matrix <- get_synapse_entity_data_in_csv(synid_file_matrix, sep = "\t", na.strings = "NA")
  
  sample_ids <- as.character(unlist(df_matrix %>% 
                                      filter(!is.na(cna) | !is.na(mutations)) %>%
                                      select(SAMPLE_ID)))
  
  df_case <- create_case_generic(cancer_study_identifier=cancer_study_identifier, 
                                 stable_id = glue("{cancer_study_identifier}_sequenced"), 
                                 case_list_name = "Sequenced Tumors", 
                                 case_list_description = "All sequenced samples", 
                                 case_list_ids = sample_ids)
}

#' Create the requested case list file.
#' 
#' @param file_type label identifying a particular case list.
#' @param synid_files_cbio named vector of strings representing Synapse IDs of cBioPortal files; names represent file labels.
#' @param cancer_study_identifier study label
#' @return vector of strings representing each row of the case file.
create_case <- function(file_type, synid_files_cbio, cancer_study_identifier) {
  
  df_type <- NULL
  
  if (file_type == "case-breast") {
    df_type <- create_case_breast(synid_files_cbio = synid_files_cbio, 
                                  cancer_study_identifier = cancer_study_identifier)
  } else if (file_type == "case-all") {
    df_type <- create_case_all(synid_files_cbio = synid_files_cbio, 
                               cancer_study_identifier = cancer_study_identifier)
  } else if (file_type == "case-cna") {
    df_type <- create_case_cna(synid_files_cbio = synid_files_cbio, 
                               cancer_study_identifier = cancer_study_identifier)
  } else if (file_type == "case-cnaseq") {
    df_type <- create_case_cnaseq(synid_files_cbio = synid_files_cbio, 
                                  cancer_study_identifier = cancer_study_identifier)
  } else if (file_type == "case-seq") {
    df_type <- create_case_seq(synid_files_cbio = synid_files_cbio, 
                               cancer_study_identifier = cancer_study_identifier)
  } else {
    message(glue("Warning: file type {file_type} not implemented.  Returnning NULL."))
  }
  
  return(df_type)
}

#' Write the vector of strings representing the rows of a case list file to a local file.
#' 
#' @param df_file data to write.
#' @param filename name of the output file.
#' @return string representing name of output file.
write_case <- function(df_file, filename) {
  write.table(df_file, row.names = F, file = filename, na = "", quote = F, col.names = F)
  return(filename)
}

# synapse login --------------------

status <- synLogin(auth = auth)

# read ----------------------------

if (verbose) { print(glue("{now()}: gathering cBioPortal files from '{get_synapse_entity_name(synid_folder_cbio)}' ({synid_folder_cbio})...")) }
synid_files_cbio <- get_cbio_release_files(synid_folder_cbio)

# main ----------------------------

for (i in 1:length(filenames)) {
  
  file_type <- names(filenames)[i]
  filename <- as.character(filenames[i]) 
  
  if (verbose) { glue("{now()}: creating '{file_type}' case file...") }
  df_file <- create_case(file_type = file_type, synid_files_cbio = synid_files_cbio, cancer_study_identifier = cancer_study_identifier)
  
  if (verbose) { glue("{now()}: writing '{file_type}' to file '{filename}'...") }
  filename <- write_case(df_file, file = filename)
  
  if (!is.na(synid_folder_output)) {
    
    # create subfolder to store case_lists
    synid_folder_case <- create_synapse_folder(name = "case_lists", parent_id = synid_folder_output)
    
    if (verbose) { print(glue("{now()}: saving '{filename}' to Synapse folder '{get_synapse_entity_name(synid_folder_case)}' ({synid_folder_case})...")) }
    
    synid_file_df <- save_to_synapse(path = filename, 
                                     parent_id = synid_folder_case, 
                                     prov_name = "ERBB2 cBioPortal file", 
                                     prov_desc = "GENIE ERBB2 Sponsored Project data in cBioPortal format", 
                                     prov_used = as.character(synid_files_cbio), 
                                     prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/main/create_cases.R")
    
    if (verbose) { print(glue("{now()}: file '{filename}' saved to {synid_file_df}.")) } 
    
    file.remove(filename)
  }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
