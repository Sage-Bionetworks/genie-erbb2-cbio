# Description: Create cBioPortal genomic information files using a sample ID list
#     and a main GENIE release from which to draw data. 
# Format: 
# - https://docs.cbioportal.org/file-formats/#mutation-data
# - https://docs.cbioportal.org/file-formats/#discrete-copy-number-data
# - https://docs.cbioportal.org/file-formats/#gene-panel-matrix-file
# - https://docs.cbioportal.org/file-formats/#fusion-data
# - https://docs.cbioportal.org/file-formats/#segmented-data
# Author: Haley Hunter-Zinck
# Date: 2022-05-03

# setup ----------------------------

tic = as.double(Sys.time())

library(glue)
library(dplyr)
library(sqldf)
library(synapser)
library(optparse)
source("shared_fxns.R")

# constants
sample_nos <- c(1:10)
filenames <- c("maf" = "data_mutations_extended.txt", 
               "cna" = "data_CNA.txt", 
               "matrix" = "data_gene_matrix.txt", 
               "fusion" = "data_fusions.txt", 
               "seg" = "genie_erbb2_data_cna_hg19.seg")

# user input ----------------------------

option_list <- list( 
  make_option(c("-i", "--synid_file_input"), type = "character",
              help="Synapse ID of uncoded REDCap export"),
  make_option(c("-g", "--synid_folder_mg"), type = "character", default = "syn26706564",
              help="Synapse ID of main GENIE release from which to draw genetic data (default: syn26706564)"),
  make_option(c("-o", "--synid_folder_output"), type = "character", default = NA,
              help="Synapse ID of output folder (default: write locally only)"),
  make_option(c("-v", "--verbose"), action="store_true", default = F, 
              help="Output script messages to the user. (default: FALSE)"),
  make_option(c("-a", "--auth"), type = "character", default = NA,
              help="Synapse personal access token or path to .synapseConfig (default: normal synapse login behavior)")
)
opt <- parse_args(OptionParser(option_list=option_list))
waitifnot(!is.null(opt$synid_file_input),
          msg = "Rscript create_genomic.R -h")

synid_file_input <- opt$synid_file_input
synid_folder_output <- opt$synid_folder_output
synid_folder_mg <- opt$synid_folder_mg
verbose <- opt$verbose
auth <- opt$auth

if (verbose) {
  print(glue("Parameters:"))
  print(glue("----------"))
  print(glue("synid_file_input:\t{synid_file_input}"))
  print(glue("synid_folder_output:\t{synid_folder_output}"))
  print(glue("synid_folder_mg:\t{synid_folder_mg}"))
  print(glue("verbose:\t\t{verbose}"))
}

# functions ----------------------------


#' Approximately determine the number of commented header rows. For used with read.csv.sql. 
#' 
#' @param filepath Path to file to examine.
#' @param n_check Number of lines to check in the file.
#' @return number of commented lines in checked lines
get_number_commented_header <- function(filepath, n_check = 100, comment = "#") {
  first_rows <- scan(filepath, what = "charcter", nlines = n_check, sep = "\n", quiet = T)
  n_skip = length(which(grepl(pattern = glue("^{comment}"), x = first_rows)))
  return(n_skip)
}

#' Extract and standardize GENIE sample IDs for downstream processing.  
#' 
#' @param df_data original data frame containing all data, including sample IDs
#' @return vector of GENIE sample IDs
get_sample_ids <- function(df_data, sample_nos) {
  
  var_sample <- paste0("sample_id_", sample_nos)
  raw_ids <- as.character(unlist(sapply(var_sample, 
                           function(x, df) {df[which(!is.na(df[x]) & df[x] != "Not available"), x]}, df = df_data)))
  trimmed_ids <- trim(raw_ids)
  unique_ids <- unique(trimmed_ids)
  genie_ids <- as.character(unlist(sapply(unique_ids, add_genie_prefix)))
  
  return(genie_ids)
}

#' Create mutation data cBioPortal file in maf extended format.  
#' 
#' @param synid_file_mg string representing Synapse ID of main GENIE sample file
#' @param sample_ids vector of string representing sample IDs with mutation data
#' @return data frame representing mutation data in cBioPortal format
create_maf <- function(synid_file_mg, sample_ids) {
  
  filepath <- synGet(synid_file_mg)$path
  n_skip = 0
  
  # check for commented headers (occurs in older releases)
  n_skip <- get_number_commented_header(filepath)
  
  # subset maf
  str_ids <- paste0(sample_ids, collapse = "','")
  query <- glue("SELECT * FROM file WHERE Tumor_Sample_Barcode IN ('{str_ids}')")
  df_subset <- read.csv.sql(filepath, sql = query, sep = "\t", skip = n_skip, colClasses = "character")
  
  # custom formatting for certain columns
  for (colname in grep(pattern = "^[nt]_", x = colnames(df_subset), value = T)) {
    df_subset[which(df_subset[,colname] == "."),colname] <- ""
  }
  df_subset[, "Validation_Status"] <- ""
    
  return(df_subset)
}

#' Create CNA data cBioPortal file.  
#' 
#' @param synid_file_mg string representing Synapse ID of main GENIE CNA file
#' @param sample_ids vector of string representing sample IDs with mutation data
#' @return data frame representing CNA data in cBioPortal format
create_cna <- function(synid_file_mg, sample_ids, n_row_per_chunk = 50) {
  
  df_subset <- c()
  filepath <- synGet(synid_file_mg)$path
  
  # determine samples in file
  header_cna <- scan(filepath, nlines = 1, what = "character", quiet = T)
  both_sample_ids <- intersect(header_cna, sample_ids)
  idx <- c(which(header_cna == "Hugo_Symbol"), which(is.element(header_cna, sample_ids)))
  
  # subset cna in chunks
  n_row_total <- as.double(strsplit(system(glue("wc -l {filepath}"), intern = T), split = " ")[[1]][1])
  df_subset <- matrix(NA, nrow = n_row_total - 1, ncol = length(idx), dimnames = list(c(), header_cna[idx]))
  for (i in seq(from = 2, to = n_row_total, by = n_row_per_chunk)) {
    chunk <- matrix(scan(filepath, nlines = n_row_per_chunk, skip = i-1, what = "character", quiet = T), ncol = length(header_cna), byrow = T)[,idx]
    df_subset[(i-1):(i+nrow(chunk)-2),] <- chunk
  }
  
  return(df_subset)
}

#' Either subset the existing cBioportal matrix file or create from the maf and CNA files.
#' 
#' @param synid_file_mg Synapse ID of main GENIE matrix file.
#' @param sample_ids Vector of sample IDs for which to filter
#' @param synid_file_maf Synapse ID of main GENIE maf file.
#' @param synid_file_cna Synapse ID of main GENIE CNA file.
#' @param synid_file_sam Synapse ID of main GENIE clinical sample file.
#' @return data frame representing data for matrix file for sample IDs. 
create_matrix <- function(synid_file_mg, sample_ids, 
                          synid_file_maf = NA, 
                          synid_file_cna = NA,
                          synid_file_sam = NA) {
  
  df_subset <- NULL
  
  # create new matrix file or subset existing matrix file
  if (is.na(synid_file_mg)) {
    df_mg_maf <- create_maf(synid_file_maf, sample_ids)
    df_mg_cna <- create_cna(synid_file_cna, sample_ids)
    df_mg_sam <- get_synapse_entity_data_in_csv(synid_file_sam, sep = "\t")
    
    maf_assay <- df_mg_maf %>% 
      filter(is.element(Tumor_Sample_Barcode, sample_ids)) %>%
      rename(SAMPLE_ID = Tumor_Sample_Barcode) %>%
      select(SAMPLE_ID) %>%
      distinct() %>%
      left_join(df_mg_sam, by = "SAMPLE_ID") %>%
      rename(mutations = SEQ_ASSAY_ID) %>%
      select(SAMPLE_ID, mutations)
    cna_assay <- data.frame(SAMPLE_ID = colnames(df_mg_cna)[2:ncol(df_mg_cna)]) %>% 
      select(SAMPLE_ID) %>%
      distinct() %>%
      left_join(df_mg_sam, by = "SAMPLE_ID") %>%
      rename(cna = SEQ_ASSAY_ID) %>%
      select(SAMPLE_ID, cna)
      
    df_subset <- data.frame(SAMPLE_ID = sample_ids) %>%
      left_join(maf_assay, by = "SAMPLE_ID") %>%
      left_join(cna_assay, by = "SAMPLE_ID") %>%
      filter(!is.na(mutations) | !is.na(cna))
  } else {
    df_mg <- get_synapse_entity_data_in_csv(synid_file_mg, sep = "\t", na.strings = "")
    df_subset <- df_mg %>% 
      filter(is.element(SAMPLE_ID, sample_ids))
  }
  
  return(df_subset)
}

#' Create fusion data cBioPortal file.  
#' 
#' @param synid_file_mg string representing Synapse ID of main GENIE fusion file
#' @param sample_ids vector of string representing sample IDs with mutation data
#' @return data frame representing fusion data in cBioPortal format
create_fusion <- function(synid_file_mg, sample_ids) {
  df_mg <- get_synapse_entity_data_in_csv(synid_file_mg, sep = "\t", na.strings = c("NA", ""))
  df_subset <- df_mg %>% 
    filter(is.element(Tumor_Sample_Barcode, sample_ids) & (!is.na(Hugo_Symbol) | !is.na(Entrez_Gene_Id)))
  return(df_subset)
}

#' Create seg data cBioPortal file.  
#' 
#' @param synid_file_mg string representing Synapse ID of main GENIE seg file
#' @param sample_ids vector of string representing sample IDs with mutation data
#' @return data frame representing seg data in cBioPortal format
create_seg <- function(synid_file_mg, sample_ids) {
  filepath <- synGet(synid_file_mg)$path
  
  # subset maf
  str_ids <- paste0(sample_ids, collapse = "','")
  query <- glue("SELECT * FROM file WHERE ID IN ('{str_ids}')")
  df_subset <- read.csv.sql(filepath, 
                            sql = query, 
                            sep = "\t",
                            colClasses = c(rep("character",2), rep("integer", 2), rep("character",2)))
  
  return(df_subset)
}

#' Create data frame for requested genomic cBioPortal file.
#' 
#' @param file_type string representing cBioportal file type label.  
#' @param synid_files_mg named vector of strings representing Synapse IDs of main GENIE cBioPortal files; names represent file labels.
#' @param sample_ids vector of string representing sample IDs in the study.
#' @return data frame
create_df <- function(file_type, synid_files_mg, sample_ids) {
  
  df_type <- NULL
  
  if (file_type == "maf") {
    df_type <- create_maf(synid_file_mg = as.character(synid_files_mg[file_type]), 
                          sample_ids = sample_ids)
  } else if (file_type == "cna") {
    df_type <- create_cna(synid_file_mg = as.character(synid_files_mg[file_type]), 
                          sample_ids = sample_ids)
  } else if (file_type == "matrix") {
    df_type <- create_matrix(synid_file_mg = as.character(synid_files_mg[file_type]), 
                             sample_ids = sample_ids, 
                             synid_file_maf = as.character(synid_files_mg["maf"]),
                             synid_file_cna = as.character(synid_files_mg["cna"]),
                             synid_file_sam = as.character(synid_files_mg["sample"]))
  } else if (file_type == "fusion") {
    df_type <- create_fusion(synid_file_mg = as.character(synid_files_mg[file_type]), 
                             sample_ids = sample_ids)
  } else if (file_type == "seg") {
    df_type <- create_seg(synid_file_mg = as.character(synid_files_mg[file_type]), 
                          sample_ids = sample_ids)
  } else if (file_type == "panel") {
    df_type <- create_panel(synid_file_mg = as.character(synid_files_mg[file_type]), sample_ids = sample_ids)
  } else {
    message(glue("Warning: file type {file_type} not implemented.  Returnning NULL."))
  }
  
  return(df_type)
}

#' Write data frame to a local file.
#' 
#' @param df_file data frame to write.
#' @param filename string representing file name.
#' @param delim delimiter for file.
write_df <- function(df_file, filename, delim = "\t") {
  write.table(df_file, row.names = F, sep = delim, file = filename, na = "NA", quote = F)
  return(filename)
}

#' Save cBioPortal file to a destination folder on Synapse, sorting into
#' cBioPortal subfolders as appropriate.  
#' 
#' @param filename path to local file
#' @param synid_folder_output Synapse ID of umbrella output folder
#' @param cleanup if true, delete local file once loaded to Synapse
#' @param prov_name provenance short description
#' @param prov_desc provenance long description
#' @param prov_used vector of Synapse IDs used to create the file
#' @param prov_exec url to script used to create the file
#' @return Synapse ID of stored file entity
save_df <- function(filename, 
                     synid_folder_output, 
                     cleanup = F, 
                     prov_name = "ERBB2 cBioPortal file",
                     prov_desc = "GENIE ERBB2 Sponsored Project data in cBioPortal format",
                     prov_used = NA,
                     prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/main/create_genomic.R") {
  
  synid_folder_dest <- synid_folder_output
  
  if (grepl(pattern = "^case", x = filename)) {
    synid_folder_dest <- create_synapse_folder(name = "case_lists", parent_id = synid_folder_output)
  }
  synid_file_dest <- save_to_synapse(path = filename, 
                  parent_id = synid_folder_dest, 
                  prov_name = prov_name, 
                  prov_desc = prov_desc, 
                  prov_used = prov_used, 
                  prov_exec = prov_exec)
  
  if (cleanup) {
    file.remove(filename)
  }
  
  return(synid_file_dest)
}

# synapse login --------------------

status <- synLogin(auth = auth)

# read ----------------------------

if (verbose) { print(glue("{now()}: reading data and mapping...")) }
df_data <- get_synapse_entity_data_in_csv(synid_file_input, na.strings = c("NA",""))

if (verbose) { print(glue("{now()}: gathering main GENIE files from '{get_synapse_entity_name(synid_folder_mg)}' ({synid_folder_mg})...")) }
synid_files_mg <- get_mg_release_files(synid_folder_mg)

# main ----------------------------

sample_ids <- get_sample_ids(df_data, sample_nos)

for (i in 1:length(filenames)) {
  
  file_type <- names(filenames)[i]
  filename <- as.character(filenames[i])
  
  if (verbose) { print(glue("{now()}: creating data frame for file type '{file_type}'...")) }
  df_file <- create_df(file_type = file_type, synid_files_mg = synid_files_mg, sample_ids = sample_ids)
  
  if (verbose) { print(glue("{now()}: writing data frame for file type '{file_type}' to file '{filename}'...")) }
  filename <- write_df(df_file, file = filename)
  
  if (!is.na(synid_folder_output)) {
    
    if (verbose) { print(glue("{now()}: saving '{filename}' to Synapse folder '{get_synapse_entity_name(synid_folder_output)}' ({synid_folder_output})...")) }
    
    prov_used <- c(synid_file_input) 
    if (!is.na(synid_files_mg[file_type])) {
      prov_used <- append(prov_used, as.character(synid_files_mg[file_type]))
    }
    
    synid_file_df <- save_to_synapse(path = filename, 
                                       parent_id = synid_folder_output, 
                                       prov_name = "ERBB2 cBioPortal file", 
                                       prov_desc = "GENIE ERBB2 Sponsored Project data in cBioPortal format", 
                                       prov_used = prov_used, 
                                       prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/main/create_genomic.R")
    
    if (verbose) { print(glue("{now()}: file '{filename}' saved to {synid_file_df}.")) } 
    
    file.remove(filename)
  }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
