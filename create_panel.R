# Description: Create panel files for ERBB2 cBioPortal files. 
# Format: 
# - https://docs.cbioportal.org/import-gene-panels/#gene-panel-file-format
# Author: Haley Hunter-Zinck
# Date: 2022-05-05

# setup ----------------------------

tic = as.double(Sys.time())

library(glue)
library(dplyr)
library(sqldf)
library(synapser)
library(optparse)
source("shared_fxns.R")

# constants
cancer_study_identifier <- "genie_erbb2"
filename <- "data_gene_panel_{panel_name}.txt"

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
  make_option(c("-i", "--synid_folder_cbio"), type = "character",
              help="Synapse ID of folder containing cBioPortal clinical and genomic files"),
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
waitifnot(!is.null(opt$synid_folder_cbio),
          msg = "Rscript create_case.R -h")

synid_folder_cbio <- opt$synid_folder_cbio
synid_folder_mg <- opt$synid_folder_mg
synid_folder_output <- opt$synid_folder_output
verbose <- opt$verbose
auth <- opt$auth

# functions ----------------------------

#' Get sequencing assay IDs for each sample.
#' 
#' @param synid_files_cbio named vector of strings representing Synapse IDs of cBioPortal files; names represent file labels.
#' @return vector of strings 
get_assays <- function(synid_files_cbio) {
  synid_file_sample <- as.character(synid_files_cbio["sample"])
  df_sample <- get_synapse_entity_data_in_csv(synid_file_sample, sep = "\t", na.strings = "")
  seq_assay_ids <- as.character(unlist(df_sample %>% 
                                         filter(!is.na(SEQ_ASSAY_ID)) %>% 
                                         select(SEQ_ASSAY_ID) %>%
                                         distinct()))
  return(seq_assay_ids)
}

#' Read consolidated genomic region information from main GENIE for sequencing 
#' assay IDs of interest.
#' 
#' @param synid_files_mg named vector of strings representing Synapse IDs of main GENIE cBioPortal files; names represent file labels.
#' @param seq_assay_ids vector of strings representing sequencing assay IDs
#' @return data frame of genomic region information for selected assays
read_gene <- function(synid_files_mg, seq_assay_ids) {
  
  synid_file_gene <- as.character(synid_files_mg["gene"])
  filepath <- synGet(synid_file_gene)$path
  
  # subset gene
  str_ids <- paste0(seq_assay_ids, collapse = "','")
  query <- glue("SELECT * FROM file WHERE SEQ_ASSAY_ID IN ('{str_ids}')")
  df_subset <- read.csv.sql(filepath, 
                            sql = query, 
                            sep = "\t",
                            colClasses = c("character",rep("integer",2), rep("character", 5)))
  
  return(df_subset)
}

#' Create a vector representing each row of a generic cBioPortal panel file given the values for all fields.
#' 
#' @param stable_id name of the gene panel
#' @param description description of the gene panel
#' @param gene_list tab separated genes, represented either by all gene symbols or all Entrez gene IDs
#' @return vector of strings representing each row of the panel file.
create_panel_generic <- function(stable_id, description, gene_list) {
  
  rows <- rep("", 3)
  rows[1] <- c(glue("stable_id: {stable_id}"))
  rows[2] <- c(glue("description: {description}"))
  rows[3] <- c(glue("gene_list: {paste0(gene_list, collapse = '\t')}"))
  
  return(rows)
}

#' Create requested gene panel file.
#' 
#' @param seq_assay_id string representing 
#' @param df_gene data frame represening relevant main GENIE genomic information.
#' @return vector of strings representing each row of the panel file.
create_panel <- function(seq_assay_id, df_gene) {
  
  gene_list <- as.character(unlist(df_gene %>% 
    filter(SEQ_ASSAY_ID == seq_assay_id) %>%
    select(Hugo_Symbol) %>%
    distinct()))
  
  rows <- create_panel_generic(stable_id = seq_assay_id, 
                               description = glue("{seq_assay_id}, Number of Genes - {length(gene_list)}"), 
                               gene_list = gene_list)
  
  return(rows)
}

#' Write the panel file content to file.
#' 
#' @param df_file data frame of data to write
#' @param filename name of the output file
#' @param delim text delimiter for output file.
#' @return name of output file.
write_panel <- function(df_file, filename, delim = "\t") {
  write.table(df_file, row.names = F, sep = delim, file = filename, na = "", quote = F, col.names = F)
  return(filename)
}

# synapse login --------------------

status <- synLogin(auth = auth)

# read ----------------------------

if (verbose) { print(glue("{now()}: gathering cBioPortal files from '{get_synapse_entity_name(synid_folder_cbio)}' ({synid_folder_cbio})...")) }
synid_files_cbio <- get_cbio_release_files(synid_folder_cbio)

if (verbose) { print(glue("{now()}: gathering main GENIE files from '{get_synapse_entity_name(synid_folder_mg)}' ({synid_folder_mg})...")) }
synid_files_mg <- get_mg_release_files(synid_folder_mg)

if (verbose) { print(glue("{now()}: read genomic information...")) }
seq_assay_ids <- get_assays(synid_files_cbio)
df_gene <- read_gene(synid_files_mg, seq_assay_ids)

# main ----------------------------

for (i in 1:length(seq_assay_ids)) {
  
  seq_assay_id <- seq_assay_ids[i]
  filename <- glue("data_gene_panel_{seq_assay_id}.txt")
  
  if (verbose) { glue("{now()}: creating panel file for '{seq_assay_id}'...") }
  df_file <- create_panel(seq_assay_id, df_gene)
  
  if (verbose) { glue("{now()}: writing panel file for '{seq_assay_id}' to '{filename}'...") }
  filename <- write_panel(df_file, file = filename)
  
  if (!is.na(synid_folder_output)) {
    
    # create subfolder to store gene panels
    synid_folder_panel <- create_synapse_folder(name = "gene_panels", parent_id = synid_folder_output)
    
    if (verbose) { print(glue("{now()}: saving '{filename}' to Synapse folder '{get_synapse_entity_name(synid_folder_panel)}' ({synid_folder_panel})...")) }
    
    synid_file_df <- save_to_synapse(path = filename, 
                                     parent_id = synid_folder_panel, 
                                     prov_name = "ERBB2 cBioPortal file", 
                                     prov_desc = "GENIE ERBB2 Sponsored Project data in cBioPortal format", 
                                     prov_used = as.character(synid_files_mg["gene"]), 
                                     prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/main/create_panel.R")
    
    if (verbose) { print(glue("{now()}: file '{filename}' saved to {synid_file_df}.")) } 
    
    file.remove(filename)
  }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
