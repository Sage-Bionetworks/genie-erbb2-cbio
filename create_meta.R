# Description: create cBioPortal meta data files for each file type.  
# Author: Haley Hunter-Zinck
# Based on: https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects/blob/develop/geniesp/metafiles.py
# Date: 2022-05-12

# setup ----------------------------

tic = as.double(Sys.time())

library(glue)
library(dplyr)
library(synapser)
library(optparse)
source("shared_fxns.R")

# constants
filenames <- c("seg" = "genie_erbb2_meta_cna_hg19_seg.txt",
               "cna" = "meta_CNA.txt",
               "patient" = "meta_clinical_patient.txt", 
               "sample" = "meta_clinical_sample.txt", 
               "fusion" = "meta_fusions.txt", 
               "matrix" = "meta_gene_matrix.txt", 
               "maf" = "meta_mutations_extended.txt",
               "study" = "meta_study.txt",
               "timeline" = "meta_timeline.txt")

# user input ----------------------------

option_list <- list( 
  make_option(c("-i", "--synid_file_input"), type = "character",
              help="Synapse ID of input file"),
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
          msg = "Rscript create_meta.R -h")

synid_file_input <- opt$synid_file_input
synid_folder_output <- opt$synid_folder_output
verbose <- opt$verbose
auth <- opt$auth


# functions ----------------------------

create_meta_clinical_generic <- function(cancer_study_identifier, genetic_alteration_type, datatype, data_filename) {
  rows <- rep("", 4)
  rows[1] <- c(glue("cancer_study_identifier: {cancer_study_identifier}"))
  rows[2] <- c(glue("genetic_alteration_type: {genetic_alteration_type}"))
  rows[3] <- c(glue("datatype: {datatype}"))
  rows[4] <- c(glue("data_filename: {data_filename}"))
  return(NULL)
}

create_meta_genomic_generic <- function(cancer_study_identifier, genetic_alteration_type, datatype, stable_id, profile_name, profile_description, data_filename) {
  
  rows <- rep(NA, 8)
  rows[1] <- glue("cancer_study_identifier: {cancer_study_identifier}")
  rows[2] <- glue("genetic_alteration_type: {genetic_alteration_type}")
  rows[3] <- glue("datatype: {datatype}")
  rows[4] <- glue("stable_id: {stable_id}")
  rows[5] <- glue("show_profile_in_analysis_tab: true")
  rows[6] <- glue("profile_name: {profile_name}")
  rows[7] <- glue("profile_description: {profile_description}")
  rows[8] <- glue("data_filename: {data_filename}")
  return(rows)
}

create_meta_seg_generic <- function(cancer_study_identifier, genetic_alteration_type, datatype, reference_genome_id, description, data_filename) {
  rows <- rep(NA, )
  rows[1] <- glue("cancer_study_identifier: {cancer_study_identifier}")
  rows[2] <- glue("genetic_alteration_type: {genetic_alteration_type}")
  rows[3] <- glue("datatype: {datatype}")
  rows[4] <- glue("reference_genome_id: {reference_genome_id}")
  rows[5] <- glue("description: {description}")
  rows[6] <- glue("data_filename: {data_filename}") 
  return(rows)
}

create_meta_study_generic <- function(type_of_cancer, cancer_study_identifier, name, description, groups, short_name) {
  rows <- rep(NA, 6)
  rows[1] <- glue("type_of_cancer: {type_of_cancer}")
  rows[2] <- glue("cancer_study_identifier: {cancer_study_identifier}")
  rows[3] <- glue("name: {name}")
  rows[4] <- glue("description: {description}")
  rows[5] <- glue("groups: {groups}")
  rows[6] <- glue("short_name: {short_name}")
  return(rows)
}

create_meta_patient <- function(cancer_study_identifier, data_filename = "meta_clinical_patient.txt") {
  df_file <- create_meta_clinical_generic(cancer_study_identifier = cancer_study_identifier, 
                               genetic_alteration_type = "CLINICAL", 
                               datatype = "PATIENT_ATTRIBUTES", 
                               data_filename = data_filename)
  return(df_file)
}

create_meta_sample <- function(cancer_study_identifier, data_filename = "data_clinical_sample.txt") {
  df_file <- create_meta_clinical_generic(cancer_study_identifier = cancer_study_identifier, 
                                          genetic_alteration_type = "CLINICAL", 
                                          datatype = "SAMPLE_ATTRIBUTES", 
                                          data_filename = data_filename)
  return(df_file)
}

create_meta_timeline <- function(cancer_study_identifier, data_filename = "data_timeline.txt") {
  df_file <- create_meta_clinical_generic(cancer_study_identifier = cancer_study_identifier, 
                                          genetic_alteration_type = "CLINICAL", 
                                          datatype = "TIMELINE", 
                                          data_filename = data_filename)
  return(df_file)
}

create_meta_matrix <- function(cancer_study_identifier, data_filename = "data_genie_matrix.txt") {
  df_file <- create_meta_clinical_generic(cancer_study_identifier = cancer_study_identifier, 
                                          genetic_alteration_type = "GENE_PANEL_MATRIX", 
                                          datatype = "GENE_PANEL_MATRIX", 
                                          data_filename = data_filename)
  return(df_file)
}

create_meta_cna <- function(cancer_study_identifier, data_filename = "data_CNA.txt") {
  df_file <- create_meta_genomic_generic(cancer_study_identifier, 
                                         genetic_alteration_type = "COPY_NUMBER_ALTERATION", 
                                         datatype = "DISCRETE", 
                                         stable_id = "cna", 
                                         profile_name = "Copy-number alterations", 
                                         profile_description = "Copy-number alterations", 
                                         data_filename)
  
  return(df_file)
}

create_meta_fusion <- function(cancer_study_identifier, data_filename = "data_CNA.txt") {
  df_file <- create_meta_genomic_generic(cancer_study_identifier = cancer_study_identifier, 
                                         genetic_alteration_type = "FUSION", 
                                         datatype = "FUSION", 
                                         stable_id = "fusion", 
                                         profile_name = "Fusions", 
                                         profile_description = "Fusions",
                                         data_filename = data_filename)
  return(df_file)
}


create_meta_maf <- function() {
  df_file <- create_meta_genomic_generic(cancer_study_identifier = cancer_study_identifier, 
                                         genetic_alteration_type = "MUTATION_EXTENDED", 
                                         datatype = "MAF", 
                                         stable_id = "mutations",
                                         profile_name = "Mutations",
                                         profile_description = "Mutation data from next-gen sequencing.",
                                         data_filename = data_filename)
  return(df_file)
}

create_meta_seg <- function(cancer_study_identifier, data_filename = "") {
  df_file <- create_meta_seg_generic(cancer_study_identifier = cancer_study_identifier, 
                                     genetic_alteration_type = "COPY_NUMBER_ALTERATION", 
                                     datatype = "SEG", 
                                     reference_genome_id = "hg19", 
                                     description = "Segment data for the genie study", 
                                     data_filename = data_filename) {
    
  return(df_file)
}

create_meta_study <- function() {
  df_file <- create_meta_study_generic(type_of_cancer = type_of_cancer, 
                                       cancer_study_identifier, 
                                       name, 
                                       description, 
                                       groups, 
                                       short_name)
  return(df_file)
}

#' Map function names to implementations.
#' 
#' @param file_type string representing cBio file type labels
#' @param file_types vector of string representing all valid file type labels
#' @return object representing function
get_fxn <- function(file_type, 
                    file_types = c("patient", "sample", "timeline", "cna", "matrix", "maf", "seg", "study")) {
  
  map_fxn_str <- list()
  names_fxn <- paste0("create_meta_", file_types)
  
  for (name_fxn in names_fxn) {
    map_fxn_str[[name_fxn]] <- get(name_fxn)
  }
  
  return(map_fxn_str[[glue("create_meta_{file_type}")]])
}

create_meta <- function(file_type) {
  
  fxn <- get_fxn(file_type)
  df_type <- fxn()
  return(df_type)
}

#' Write the vector of strings representing the rows of a meta file to a local file.
#' 
#' @param df_file data to write.
#' @param filename name of the output file.
#' @return string representing name of output file.
write_meta <- function(df_file, filename) {
  write.table(df_file, row.names = F, file = filename, na = "", quote = F, col.names = F)
  return(filename)
}

# synapse login --------------------

status <- synLogin(auth = auth)

# read ----------------------------


# main ----------------------------

for (i in 1:length(filenames)) {
  file_name <- filenames[i]
  file_type <- names(filenames)[i]
  
  df_file <- create_meta(file_type)
  write_meta(df_file, filename) 
  
  if (!is.na(synid_folder_output)) {
    synid_file_output <- save_to_synapse(path = filename, 
                                          parent_id = synid_folder_output, 
                                          file_name = filename, 
                                          prov_name = "cBioPortal meta file", 
                                          prov_desc = "cBioPortal meta files for annotating provided cBioPortal data files", 
                                          prov_used = NA, 
                                          prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/develop/create_meta.R")
                        }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
