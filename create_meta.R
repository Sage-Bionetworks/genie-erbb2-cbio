# Description: create cBioPortal meta data files for each file type.  
# Author: Haley Hunter-Zinck
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

create_meta_clinical <- function() {
  return(NULL)
}

create_meta_genomic <- function() {
  return(NULL)
}

create_meta_seg <- function() {
  return(NULL)
}

create_meta_study <- function() {
  return(NULL)
}

create_meta <- function(file_type) {
  
  df_type <- NULL
  
  if (file_type == "study") {
    df_file <- create_meta_study()
  } else if (file_type == "seg") {
    df_file <- create_meta_seg()
  } else if (is.element(file_type, c("patient", "sample", "timeline"))) {
    df_file <- create_meta_clinical()
  } else if (is.element(file_type, c("cna", "fusion", "matrix", "maf"))) {
    df_file <- create_meta_genomic()
  }
  
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
                                          prov_desc = "", 
                                          prov_used = NA, 
                                          prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/develop/create_meta.R")
                        }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
