# Description: Add non-GENIE data to cBioPortal files for ERBB2.
# Author: Haley Hunter-Zinck
# Date: 2022-05-06

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
  make_option(c("-i", "--synid_folder_cbio"), type = "character",
              help="Synapse ID of folder containing cBioPortal files"),
  make_option(c("-o", "--synid_folder_output"), type = "character",
              help="Synapse ID of output folder"),
  make_option(c("-s", "--synid_file_ng_sam"), type = "character", default = "syn30041987",
              help="Synapse ID of file containing non-GENIE clinical sample info (default: syn30041987)"),
  make_option(c("-b", "--synid_file_ng_bed"), type = "character", default = "syn30135792",
              help="Synapse ID of file containing non-GENIE bed info  (default: syn30135792)"),
  make_option(c("-m", "--synid_file_ng_maf"), type = "character", default = "syn30041988",
              help="Synapse ID of file containing non-GENIE mutation info  (default: syn30041988)"),
  make_option(c("-v", "--verbose"), action="store_true", default = FALSE, 
              help="Output script messages to the user."),
  make_option(c("-a", "--auth"), 
              type = "character",
              default = NA,
              help="Synapse personal access token or path to .synapseConfig (default: normal synapse login behavior)")
)
opt <- parse_args(OptionParser(option_list=option_list))
waitifnot(!is.null(opt$synid_file_input) && !is.null(opt$synid_folder_output),
          msg = "Rscript add_nongenie.R -h")

synid_folder_cbio <- opt$synid_folder_cbio
synid_folder_output <- opt$synid_folder_output
synid_file_ng_sam <- opt$synid_file_ng_sam
synid_file_ng_bed <- opt$synid_file_ng_bed
synid_file_ng_maf <- opt$synid_file_ng_maf
verbose <- opt$verbose
auth <- opt$auth

# functions ----------------------------

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

#' Get all child entities of a synapse folder.
#' 
#' @param synapse_id Synapse ID of the folder
#' @param include_types Types of child entities to return
#' @return Vector with values as Synapse IDs and names as entity names.
get_synapse_folder_children <- function(synapse_id, 
                                        include_types=list("folder", "file", "table", "link", "entityview", "dockerrepo")) {
  
  ent <- as.list(synGetChildren(synapse_id, includeTypes = include_types))
  
  children <- c()
  
  if (length(ent) > 0) {
    for (i in 1:length(ent)) {
      children[ent[[i]]$name] <- ent[[i]]$id
    }
  }
  
  return(children)
}

#' Gather cBioPortal files from a designated release folder.
#' 
#' @param synid_folder_cbio Synapse ID of main GENIE release
#' @return named vector of Synapse IDs corresponding to labeled files. 
get_cbio_release_files <- function(synid_folder_cbio) {
  synid_folder_children <- get_synapse_folder_children(synapse_id = synid_folder_cbio, 
                                                       include_types=list("file"))
  synid_files_mg <- c("maf" = as.character(synid_folder_children["data_mutations_extended.txt"]), 
                      "cna" = as.character(synid_folder_children["data_CNA.txt"]), 
                      "matrix" = as.character(synid_folder_children["data_gene_matrix.txt"]), 
                      "fusion" = as.character(synid_folder_children["data_fusions.txt"]),
                      "seg" = as.character(synid_folder_children["genie_erbb2_data_cna_hg19.seg"]),
                      "patient" = as.character(synid_folder_children["data_clinical_patient.txt"]),
                      "sample" = as.character(synid_folder_children["data_clinical_sample.txt"]),
                      "timeline" = as.character(synid_folder_children["data_timeline.txt"]))
  
  return(synid_files_mg)
}

read_gene <- function(synid_file_gene) {
  df_raw <- get_synapse_entity_data_in_csv(synid_file_ng_bed, sep = "\t", na.strings = c("NA",""), header = F)
  colnames(df_raw) <- c("chr", "Start_Position", "End_Position", "label", "SEQ_ASSAY_ID")
  df_gene <- df_raw %>%
    mutate(Chromosome = gsub(chr, pattern = "chr", replacement = "")) %>%
    mutate(ID = unlist(lapply(strsplit(x = label, split = "_"), head, n = 1))) %>%
    select(Chromosome, Start_Position, End_Position, ID, SEQ_ASSAY_ID)
  return(df_raw)
}

modify_sample <- function(df_mg_sam, df_ng_sam) {
  return(NULL)
}

modify_maf <- function(df_mg_maf, df_ng_maf) {
  return(NULL)
}

modify_matrix <- function(df_mg_matrix, df_ng_maf) {
  return(NULL)
}

create_panel_generic <- function(stable_id, description, gene_list) {
  
  rows <- rep("", 3)
  rows[1] <- c(glue("stable_id: {stable_id}"))
  rows[2] <- c(glue("description: {description}"))
  rows[3] <- c(glue("gene_list: {paste0(gene_list, collapse = '\t')}"))
  
  return(rows)
}

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

# synapse login --------------------

status <- synLogin(auth = auth)

# read ----------------------------

if (verbose) { print(glue("{now()}: reading non-GENIE files...")) }

df_ng_sam <- get_synapse_entity_data_in_csv(synid_file_ng_sam, sep = "\t", na.strings = c("NA",""))
df_ng_gene <- read_gene(synid_file_ng_bed)
df_ng_maf <- get_synapse_entity_data_in_csv(synid_file_ng_maf, sep = "\t", na.strings = c("NA",""))

if (verbose) { print(glue("{now()}: gathering GENIE cBioPortal files from '{get_synapse_entity_name(synid_folder_cbio)}' ({synid_folder_cbio})...")) }
synid_files_cbio <- get_cbio_release_files(synid_folder_cbio)

if (verbose) { print(glue("{now()}: reading relevant GENIE cBioPortal files...")) }
df_mg_sam <- get_synapse_entity_data_in_csv(synid_files_cbio["sample"], sep = "\t")
df_mg_maf <- get_synapse_entity_data_in_csv(synid_files_cbio["maf"], sep = "\t")
df_mg_matrix <- get_synapse_entity_data_in_csv(synid_files_cbio["maf"], sep = "\t")

# main ----------------------------

# modify clinical sample, maf, gene matrix files
df_sam <- modify_sample(df_mg_sam = df_mg_sam, df_ng_sam = df_ng_sam)
df_sam <- modify_maf(df_mg_maf = df_mg_maf, df_ng_maf = df_ng_maf)
df_sam <- modify_matrix(df_mg_matrix = df_mg_matrix, df_ng_maf = df_ng_maf)

# create new gene panel files
df_gene <- create_panel(seq_assay_id, df_gene)

# write modified files locally

# if storing to synapse
#   copy GENIE cBioPortal entities to new folder
#   save modified files to synapse
#   remove local versions

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
