# Description: Create case files for ERBB2 cBioPortal files. 
# Author: Haley Hunter-Zinck
# Date: 2022-05-05

# setup ----------------------------

tic = as.double(Sys.time())

library(glue)
library(dplyr)
library(synapser)
library(optparse)

# constants
cancer_study_identifier <- "genie_erbb2"
filenames <- c("case-breast" = "cases_Breast_Cancer.txt", 
               "case-all" = "cases_all.txt", 
               "case-cna" = "cases_cna.txt", 
               "case-cnaseq" = "cases_cnaseq.txt", 
               "case-seq" = "cases_sequenced.txt")

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

#' Remove leading and trailing whitespace from a string.
#' @param str String
#' @return String without leading or trailing whitespace
trim <- function(str) {
  front <- gsub(pattern = "^[[:space:]]+", replacement = "", x = str)
  back <- gsub(pattern = "[[:space:]]+$", replacement = "", x = front)
  
  return(back)
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

#' Get the name of a Synapse entity. 
#' 
#' @param synapse_id Synapse ID string
#' @return String representing entity name
#' @example get_synapse_entity_name("syn12345")
get_synapse_entity_name <- function(synapse_id) {
  return(synGet(synapse_id, downloadFile = F)$properties$name)
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

#' Create a new Synape Folder entity. 
#' 
#' @param name Name of the Synapse Folder entity
#' @param parentId Synapse ID of Project or Folder in which to create the new Folder
#' @return Synapse ID of the new Synapse Folder entity
create_synapse_folder <- function(name, parent_id) {
  
  # check if folder already exists
  children <- get_synapse_folder_children(parent_id, include_types = list("folder"))
  if(is.element(name, names(children))) {
    return(as.character(children[name]))
  }
  
  concreteType <- "org.sagebionetworks.repo.model.Folder"
  uri <- "/entity"
  payload <- paste0("{", glue("'name':'{name}', 'parentId':'{parent_id}', 'concreteType':'{concreteType}'"), "}")
  ent <- synRestPOST(uri = uri, body = payload)
  return(ent$id)
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


create_case_generic <- function(cancer_study_identifier, stable_id, case_list_name, case_list_description, case_list_ids) {
  
  rows <- rep("", 5)
  rows[1] <- c(glue("cancer_study_identifier: {cancer_study_identifier}"))
  rows[2] <- c(glue("stable_id: {stable_id}"))
  rows[3] <- c(glue("case_list_name: {case_list_name}"))
  rows[4] <- c(glue("case_list_description: {case_list_description}"))
  rows[5] <- c(glue("case_list_ids: {paste0(case_list_ids, collapse = '\t')}"))
  
  return(rows)
}

#' All samples.
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

#' Samples with a breast cancer diagnosis.
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

#' Samples with CNA data
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

#' Samples that have mutation AND CNA data.
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

#' Samples that have mutation OR CNA data.
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

write_case <- function(df_file, filename, delim = "\t") {
  write.table(df_file, row.names = F, sep = delim, file = filename, na = "", quote = F, col.names = F)
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
  }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
