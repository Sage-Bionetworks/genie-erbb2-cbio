# Description: Transform uncoded ERBB2 REDCap export to cBioPortal format
#               to create patients, sample, and timeline files. 
# Author: Haley Hunter-Zinck
# Date: 2022-05-03

# setup ----------------------------

tic = as.double(Sys.time())

library(glue)
library(dplyr)
library(sqldf)
library(synapser)
library(optparse)

# constants
sample_nos <- c(1:10)
filenames <- c("maf" = "data_mutations_extended.txt", 
               "cna" = "data_CNA.txt", 
               "matrix" = "data_gene_matrix.txt", 
               "fusion" = "data_fusions.txt", 
               "seg" = "genie_erbb2_data_cna_hg19.seg")

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

#' Gather main GENIE files from a designated release folder.
#' 
#' @param synid_folder_mg Synapse ID of main GENIE release
#' @return named vector of Synapse IDs corresponding to labeled files. 
get_mg_release_files <- function(synid_folder_mg) {
  synid_folder_children <- get_synapse_folder_children(synapse_id = synid_folder_mg, 
                                                       include_types=list("file"))
  synid_files_mg <- c("maf" = as.character(synid_folder_children["data_mutations_extended.txt"]), 
                      "cna" = as.character(synid_folder_children["data_CNA.txt"]), 
                      "matrix" = as.character(synid_folder_children["data_gene_matrix.txt"]), 
                      "fusion" = as.character(synid_folder_children["data_fusions.txt"]),
                      "seg" = as.character(synid_folder_children["genie_data_cna_hg19.seg"]),
                      "gene" = as.character(synid_folder_children["genomic_information.txt"]))
  
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

#' If a string does not begin with the prefix "GENIE-", add it to
#' the beginning of the string.
#' 
#' @param x string
#' @return string with prefix "GENIE-"
add_genie_prefix <- function(x) {
  if (!grepl(pattern = "^GENIE-", x = x)) {
    return(glue("GENIE-{x}"))
  }
  return(x)
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

create_maf <- function(synid_file_mg, sample_ids) {
  
  filepath <- synGet(synid_file_mg)$path
  
  # subset maf
  str_ids <- paste0(sample_ids, collapse = "','")
  query <- glue("SELECT * FROM file WHERE Tumor_Sample_Barcode IN ('{str_ids}')")
  df_subset <- read.csv.sql(filepath, sql = query, sep = "\t")
  
  # custom formatting for certain columns
  for (colname in grep(pattern = "^[nt]_", x = colnames(df_subset), value = T)) {
    df_subset[which(df_subset[,colname] == "."),colname] <- ""
  }
  df_subset[, "Validation_Status"] <- ""
    
  return(df_subset)
}

create_cna <- function(synid_file_mg, sample_ids, n_row_per_chunk = 50) {
  
  df_subset <- c()
  filepath <- synGet(synid_file_mg)$path
  
  # determine samples in file
  header_cna <- scan(filepath, nlines = 1, what = "character", quiet = T)
  both_sample_ids <- intersect(header_cna, sample_ids)
  idx <- c(which(header_cna == "Hugo_Symbol"), which(is.element(header_cna, sample_ids)))
  
  # subset cna in chunks
  n_row_total <- as.double(strsplit(system(glue("wc -l {filepath}"), intern = T), split = " ")[[1]][1])
  df_subset <- matrix("NA", nrow = n_row_total - 1, ncol = length(idx), dimnames = list(c(), header_cna[idx]))
  for (i in seq(from = 2, to = n_row_total, by = n_row_per_chunk)) {
    chunk <- matrix(scan(filepath, nlines = n_row_per_chunk, skip = i-1, what = "character", quiet = T), ncol = length(header_cna), byrow = T)[,idx]
    df_subset[(i-1):(i+nrow(chunk)-2),] <- chunk
  }
  
  return(df_subset)
}

create_matrix <- function(synid_file_mg, sample_ids) {
  df_mg <- get_synapse_entity_data_in_csv(synid_file_mg, sep = "\t", na.strings = "")
  df_subset <- df_mg %>% 
    filter(is.element(SAMPLE_ID, sample_ids))
  df_subset[which(is.na(df_subset), arr.ind = T)] <- "NA"
  return(df_subset)
}

create_fusion <- function(synid_file_mg, sample_ids) {
  df_mg <- get_synapse_entity_data_in_csv(synid_file_mg, sep = "\t")
  df_subset <- df_mg %>% 
    filter(is.element(Tumor_Sample_Barcode, sample_ids))
  return(df_subset)
}

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

create_panel <- function(synid_file_mg, sample_ids) {
  return(NULL)
}

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
                             sample_ids = sample_ids)
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

write_df <- function(df_file, filename, delim = "\t") {
  write.table(df_file, row.names = F, sep = delim, file = filename, na = "", quote = F)
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
  }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
