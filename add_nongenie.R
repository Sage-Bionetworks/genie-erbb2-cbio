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
  make_option(c("-o", "--synid_folder_output"), type = "character", default = NA,
              help="Synapse ID of output folder (default: write locally only)"),
  make_option(c("-s", "--synid_file_ng_sam"), type = "character", 
              help="Synapse ID of file containing non-GENIE clinical sample info"),
  make_option(c("-b", "--synid_file_ng_bed"), type = "character", 
              help="Synapse ID of file containing non-GENIE bed info"),
  make_option(c("-m", "--synid_file_ng_maf"), type = "character",
              help="Synapse ID of file containing non-GENIE mutation info "),
  make_option(c("-v", "--verbose"), action="store_true", default = FALSE, 
              help="Output script messages to the user."),
  make_option(c("-a", "--auth"), 
              type = "character",
              default = NA,
              help="Synapse personal access token or path to .synapseConfig (default: normal synapse login behavior)")
)
opt <- parse_args(OptionParser(option_list=option_list))
waitifnot(!is.null(opt$synid_folder_cbio) && 
            !is.null(opt$synid_file_ng_sam) && 
            !is.null(opt$synid_file_ng_bed) &&
            !is.null(opt$synid_file_ng_maf),
          msg = "Rscript add_nongenie.R -h")

synid_folder_cbio <- opt$synid_folder_cbio
synid_folder_output <- opt$synid_folder_output
verbose <- opt$verbose
auth <- opt$auth

synid_files_ng <- c("sample" = opt$synid_file_ng_sam,
                    "bed" = opt$synid_file_ng_bed,
                    "maf" = opt$synid_file_ng_maf)

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

#' Get the name of a Synapse entity. 
#' 
#' @param synapse_id Synapse ID string
#' @return String representing entity name
#' @example get_synapse_entity_name("syn12345")
get_synapse_entity_name <- function(synapse_id) {
  return(synGet(synapse_id, downloadFile = F)$properties$name)
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

#' Read in the non-GENIE clinical sample file and filter for duplicates.
#' 
#' @param synid_file_sam Synapse ID of clinical sample file.
#' @return Data frame representation of file's data with duplicate rows removed.
create_merged_sample <- function(synid_file_mg, synid_file_ng) {
  
  df_mg <- get_synapse_entity_data_in_csv(synid_file_mg, sep = "\t", na.strings = c("NA",""))
  df_raw <- get_synapse_entity_data_in_csv(synid_file_ng, sep = "\t", na.strings = c("NA",""))
  df_raw <- df_raw %>%
    distinct()
  
  df_mod <- df_mg %>%
    left_join(df_raw, by = c("SAMPLE_ID", "PATIENT_ID")) %>%
    mutate(SEQ_ASSAY_ID = case_when (!is.na(SEQ_ASSAY_ID.y) ~ SEQ_ASSAY_ID.y,
                                     TRUE ~ SEQ_ASSAY_ID.x)) %>%
    select(-SEQ_ASSAY_ID.x, -SEQ_ASSAY_ID.y)
  
  return(df_mod)
}

#' Append the non-GENIE maf file data to the GENIE maf file data
#' for the ERBB2 samples. 
#' 
#' @param synid_file_mg Synapse ID of ERBB2 maf file derived from main GENIE data.
#' @param synid_file_ng Synapse ID of ERBB2 maf file derived from non-GENIE data.
#' @return data frame representing ERBB2 maf file for main GENIE and non-GENIE data.
create_merged_maf <- function(synid_file_mg, synid_file_ng) {
  
  df_mg <- get_synapse_entity_data_in_csv(synid_file_mg, sep = "\t", na.strings = c("NA",""))
  df_ng <- get_synapse_entity_data_in_csv(synid_file_ng, sep = "\t", na.strings = c("NA",""))
  
  df_maf <- df_mg %>% bind_rows(df_ng)
  
  # custom formatting for certain columns
  for (colname in grep(pattern = "^[nt]_", x = colnames(df_maf), value = T)) {
    df_maf[which(df_maf[,colname] == "."),colname] <- ""
  }
  df_maf[, "Validation_Status"] <- ""
  
  return(df_maf)
}

create_merged_matrix <- function(df_maf, df_cna, df_sam) {
  maf_assay <- df_maf %>% 
    rename(SAMPLE_ID = Tumor_Sample_Barcode) %>%
    select(SAMPLE_ID) %>%
    distinct() %>%
    left_join(df_sam, by = "SAMPLE_ID") %>%
    rename(mutations = SEQ_ASSAY_ID) %>%
    select(SAMPLE_ID, mutations)
  cna_assay <- data.frame(SAMPLE_ID = colnames(df_cna)[2:ncol(df_cna)]) %>% 
    select(SAMPLE_ID) %>%
    distinct() %>%
    left_join(df_sam, by = "SAMPLE_ID") %>%
    rename(cna = SEQ_ASSAY_ID) %>%
    select(SAMPLE_ID, cna)
  
  sample_ids <- union(maf_assay$SAMPLE_ID, cna_assay$SAMPLE_ID)
  df_matrix <- data.frame(SAMPLE_ID = sample_ids) %>%
    left_join(maf_assay, by = "SAMPLE_ID") %>%
    left_join(cna_assay, by = "SAMPLE_ID") %>%
    filter(!is.na(mutations) | !is.na(cna))
  
  return(df_matrix)
}

#' Get SEQ ASSAY ID given that it is the first part of the file.
#' @param synid_file_ng Synapse ID of the bed file.  
#' @return character string matching a SEQ_ASSAY_ID in clinical sample file. 
get_seq_assay_id <- function(synid_file_ng) {
  filename <- get_synapse_entity_name(synid_file_ng)
  return(gsub(x = filename, pattern = "-ASCII.bed$", replacement = ""))
}

create_nongenie_bed <- function(synid_file_ng, seq_assay_id) {
  df_raw <- get_synapse_entity_data_in_csv(synid_file_ng, sep = "\t", na.strings = c("NA",""), header = F)
  colnames(df_raw) <- c("chr", "Start_Position", "End_Position", "label", "unknown")
  df_bed <- df_raw %>%
    mutate(Chromosome = gsub(chr, pattern = "chr", replacement = "")) %>%
    mutate(Hugo_Symbol = unlist(lapply(strsplit(x = label, split = "_"), head, n = 1))) %>%
    mutate(SEQ_ASSAY_ID = seq_assay_id) %>%
    select(Chromosome, Start_Position, End_Position, Hugo_Symbol, SEQ_ASSAY_ID)
  return(df_bed)
}

create_panel_generic <- function(stable_id, description, gene_list) {
  
  rows <- rep("", 3)
  rows[1] <- c(glue("stable_id: {stable_id}"))
  rows[2] <- c(glue("description: {description}"))
  rows[3] <- c(glue("gene_list: {paste0(gene_list, collapse = '\t')}"))
  
  return(rows)
}

create_panel <- function(seq_assay_id, df_bed) {
  
  gene_list <- as.character(unlist(df_bed %>% 
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

get_cbio_sample_header <- function(synid_file, nlines = 5, delim = "\t") {
  
  filepath <- synGet(synid_file)$path
  header <- c()
  
  for (i in 1:nlines) {
    header <- rbind(header, scan(filepath, nlines = 1, skip = i-1, what = "character", sep = "\t", quote = "", quiet = T))
  }
  
  return(header)
}

write_cbio_clinical <- function(df_file, file_name, header, delim = "\t") {
  
  colnames(header) <- header[5,]
  df_write <- rbind(header, df_file[colnames(header)])
  write.table(df_write, file = file_name, sep = delim, na = "",
              col.names = F, row.names = F, quote = F)
  
  return(file_name)
}

write_df <- function(df_file, filename, delim = "\t", na = "NA") {
  write.table(df_file, row.names = F, sep = delim, file = filename, na = na, quote = F)
  return(filename)
}

write_panel <- function(df_file, filename, delim = "\t") {
  write.table(df_file, row.names = F, sep = delim, file = filename, na = "", quote = F, col.names = F)
  return(filename)
}

# synapse login --------------------

status <- synLogin(auth = auth)

# read ----------------------------

if (verbose) { print(glue("{now()}: gathering new cBioPortal files from '{get_synapse_entity_name(synid_folder_cbio)}' ({synid_folder_cbio})...")) }
synid_files_cbio <- get_cbio_release_files(synid_folder_cbio)

# main ----------------------------

dfs <- list()

if (verbose) { print(glue("{now()}: merging non-GENIE sample info...")) }

dfs[["sample"]] <- create_merged_sample(synid_file_mg = as.character(synid_files_cbio["sample"]),
                               synid_file_ng = as.character(synid_files_ng["sample"]))

if (verbose) { print(glue("{now()}: merging non-GENIE maf file...")) }

dfs[["maf"]] <- create_merged_maf(synid_file_mg = as.character(synid_files_cbio["maf"]),
                               synid_file_ng = as.character(synid_files_ng["maf"]))

if (verbose) { print(glue("{now()}: merging non-GENIE matrix info...")) }

df_cna <- get_synapse_entity_data_in_csv(as.character(synid_files_cbio["cna"]), sep = "\t", na.strings = c("NA",""))
dfs[["matrix"]] <- create_merged_matrix(df_maf = dfs[["maf"]], 
                                  df_cna = df_cna, 
                                  df_sam = dfs[["sample"]])

if (verbose) { print(glue("{now()}: creating non-GENIE panel file...")) }

# create new gene panel files
seq_assay_id <- get_seq_assay_id(as.character(synid_files_ng["bed"]))
df_bed <- create_nongenie_bed(synid_file_ng = as.character(synid_files_ng["bed"]), 
                                seq_assay_id = seq_assay_id)
dfs[["panel"]] <- create_panel(seq_assay_id, df_bed)

# write ---------------------

if (verbose) { print(glue("{now()}: writing files locally...")) }

# write modified files locally
outfiles <- c()
for (file_type in setdiff(names(dfs), "panel")) {
  
  if (file_type == "panel") {
    outfiles[file_type] <- write_panel(df_file = df_panel, filename = glue("data_gene_panel_{seq_assay_id}.txt"))
  } else if (file_type == "sample") {
    entity <- synGet(as.character(synid_files_cbio[[file_type]]))
    file_name <- entity$properties$name
    header <- get_cbio_sample_header(entity$properties$id)
    
    outfiles[file_type] <- write_cbio_clinical(df_file = dfs[[file_type]], 
                                               file_name = file_name, 
                                               header = header)
  } else if (file_type == "maf") {
    outfiles[file_type] <- write_df(df_file = dfs[[file_type]], 
                                    filename = get_synapse_entity_name(as.character(synid_files_cbio[[file_type]])), na = "")
  } else {
    outfiles[file_type] <- write_df(df_file = dfs[[file_type]], 
                                    filename = get_synapse_entity_name(as.character(synid_files_cbio[[file_type]])))
  }
}

# store ---------------------

# if requested, store to synapse
if (!is.na(synid_folder_output)) {
  for (file_type in names(outfiles)) {
    
    synid_folder_dest <- ""
    if (file_type == "panel") {
      synid_folder_dest <- create_synapse_folder(name = "gene_panels", parent_id = synid_folder_output)
      used <- c(synid_files_ng["bed"])
    } else {
      synid_folder_dest <- synid_folder_output
      if (!is.na(synid_files_ng[file_type])) {
        used <- c(synid_files_cbio[file_type], synid_files_ng[file_type])
      } else {
        used <- c(synid_files_cbio[file_type])
      }
    }
    
    if (verbose) { print(glue("{now()}: saving '{outfiles[file_type]}' to Synapse folder '{get_synapse_entity_name(synid_folder_dest)}' ({synid_folder_dest})...")) }
    
    synid_file_df <- save_to_synapse(path = as.character(outfiles[file_type]), 
                                     parent_id = synid_folder_dest, 
                                     prov_name = "ERBB2 cBioPortal file plus", 
                                     prov_desc = "GENIE ERBB2 Sponsored Project data in cBioPortal format with non-GENIE data", 
                                     prov_used = as.character(used), 
                                     prov_exec = "https://github.com/Sage-Bionetworks/genie-erbb2-cbio/blob/main/add_nongenie.R")
    
    if (verbose) { print(glue("{now()}: file '{outfiles[file_type]}' saved to {synid_file_df}.")) } 
    
    file.remove(outfiles[file_type])
  }
}

# close out ----------------------------

toc = as.double(Sys.time())
print(glue("Runtime: {round(toc - tic)} s"))
