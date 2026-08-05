#!/usr/bin/env Rscript
# DATE:       20250730
# AUTHOR:     Maximiliano Zuluaga Forero
# SCRIPT:     03_compile_all_DMRs_V04.R
# VERSION:    04

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DESCRIPTION:
#             This script recursively searches the methylasso DMR output
#             folders to collect all *_dmrs.tsv files and compiles them 
#             into a unified file named "dmr_all.tsv". It adds metadata
#             columns: condition, comparison/group, source file and chromosome.
# 
#             
# VERSIONS:
#             V02:  Adds Chromosome type column 
#             V03:  Do not use 'parameter as part of the folder structure'
#             V04:  works! 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load libraries
library(dplyr)
library(readr)
library(stringr)
library(purrr)

# Set path to DMR output data folder
base_dir <- "/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/02_DMR/data/01_methylasso/methylasso_DMRs"
setwd(base_dir)
# Get list of all *_dmrs.tsv files in the nested folder structure
dmr_files <- list.files(
  path = base_dir,
  pattern = "_dmrs.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

# Infer the parameter name from the first file path
#parameter <- str_split(dmr_files[1], "/")[[1]] %>%
#  .[which(. == "dmr_output") + 1]  # grabs folder immediately after 'dmr_output'

# Confirm what was extracted
#cat("Parameter folder detected:", parameter, "\n")



# Function to read and annotate each DMR file
process_dmr_file <- function(file_path) {
  # Remove base_dir prefix and split path
  rel_path <- str_remove(file_path, fixed(paste0(base_dir, "/")))
  path_parts <- unlist(str_split(rel_path, "/"))
  
  # Check expected folder structure: parameter/condition/comparison/filename
  if (length(path_parts) < 3) {
    warning("⚠️ Skipping file due to unexpected folder structure:", file_path)
    return(NULL)
  }
  
  Condition   <- path_parts[1]
  Comparison  <- path_parts[2]
  
  # Read the file
  df <- tryCatch(
    read_tsv(file_path, show_col_types = FALSE),
    error = function(e) {
      warning(paste("Failed to read file:", file_path))
      return(NULL)
    }
  )
  
  if (!is.null(df)) {
    # Expected columns in order
    expected_cols <- c(
      "chr", "start", "end", "num.cpgs1", "num.cpgs2", "cov.score",
      "meth1", "meth2", "diff", "pvalue", "FDR", "annotation"
    )
    
    # Ensure all expected columns exist
    if (!all(expected_cols %in% colnames(df))) {
      warning("⚠️ Skipping file due to missing expected columns: ", file_path)
      return(NULL)
    }
    
    # Subset and reorder columns
    df <- df[, expected_cols]
    
    # Add metadata columns
    df <- df %>%
      mutate(
        condition = Condition,
        comparison = Comparison,
        chromosome_type = ifelse(chr == "Parae_12", "Sex_Ch", "Autosome")
      )
  }
  
  return(df)
}




# Process all files and bind into one data frame
dmr_all <- map_dfr(dmr_files, process_dmr_file)

# Save to file inside the parameter folder
output_path <- file.path(base_dir, "dmr_all.tsv")
write_tsv(dmr_all, output_path)
cat("DMR summary file saved to:", output_path, "\n")


# sanity check of q.value
summary(dmr_all$FDR)
any(dmr_all$FDR > 0.05)  # Should return FALSE

summary(dmr_all$annotation)
