#!/usr/bin/env Rscript
# DATE:       2026-03-23
# AUTHOR:     MZF & ChatGPT
# SCRIPT:     07_flag_consensus_DMR_classes_V01.R
# VERSION:    01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Classify consensus DMRs from a morph-mapped consensus table according to
#   which groups contribute non-NA methylation values in the morph-specific
#   meth_* columns.
#
# INPUT:
#   A TSV file like:
#     dmr_real_consensus.morph_mapped.tsv
#
#   Required methylation columns:
#     meth_female
#     meth_immaculata
#     meth_parae
#     meth_yellow
#
# LOGIC:
#   Presence is defined as meth_* != NA.
#
#   New flags added:
#     - has_female
#     - has_immaculata
#     - has_parae
#     - has_yellow
#     - n_male_morphs_present
#     - n_groups_present
#     - sex_cDMR
#     - morph_cDMR
#     - sex_cDMR_1   = female + exactly 1 male morph
#     - sex_cDMR_2   = female + exactly 2 male morphs
#     - sex_cDMR_3   = female + exactly 3 male morphs
#     - morph_cDMR_2 = exactly 2 male morphs and no female
#     - morph_cDMR_3 = exactly 3 male morphs and no female
#
#   Definitions used:
#     sex_cDMR   = female present AND at least 1 male morph present
#     morph_cDMR = female absent  AND at least 2 male morphs present
#
# OUTPUTS:
#   1) Full table with all original columns + new flag columns
#   2) sex_cDMR-only table
#   3) morph_cDMR-only table
#
# USAGE:
#   Rscript 01_flag_consensus_DMR_classes.R \
#     dmr_real_consensus.morph_mapped.tsv \
#     results_dir
#
# NOTES:
#   - This script keeps all original columns in every output.
#   - The specialized classes are nested within the broader classes:
#       sex_cDMR_1 / _2 / _3 are subsets of sex_cDMR
#       morph_cDMR_2 is a subset of morph_cDMR
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# 1) LOAD INPUT
# =========================

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) args[1] else "dmr_real_consensus.morph_mapped.tsv"
outdir     <- if (length(args) >= 2) args[2] else "classes_consensus_DMRs"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dt <- fread(input_file, sep = "\t", header = TRUE, na.strings = c("NA", "", "NaN"))

required_cols <- c(
  "meth_female",
  "meth_immaculata",
  "meth_parae",
  "meth_yellow"
)

missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# =========================
# 2) PROCESS DATAFRAME
# =========================

# Presence flags based on non-NA methylation values
dt[, has_female      := !is.na(meth_female)]
dt[, has_immaculata  := !is.na(meth_immaculata)]
dt[, has_parae       := !is.na(meth_parae)]
dt[, has_yellow      := !is.na(meth_yellow)]

# Number of male morphs represented
dt[, n_male_morphs_present := 
     as.integer(has_immaculata) +
     as.integer(has_parae) +
     as.integer(has_yellow)]

# Total number of groups represented
dt[, n_groups_present :=
     as.integer(has_female) +
     as.integer(has_immaculata) +
     as.integer(has_parae) +
     as.integer(has_yellow)]

# Broad classes
dt[, sex_cDMR   := has_female & (n_male_morphs_present >= 1)]
dt[, morph_cDMR := (!has_female) & (n_male_morphs_present >= 2)]

# Specific nested classes
dt[, sex_cDMR_1   := has_female & (n_male_morphs_present == 1)]
dt[, sex_cDMR_2   := has_female & (n_male_morphs_present == 2)]
dt[, sex_cDMR_3   := has_female & (n_male_morphs_present == 3)]
dt[, morph_cDMR_2 := (!has_female) & (n_male_morphs_present == 2)]
dt[, morph_cDMR_3 := (!has_female) & (n_male_morphs_present == 3)]

# Optional stable ordering: keep original columns first, append new columns last
new_flag_cols <- c(
  "has_female",
  "has_immaculata",
  "has_parae",
  "has_yellow",
  "n_male_morphs_present",
  "n_groups_present",
  "sex_cDMR",
  "morph_cDMR",
  "sex_cDMR_1",
  "sex_cDMR_2",
  "sex_cDMR_3",
  "morph_cDMR_2",
  "morph_cDMR_3"
)

original_cols <- setdiff(names(dt), new_flag_cols)
setcolorder(dt, c(original_cols, new_flag_cols))

# =========================
# 3) SAVE OUTPUTS
# =========================

base_name <- tools::file_path_sans_ext(basename(input_file))

full_out  <- file.path(outdir, paste0(base_name, ".flagged.tsv"))
sex_out   <- file.path(outdir, paste0(base_name, ".sex_cDMR.tsv"))
morph_out <- file.path(outdir, paste0(base_name, ".morph_cDMR.tsv"))

fwrite(dt, full_out, sep = "\t", quote = FALSE, na = "NA")

fwrite(
  dt[sex_cDMR == TRUE],
  sex_out,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  dt[morph_cDMR == TRUE],
  morph_out,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Done.\n")
cat("Full flagged table:  ", full_out,  "\n", sep = "")
cat("sex_cDMR table:      ", sex_out,   "\n", sep = "")
cat("morph_cDMR table:    ", morph_out, "\n", sep = "")