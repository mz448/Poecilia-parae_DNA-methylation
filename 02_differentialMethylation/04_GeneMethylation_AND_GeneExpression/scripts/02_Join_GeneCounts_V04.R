#!/usr/bin/env Rscript
# DATE: 20260317
# AUTHOR: MZF
# SCRIPT: 04_Join_GeneCounts_V04.R
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Build a reusable joined table linking morph-level methylation and
#   morph-level gene expression by gene and morph.
#
# DESCRIPTION:
#   This script reads:
#     1) per_morph_gene_means.tsv
#        - long format
#        - one row per gene x morph
#
#     2) meth.combined.tsv
#        - wide format
#        - one row per gene x feature
#        - morph-level methylation stored in separate columns:
#            meth.morph.female
#            meth.morph.immaculata
#            meth.morph.parae
#            meth.morph.yellow
#
#   It then:
#     A) standardizes gene identifiers
#     B) reshapes methylation wide -> long
#     C) standardizes morph names across both tables
#     D) joins the tables by:
#           gene_id + Morph_join
#     E) writes a reusable joined table for downstream analyses
#
# MORPH MATCHING RULES:
#   per_morph_gene_means.tsv   meth.combined.tsv   Morph_join
#   female                     female              female
#   immac                      immaculata          immaculata
#   mel                        yellow              yellow
#   parae                      parae               parae
#
# JOIN RULES:
#   - Keep all rows from the methylation table after reshaping to long format
#   - Join CPM / expression information onto those rows
#   - One output row corresponds to:
#         gene x feature x morph
#
# EXPECTED INPUT FILES:
#   - ./scDMR_geneSubsets_cpm/per_morph_gene_means.tsv
#   - ./methylation_means/final_tables/meth.combined.tsv
#
# MAIN OUTPUT FILES:
#   - ./scDMR_geneSubsets_cpm/methPerFeature_AND_cpm.long.tsv
#       Reusable joined long table
#
#   - ./scDMR_geneSubsets_cpm/methPerFeature_AND_cpm.morph_check.tsv
#       Morph-matching diagnostic summary
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
cpm_file  <- "./scDMR_geneSubsets_cpm/per_morph_gene_means.tsv"
meth_file <- "./methylation_means/final_tables/meth.combined.tsv"

# ------------------------------------------------------------
# Outputs
# ------------------------------------------------------------
out_long_file   <- "./scDMR_geneSubsets_cpm/methPerFeature_AND_cpm.long.tsv"
out_morph_check <- "./scDMR_geneSubsets_cpm/methPerFeature_AND_cpm.morph_check.tsv"

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------
meth_dt <- fread(meth_file, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))
cpm_dt  <- fread(cpm_file,  sep = "\t", header = TRUE, na.strings = c("NA", "", "."))

cat("Read methylation rows: ", nrow(meth_dt), "\n", sep = "")
cat("Read CPM rows: ", nrow(cpm_dt), "\n", sep = "")

# ------------------------------------------------------------
# Basic input checks
# ------------------------------------------------------------
if (!"gene" %in% names(meth_dt)) stop("Missing column 'gene' in methylation table.")
if (!"feature" %in% names(meth_dt)) stop("Missing column 'feature' in methylation table.")

if (!"ID" %in% names(cpm_dt)) stop("Missing column 'ID' in CPM table.")
if (!"Morph" %in% names(cpm_dt)) stop("Missing column 'Morph' in CPM table.")

# ------------------------------------------------------------
# Standardize gene column names
# ------------------------------------------------------------
setnames(meth_dt, "gene", "gene_id")
setnames(cpm_dt,  "ID",   "gene_id")

# ------------------------------------------------------------
# Identify methylation columns
# ------------------------------------------------------------
base_cols <- c("chr", "start", "end", "gene_id", "feature", "chromosome_type")

missing_base <- setdiff(base_cols, names(meth_dt))
if (length(missing_base) > 0) {
  stop(
    "Missing expected columns in methylation table: ",
    paste(missing_base, collapse = ", ")
  )
}

meth_cols <- grep("^meth\\.morph\\.", names(meth_dt), value = TRUE)
cov_cols  <- grep("^cov\\.morph\\.",  names(meth_dt), value = TRUE)

if (length(meth_cols) == 0) stop("No columns matched ^meth\\.morph\\.")
if (length(cov_cols)  == 0) stop("No columns matched ^cov\\.morph\\.")

meth_morphs <- sub("^meth\\.morph\\.", "", meth_cols)
cov_morphs  <- sub("^cov\\.morph\\.",  "", cov_cols)

if (!setequal(meth_morphs, cov_morphs)) {
  stop("Morph names in meth.morph.* and cov.morph.* do not match.")
}

# preserve the stable order found in methylation columns
morph_levels_meth <- meth_morphs

# ------------------------------------------------------------
# Reduce methylation table to a manageable set of columns
# ------------------------------------------------------------
meth_keep  <- unique(c(base_cols, cov_cols, meth_cols))
meth_small <- meth_dt[, ..meth_keep]

# ------------------------------------------------------------
# Reshape methylation wide -> long
# Output: one row per gene x feature x morph
# ------------------------------------------------------------
meth_long_list <- lapply(morph_levels_meth, function(m) {
  cov_col  <- paste0("cov.morph.", m)
  meth_col <- paste0("meth.morph.", m)
  
  meth_small[, .(
    gene_id         = gene_id,
    feature         = feature,
    chr             = chr,
    start           = start,
    end             = end,
    chromosome_type = chromosome_type,
    Morph_meth      = m,
    cov_morph       = get(cov_col),
    meth_morph      = get(meth_col)
  )]
})

meth_long <- rbindlist(meth_long_list, use.names = TRUE)

# ------------------------------------------------------------
# Prepare CPM table
# ------------------------------------------------------------
cpm_keep <- c("gene_id", "Morph", "InSubset", "mean_CPM", "mean_logCPM1", "n_rep")
missing_cpm <- setdiff(cpm_keep, names(cpm_dt))

if (length(missing_cpm) > 0) {
  stop("Missing expected CPM columns: ", paste(missing_cpm, collapse = ", "))
}

cpm_small <- cpm_dt[, ..cpm_keep]
setnames(cpm_small, "Morph", "Morph_cpm")

# ------------------------------------------------------------
# Standardize morph names for joining
# ------------------------------------------------------------
map_cpm <- c(
  "female" = "female",
  "immac"  = "immaculata",
  "mel"    = "yellow",
  "parae"  = "parae"
)

map_meth <- c(
  "female"      = "female",
  "immaculata"  = "immaculata",
  "yellow"      = "yellow",
  "parae"       = "parae"
)

cpm_small[, Morph_join := unname(map_cpm[Morph_cpm])]
meth_long[, Morph_join := unname(map_meth[Morph_meth])]

# ------------------------------------------------------------
# Safety checks for morph mapping
# ------------------------------------------------------------
if (any(is.na(cpm_small$Morph_join))) {
  bad <- unique(cpm_small[is.na(Morph_join), Morph_cpm])
  stop("Unmatched morph names in CPM table: ", paste(bad, collapse = ", "))
}

if (any(is.na(meth_long$Morph_join))) {
  bad <- unique(meth_long[is.na(Morph_join), Morph_meth])
  stop("Unmatched morph names in methylation table: ", paste(bad, collapse = ", "))
}

# ------------------------------------------------------------
# Join by gene + standardized morph
# Keep all methylation rows
# ------------------------------------------------------------
joined_dt <- merge(
  meth_long,
  cpm_small,
  by = c("gene_id", "Morph_join"),
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------
# Reorder columns
# ------------------------------------------------------------
wanted_order <- c(
  "gene_id", "Morph_join", "Morph_meth", "Morph_cpm",
  "feature", "chr", "start", "end", "chromosome_type",
  "cov_morph", "meth_morph",
  "InSubset", "mean_CPM", "mean_logCPM1", "n_rep"
)

present_order <- wanted_order[wanted_order %in% names(joined_dt)]
setcolorder(joined_dt, present_order)

# ------------------------------------------------------------
# Write reusable joined long table
# ------------------------------------------------------------
fwrite(
  joined_dt,
  out_long_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote full joined long table: ", out_long_file, "\n", sep = "")
cat("Rows in methylation-long table: ", nrow(meth_long), "\n", sep = "")
cat("Rows in joined table: ", nrow(joined_dt), "\n", sep = "")

# ------------------------------------------------------------
# Write morph matching diagnostic
# ------------------------------------------------------------
morph_check_dt <- joined_dt[
  , .N,
  by = .(Morph_join, Morph_meth, Morph_cpm)
][order(Morph_join, Morph_meth, Morph_cpm)]

fwrite(
  morph_check_dt,
  out_morph_check,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote morph matching diagnostic table: ", out_morph_check, "\n", sep = "")
cat("Done.\n")