#!/usr/bin/env Rscript
# DATE: 20260316
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Join VMR annotations to the differential gene expression table.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# JOIN RULES:
#   - Keep all rows from the differential gene expression table
#   - Join by:
#       DEG$gene_id  <->  DMR$gene
#   - Fill unmatched VMR information with NA
#
# IMPORTANT:
#   - If one gene is associated with multiple VMRs, the DEG row will be
#     repeated once per matching VMR.
#   - If one gene has no matching VMR, it will still be retained with NA
#     in the VMR-derived columns.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# V04: Fixes name variability in input files. 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

# SELECT ONE: to join to the DEG database
  FILE <- "closest_mRNA"      # all the VMRs 
  # FILE <- "overlap_to_mRNA" # Only the VMRs that overlap with Genes
  # FILE <- "intergenic"      # Only the VMRs that do NOT overlap with Genes
  

DIFF_THRESH <- "0"
plots_folder <- paste0("../plots")
vmr_file <- paste0("VMR_diff", DIFF_THRESH,".",FILE,".tsv")
out_file <- paste0("VMR_diff", DIFF_THRESH,".",FILE,".joined-DEG.tsv")
deg_file <- "Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv"

# dir_create(out_file, recurse = TRUE)

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------
dmr_dt <- fread(vmr_file, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))
deg_dt <- fread(deg_file, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))

# ------------------------------------------------------------
# Prepare DMR table
# ------------------------------------------------------------
# keep all DMR columns except the helper column x, if present
if ("x" %in% names(dmr_dt)) {
  dmr_dt[, x := NULL]
}

# rename gene -> gene_id so the join key matches DEG
setnames(dmr_dt, "gene", "gene_id")

# ------------------------------------------------------------
# Left join: keep all DEG rows
# ------------------------------------------------------------
joined_dt <- merge(
  deg_dt,
  dmr_dt,
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------
fwrite(
  joined_dt,
  file=out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote: ", out_file, "\n", sep = "")
cat("Rows in DEG table: ", nrow(deg_dt), "\n", sep = "")
cat("Rows in joined table: ", nrow(joined_dt), "\n", sep = "")