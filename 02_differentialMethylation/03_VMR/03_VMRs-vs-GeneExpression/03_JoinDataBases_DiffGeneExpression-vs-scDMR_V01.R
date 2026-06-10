#!/usr/bin/env Rscript
# DATE: 20260316
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Join cDMR annotations to the differential gene expression table.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# JOIN RULES:
#   - Keep all rows from the differential gene expression table
#   - Join by:
#       DEG$gene_id  <->  DMR$gene
#   - Fill unmatched cDMR information with NA
#
# IMPORTANT:
#   - If one gene is associated with multiple cDMRs, the DEG row will be
#     repeated once per matching cDMR.
#   - If one gene has no matching cDMR, it will still be retained with NA
#     in the cDMR-derived columns.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
DIFF_THRESH <- "0"
plots_folder <- paste0("../plots_", DIFF_THRESH)


dmr_file <- paste0(plots_folder,"/cDMRs_diff", DIFF_THRESH, ".overlap_to_mRNA.tsv")
# dmr_file <- paste0(plots_folder,"/cDMRs_diff", DIFF_THRESH, ".intergenic.tsv")
# dmr_file <- paste0(plots_folder,"/cDMRs_diff", DIFF_THRESH, ".closest_mRNA.tsv")


deg_file <- "Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv"

out_file <- paste0(plots_folder,"/cDMRs_diff", DIFF_THRESH,".overlap_to_mRNA.joined_DEGs.tsv")
# out_file <- paste0(plots_folder,"/cDMRs_diff", DIFF_THRESH,".intergenic.joined_DEGs.tsv")
# out_file <- paste0(plots_folder,"/cDMRs_diff", DIFF_THRESH,".closest_mRNA.joined_DEGs.tsv")




# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------
dmr_dt <- fread(dmr_file, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))
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
  out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote: ", out_file, "\n", sep = "")
cat("Rows in DEG table: ", nrow(deg_dt), "\n", sep = "")
cat("Rows in joined table: ", nrow(joined_dt), "\n", sep = "")