#!/usr/bin/env Rscript
# DATE: 20260316
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#       Plot the expression per morph given a gene list. 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Run this command to obtain the list of genes 
#
# awk 'NR >1 {print $10}' cDMRs_diff0.overlap_to_mRNA.tsv > names.gene-w-cDMR.txt
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# DESCRIPTION: (Pipeline Step 1/2):
#   Load sample metadata and raw counts, filter samples to Tissue=Muscle &
#   Treatment=veh, compute CPM, flag genes in the subset list, and export
#   tidy tables for downstream stats and plotting.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# INPUT:
#   1) combined_Counts.tsv          (tab, first column "ID" = gene)
#   2) sample_Info.tsv              (tab, columns: Sample, Tissue, Treatment, Morph, ...)
#   3) names.gene-w-UMR-DMV.txt     (one gene ID per line)
# OUTPUT (to --out_dir):
#   - filtered_counts.tsv           (genes x selected samples)
#   - filtered_cpm.tsv              (genes x selected samples, CPM)
#   - tidy_cpm_long.tsv             (long format with metadata, CPM and log2CPM+1)
#   - sample_selection.tsv          (the exact samples used)
# USAGE:
#   Rscript 01_cpm_and_flag.R \
#     --counts combined_Counts.tsv \
#     --samples sample_Info.tsv \
#     --subset names.gene-w-cDMR.txt\
#     --out_dir ./out_cpm_flag
# NOTE:
#   Pure CPM (no gene lengths) is appropriate for 3' RNA-seq. No DESeq2 here.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressPackageStartupMessages({
  # BiocManager::install("optparse")
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

# --------------------------- CLI ---------------------------
DIFF_THRESH <- "0"
plots_folder <- paste0("../plots_", DIFF_THRESH)


opt <- OptionParser()
opt <- add_option(opt, c("--counts"),  type="character", help="Path to combined count table (tsv)", default = "./combined_Counts.tsv")
opt <- add_option(opt, c("--samples"), type="character", help="Path to sample info table (tsv)", default = "./sample_Info.tsv")
opt <- add_option(opt, c("--subset"),  type="character", help="Path to gene list (one per line)", default = "./names.gene-w-cDMR.txt")
# opt <- add_option(opt, c("--out_dir"), type="character", help="Output directory", default = "./out_cpm_flag")
opt <- add_option(opt, c("--out_dir"), type="character", help="Output directory", default = paste0(plots_folder,"/cDMR_geneSubsets_cpm"))
args <- parse_args(opt)

stopifnot(!is.null(args$counts), !is.null(args$samples), !is.null(args$subset))

out_dir <- args$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------- Load inputs ------------------------
message("[i] Reading sample info: ", args$samples)
sinfo <- read_tsv(args$samples, show_col_types = FALSE)

required_cols <- c("Sample","Tissue","Treatment","Morph")
missing_cols <- setdiff(required_cols, names(sinfo))
if (length(missing_cols) > 0) {
  stop("sample_Info.tsv is missing required columns: ", paste(missing_cols, collapse=", "))
}

# Filter to Tissue=Muscle & Treatment=veh
s_keep <- sinfo %>%
  filter(Tissue == "Muscle", Treatment == "veh") %>%
  distinct(Sample, .keep_all = TRUE)

if (nrow(s_keep) == 0) stop("No samples matched Tissue=Muscle & Treatment=veh")

# Load counts
message("[i] Reading counts: ", args$counts)
counts <- read_tsv(args$counts, show_col_types = FALSE)

if (!"ID" %in% names(counts)) stop("Counts file must have a column named 'ID' with gene IDs")

# Intersect columns with selected samples
sample_cols <- intersect(s_keep$Sample, setdiff(names(counts), "ID"))
if (length(sample_cols) == 0) {
  stop("None of the selected samples are present as columns in counts file.")
}

# Keep only ID + selected samples
counts_f <- counts %>% select(ID, all_of(sample_cols))

# Ensure numeric matrix (handle any non-numeric)
counts_mat <- counts_f %>%
  mutate(across(-ID, ~ suppressWarnings(as.numeric(.)))) 

# Load subset gene list
message("[i] Reading subset gene list: ", args$subset)
subset_genes <- read_tsv(args$subset, col_names = FALSE, show_col_types = FALSE)[[1]] %>% unique()

# ---------------------- Compute CPM ------------------------
# CPM_j(g) = 1e6 * counts_j(g) / sum_g counts_j(g)
lib_sizes <- colSums(as.matrix(counts_mat %>% select(-ID)), na.rm = TRUE)

zero_lib <- names(lib_sizes)[lib_sizes == 0]
if (length(zero_lib) > 0) {
  warning("These samples have zero total counts and will be dropped: ", paste(zero_lib, collapse=", "))
  keep_cols <- setdiff(names(lib_sizes), zero_lib)
  counts_mat <- counts_mat %>% select(ID, all_of(keep_cols))
  lib_sizes <- lib_sizes[keep_cols]
}

cpm <- counts_mat
for (nm in names(lib_sizes)) {
  cpm[[nm]] <- 1e6 * (counts_mat[[nm]] / lib_sizes[[nm]])
}

# ------------------ Export wide summaries ------------------
write_tsv(counts_mat, file.path(out_dir, "filtered_counts.tsv"))
write_tsv(cpm,        file.path(out_dir, "filtered_cpm.tsv"))
write_tsv(s_keep %>% filter(Sample %in% names(lib_sizes)),
          file.path(out_dir, "sample_selection.tsv"))

# ------------------ Build tidy long table ------------------
long <- cpm %>%
  pivot_longer(-ID, names_to = "Sample", values_to = "CPM") %>%
  left_join(s_keep %>% select(Sample, Morph, Tissue, Treatment), by = "Sample") %>%
  mutate(logCPM1 = log2(CPM + 1),
         InSubset = ID %in% subset_genes)

write_tsv(long, file.path(out_dir, "tidy_cpm_long.tsv"))

message("[✓] Done. Outputs written to: ", out_dir)
