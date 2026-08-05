#!/usr/bin/env Rscript
# DATE:       2025-09-16
# AUTHOR:     MZF & ChatGPT
# SCRIPT:     02_join_intersect_with_DEGs_V01.R
# VERSION:    01
#
# GOAL:
#   Load the bedtools intersect output (features vs DMRs) and the gene expression
#   database, then merge them using:
#     - intersect column: gene
#     - DE table column:  gene_id
#   Output a compiled table for downstream plotting/correlation.
#
# USAGE:
# Rscript ../scripts/02_join_intersect_with_DEGs_V01.R \
#     intersect_features_vs_dmrs/intersect_features_vs_dmrs.tsv \
#     Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv \
#     compiled_intersect_DEGs_CPMfilter_labeled.tsv
#
# NOTES:
#   - Assumes the intersect file was created by 01_intersect_features_vs_dmrs_V02.sh
#     (so it has feature columns + dmr_* columns).
#   - Keeps all original columns; does an inner join (only genes present in both).
#   - No extra logic added (e.g., mapping comparison directions)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 02_join_intersect_with_DEGs_V01.R <intersect.tsv> <DEGs.tsv> <out.tsv>")
}
intersect_file <- args[1]
deg_file       <- args[2]
out_file       <- args[3]

# --- 1) Read inputs -----------------------------------------------------------

# Intersect output: tab-delimited, header present (from our bash script)
# We won’t guess types too aggressively; let readr infer and keep strings.
x <- readr::read_tsv(intersect_file, show_col_types = FALSE)

# DEGs table: may contain blank cells; treat empty strings as NA.
deg <- readr::read_tsv(
  deg_file,
  na = c("", "NA"),
  show_col_types = FALSE
)

# --- 2) Normalize join keys ---------------------------------------------------

# Ensure the keys exist
if (!("gene" %in% names(x))) {
  stop("Intersect file must contain a 'gene' column (from the features input).")
}
if (!("gene_id" %in% names(deg))) {
  stop("DEGs file must contain a 'gene_id' column.")
}

# Trim whitespace and keep them as plain character
x <- x %>% mutate(gene     = str_trim(as.character(gene)))
deg <- deg %>% mutate(gene_id = str_trim(as.character(gene_id)))

# --- 3) Join (inner) ----------------------------------------------------------

compiled <- x %>%
  inner_join(deg, by = c("gene" = "gene_id"))

# (Optional sanity messages)
message("Rows in intersect: ", nrow(x))
message("Rows in DEGs:      ", nrow(deg))
message("Rows after join:   ", nrow(compiled))

# --- 4) Save ------------------------------------------------------------------

readr::write_tsv(compiled, out_file)
message("Wrote: ", out_file)

# --- 5) Tiny sanity peek (won’t fail if columns absent) ----------------------

if ("dmr_diff" %in% names(compiled)) {
  message("Example: dmr_diff summary:")
  print(summary(compiled$dmr_diff))
}
if ("parae_vs_mel_l2fc" %in% names(compiled)) {
  message("Example: parae_vs_mel_l2fc summary:")
  print(summary(compiled$parae_vs_mel_l2fc))
}

# End
