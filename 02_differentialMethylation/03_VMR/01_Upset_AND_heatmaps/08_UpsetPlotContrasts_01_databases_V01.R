#!/usr/bin/env Rscript
# DATE: 20260401
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Build a consensus-DMR x contrast membership table for UpSet plotting.
#
#   INPUTS:
#     1) dmr_real_consensus.bed
#        cols: chr, start, end, cons_id, dmr_ids
#     2) dmr_real.withIDs.tsv
#        contains: dmr_id, comparison, condition, ...
#
#   OUTPUTS:
#     A) long table:
#        one row per cons_id x contributing DMR x contributing contrast
#     B) binary wide table:
#        one row per cons_id, one 0/1 column per contrast
#
#   LOGIC:
#     - split collapsed dmr_ids from the consensus BED
#     - join each dmr_id back to dmr_real.withIDs.tsv
#     - recover the contrast(s) contributing to each consensus DMR
#     - collapse to unique cons_id x contrast combinations
#     - cast to wide binary matrix for UpSet plotting
#
#   NOTE:
#     By default, trailing replicate suffixes are removed:
#       i_vs_y_01 -> i_vs_y
#     This is usually what you want for a contrast-sharing UpSet.
#     If you want the original names instead, replace comparison_clean with
#     comparison in the marked lines below.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

# -----------------------------
# Defaults
# -----------------------------
default_outdir    <- "."
default_cons_bed  <- "dmr_real_consensus.bed"
default_dmr_tsv   <- "dmr_real.withIDs.tsv"
default_prefix    <- "dmr_real_consensus.contrastSharing"

outdir   <- default_outdir
cons_bed <- file.path(outdir, default_cons_bed)
dmr_tsv  <- file.path(outdir, default_dmr_tsv)
prefix   <- file.path(outdir, default_prefix)

if (length(args) == 1) {
  outdir   <- args[1]
  cons_bed <- file.path(outdir, default_cons_bed)
  dmr_tsv  <- file.path(outdir, default_dmr_tsv)
  prefix   <- file.path(outdir, default_prefix)
} else if (length(args) >= 3) {
  cons_bed <- args[1]
  dmr_tsv  <- args[2]
  prefix   <- args[3]
  outdir   <- dirname(prefix)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

out_long   <- paste0(prefix, ".long.tsv")
out_binary <- paste0(prefix, ".binary.tsv")
out_counts <- paste0(prefix, ".contrastCounts.tsv")

# -----------------------------
# Read inputs
# -----------------------------
cons <- fread(cons_bed, header = FALSE)
if (ncol(cons) < 5) {
  stop("Consensus BED must have at least 5 columns: chr, start, end, cons_id, dmr_ids")
}
setnames(cons, 1:5, c("chr", "start", "end", "cons_id", "dmr_ids"))
cons[, start := as.integer(start)]
cons[, end   := as.integer(end)]

dmr <- fread(dmr_tsv)

required_cols <- c("dmr_id", "comparison")
missing_cols <- setdiff(required_cols, names(dmr))
if (length(missing_cols) > 0) {
  stop("Missing required columns in dmr_real.withIDs.tsv: ",
       paste(missing_cols, collapse = ", "))
}

if ("condition" %in% names(dmr)) {
  bad_cond <- dmr[!grepl("^real_", condition), unique(condition)]
  if (length(bad_cond) > 0) {
    warning("Found non-real_ conditions in dmr table: ",
            paste(bad_cond, collapse = ", "),
            "\nThe script will keep all rows present in the file you supplied.")
  }
}

# -----------------------------
# Clean comparison names
# -----------------------------
dmr[, comparison_clean := sub("_[0-9]+$", "", comparison)]

# keep one row per dmr_id
dmr_map <- unique(dmr[, .(dmr_id, comparison, comparison_clean)])

# -----------------------------
# Expand consensus mapping:
# one row per cons_id x dmr_id
# -----------------------------
map <- cons[, .(
  dmr_id = unlist(strsplit(dmr_ids, ",", fixed = TRUE))
), by = .(chr, start, end, cons_id, dmr_ids)]

expected_pairs <- cons[, sum(lengths(strsplit(dmr_ids, ",", fixed = TRUE)))]
if (nrow(map) != expected_pairs) {
  stop("Mapping expansion mismatch: expected ", expected_pairs,
       " cons_id-to-dmr_id pairs but got ", nrow(map))
}

# -----------------------------
# Join dmr_id back to contrast info
# -----------------------------
joined <- merge(
  map,
  dmr_map,
  by = "dmr_id",
  all.x = TRUE,
  sort = FALSE
)

missing_ids <- joined[is.na(comparison), unique(dmr_id)]
if (length(missing_ids) > 0) {
  cat("ERROR: some dmr_id values from consensus BED were not found in dmr_real.withIDs.tsv\n")
  cat("First missing IDs:\n")
  print(head(missing_ids, 20))
  stop("Fix the consistency between dmr_real_consensus.bed and dmr_real.withIDs.tsv")
}

# derive chromosome type from consensus coordinates
joined[, chromosome_type := ifelse(chr == "Parae_12", "Sex_Ch", "Autosome")]

# -----------------------------
# Long table:
# one row per unique cons_id x contrast
# -----------------------------
contrast_long <- unique(
  joined[, .(
    chr,
    start,
    end,
    cons_id,
    dmr_ids,
    chromosome_type,
    dmr_id,
    comparison,
    comparison_clean
  )]
)

# -----------------------------
# Binary table for UpSet
# use cleaned contrast names by default
# If you want suffix-aware columns, change comparison_clean -> comparison
# in the two marked places below.
# -----------------------------
binary_long <- unique(
  contrast_long[, .(
    chr,
    start,
    end,
    cons_id,
    dmr_ids,
    chromosome_type,
    contrast = comparison_clean,   # <----- use comparison here if desired
    present = 1L
  )]
)

wide <- dcast(
  binary_long,
  chr + start + end + cons_id + dmr_ids + chromosome_type ~ contrast,
  value.var = "present",
  fill = 0
)

# -----------------------------
# Order contrast columns
# -----------------------------
meta_cols <- c("chr", "start", "end", "cons_id", "dmr_ids", "chromosome_type")
contrast_cols <- setdiff(names(wide), meta_cols)

preferred_order <- c("f_vs_i", "f_vs_p", "f_vs_y", "i_vs_y", "p_vs_i", "y_vs_p")
contrast_cols <- c(
  intersect(preferred_order, contrast_cols),
  setdiff(sort(contrast_cols), preferred_order)
)

setcolorder(wide, c(meta_cols, contrast_cols))

# -----------------------------
# Add useful summary columns
# -----------------------------
wide[, n_contrasts := rowSums(.SD), .SDcols = contrast_cols]

wide[, contrast_list := apply(.SD, 1, function(x) {
  hits <- contrast_cols[as.logical(x)]
  if (length(hits) == 0) return("none")
  paste(hits, collapse = ",")
}), .SDcols = contrast_cols]

wide[, contrast_combo := apply(.SD, 1, function(x) {
  hits <- contrast_cols[as.logical(x)]
  if (length(hits) == 0) return("none")
  paste(hits, collapse = "&")
}), .SDcols = contrast_cols]

# move summary columns next to metadata
setcolorder(
  wide,
  c("chr", "start", "end", "cons_id", "dmr_ids", "chromosome_type",
    "n_contrasts", "contrast_list", "contrast_combo",
    contrast_cols)
)

# -----------------------------
# Simple counts per contrast
# -----------------------------
contrast_counts <- data.table(
  contrast = contrast_cols,
  n_consensus_dmrs = vapply(contrast_cols, function(cc) sum(wide[[cc]] == 1L), numeric(1))
)

# -----------------------------
# QC
# -----------------------------
cat("=== SETTINGS ===\n")
cat("Consensus BED:   ", normalizePath(cons_bed, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("DMR TSV:         ", normalizePath(dmr_tsv,  winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Output prefix:   ", normalizePath(prefix,   winslash = "/", mustWork = FALSE), "\n\n", sep = "")

cat("=== QC: counts ===\n")
cat("Consensus regions:                 ", nrow(cons), "\n", sep = "")
cat("Expanded cons_id -> dmr_id pairs:  ", nrow(map), "\n", sep = "")
cat("Unique cons_id x contrast rows:    ", nrow(binary_long), "\n", sep = "")
cat("Binary matrix rows:                ", nrow(wide), "\n\n", sep = "")

cat("=== QC: contrasts recovered ===\n")
print(contrast_counts[order(-n_consensus_dmrs)])
cat("\n")

cat("=== QC: number of contrasts per consensus region ===\n")
print(wide[, .N, by = n_contrasts][order(n_contrasts)])
cat("\n")

cat("=== QC: chromosome breakdown ===\n")
print(wide[, .N, by = chromosome_type][order(-N)])
cat("\n")

stopifnot(anyDuplicated(wide$cons_id) == 0)

# -----------------------------
# Write outputs
# -----------------------------
fwrite(contrast_long,   out_long,   sep = "\t", quote = FALSE, na = "NA")
fwrite(wide,            out_binary, sep = "\t", quote = FALSE, na = "NA")
fwrite(contrast_counts, out_counts, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote:\n")
cat("  ", normalizePath(out_long,   winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(out_binary, winslash = "/", mustWork = FALSE), "\n", sep = "")
cat("  ", normalizePath(out_counts, winslash = "/", mustWork = FALSE), "\n", sep = "")
