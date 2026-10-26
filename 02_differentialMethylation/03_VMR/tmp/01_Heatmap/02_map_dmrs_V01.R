#!/usr/bin/env Rscript
# DATE: 20260220
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   create a per-DMR table where meth1/meth2 (and cov.score) are mapped into
#   morph-specific columns:
#     meth_female, meth_parae, meth_immaculata, meth_yellow
#     cov.score_female,  cov.score_parae,  cov.score_immaculata,  cov.score_yellow
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Usage:
#   Defaults (if user provides no args):
#     outdir = "."
#     input  = "<outdir>/dmr_real.withIDs.tsv"
#     output = "<outdir>/dmr_real.morph_maped.tsv"
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Args supported:
#   (none)                 -> use defaults above
#   (outdir)               -> use outdir + default input/output file names
#   (input_tsv output_tsv) -> use explicit paths; outdir inferred from output_tsv
#
# Note:
#   - comparison naming: f_vs_i means meth1 corresponds to f, meth2 to i
#   - for real_morph comparisons, remove trailing suffix _## (e.g., i_vs_y_01 -> i_vs_y)
#   - cov.score exists as a single value per DMR; we duplicate it into both participating morphs
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)


default_outdir <- "."
default_infile <- "dmr_real.withIDs.tsv"
default_outfile <- "dmr_real.morph_mapped.tsv"

outdir <- default_outdir
in_tsv <- file.path(outdir, default_infile)
out_tsv <- file.path(outdir, default_outfile)

if (length(args) == 1) {
  # Treat as outdir
  outdir <- args[1]
  in_tsv <- file.path(outdir, default_infile)
  out_tsv <- file.path(outdir, default_outfile)
} else if (length(args) >= 2) {
  # Treat as input and output
  in_tsv <- args[1]
  out_tsv <- args[2]
  outdir <- dirname(out_tsv)
}

# Ensure output directory exists and set working directory as predefined output
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

# After setwd(outdir), rewrite in/out paths to be relative if they are inside outdir
# (keeps behavior consistent whether user passed relative or absolute paths)
in_tsv_abs  <- normalizePath(in_tsv, winslash = "/", mustWork = FALSE)
out_tsv_abs <- normalizePath(out_tsv, winslash = "/", mustWork = FALSE)
outdir_abs  <- normalizePath(outdir, winslash = "/", mustWork = TRUE)

# If input/output are within outdir, use basenames after setwd
if (startsWith(in_tsv_abs, outdir_abs))  in_tsv  <- basename(in_tsv_abs)
if (startsWith(out_tsv_abs, outdir_abs)) out_tsv <- basename(out_tsv_abs)

if (!file.exists(in_tsv)) {
  stop("Input substrate TSV not found: ", file.path(outdir_abs, in_tsv))
}

dt <- fread(in_tsv)

# -----------------------------
# Basic validation
# -----------------------------
required_cols <- c("chr","start","end","cov.score","meth1","meth2","condition","comparison","dmr_id")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Ensure only real_ conditions are present (should be true if you use Script 1 substrate)
bad_cond <- dt[!grepl("^real_", condition), unique(condition)]
if (length(bad_cond) > 0) {
  stop("Found non-real_ conditions in substrate: ", paste(bad_cond, collapse = ", "))
}

# -----------------------------
# Clean comparison names and parse groups
# -----------------------------
# remove trailing _digits (e.g., i_vs_y_01 -> i_vs_y)
dt[, comparison_clean := sub("_[0-9]+$", "", comparison)]

# split on _vs_
parts <- tstrsplit(dt$comparison_clean, "_vs_", fixed = TRUE, keep = c(1,2))
dt[, group1 := parts[[1]]]
dt[, group2 := parts[[2]]]

# mapping of group code -> morph name
code_to_morph <- c(
  f = "female",
  p = "parae",
  i = "immaculata",
  y = "yellow"
)

dt[, morph1 := unname(code_to_morph[group1])]
dt[, morph2 := unname(code_to_morph[group2])]

# QC: verify all groups are recognizable
bad_groups1 <- dt[is.na(morph1), unique(group1)]
bad_groups2 <- dt[is.na(morph2), unique(group2)]
if (length(bad_groups1) > 0 || length(bad_groups2) > 0) {
  cat("Unrecognized group codes detected.\n")
  if (length(bad_groups1) > 0) cat("  bad group1 codes: ", paste(bad_groups1, collapse=", "), "\n", sep="")
  if (length(bad_groups2) > 0) cat("  bad group2 codes: ", paste(bad_groups2, collapse=", "), "\n", sep="")
  stop("Fix comparison parsing/mapping before proceeding.")
}

# -----------------------------
# Build morph-wide meth/cov columns via long->wide
# -----------------------------
# long tables: each DMR contributes two rows (one per participating morph)
dt1 <- dt[, .(dmr_id, morph = morph1, meth = meth1, cov = cov.score)]
dt2 <- dt[, .(dmr_id, morph = morph2, meth = meth2, cov = cov.score)]
long <- rbind(dt1, dt2, use.names = TRUE)

# QC: for each dmr_id, morph should be unique (no duplicates)
dup_pairs <- long[, .N, by = .(dmr_id, morph)][N > 1]
if (nrow(dup_pairs) > 0) {
  cat("ERROR: duplicate (dmr_id, morph) pairs found. First few:\n")
  print(head(dup_pairs, 10))
  stop("This suggests malformed comparisons (e.g., same morph on both sides) or duplicated rows.")
}

# wide casts
meth_w <- dcast(long, dmr_id ~ morph, value.var = "meth")
cov_w  <- dcast(long, dmr_id ~ morph, value.var = "cov")

# ensure all morph columns exist
morphs <- c("female","parae","immaculata","yellow")
for (m in morphs) {
  if (!m %in% names(meth_w)) meth_w[, (m) := as.numeric(NA)]
  if (!m %in% names(cov_w))  cov_w[,  (m) := as.numeric(NA)]
}

# rename to requested headers
setnames(meth_w, morphs, paste0("meth_", morphs))
setnames(cov_w,  morphs, paste0("cov.score_", morphs))

# merge back to main table
out <- merge(dt, meth_w, by = "dmr_id", all.x = TRUE, sort = FALSE)
out <- merge(out, cov_w, by = "dmr_id", all.x = TRUE, sort = FALSE)

# -----------------------------
# QC checks (prints)
# -----------------------------
cat("=== SETTINGS ===\n")
cat("Working directory (outdir): ", getwd(), "\n", sep = "")
cat("Input substrate TSV:        ", normalizePath(in_tsv, winslash="/"), "\n", sep = "")
cat("Output TSV:                 ", normalizePath(out_tsv, winslash="/", mustWork=FALSE), "\n\n", sep = "")

cat("=== QC: rows ===\n")
cat("Input substrate DMR rows:", nrow(dt), "\n")
cat("Output rows:", nrow(out), "\n\n")

cat("=== QC: comparisons (cleaned) ===\n")
print(out[, .N, by = comparison_clean][order(-N)])
cat("\n")

cat("=== QC: each DMR should have exactly 2 non-NA meth_* values ===\n")
meth_cols <- paste0("meth_", morphs)
out[, n_meth_nonNA := rowSums(!is.na(.SD)), .SDcols = meth_cols]
print(out[, .N, by = n_meth_nonNA][order(n_meth_nonNA)])
cat("\n")
if (any(out$n_meth_nonNA != 2)) {
  cat("WARNING: Some DMRs do not have exactly 2 meth_* values.\n")
  cat("First few problematic rows (dmr_id, comparison_clean, n_meth_nonNA):\n")
  print(head(out[n_meth_nonNA != 2, .(dmr_id, comparison_clean, n_meth_nonNA)], 20))
  cat("\n")
}

cat("=== QC: morph mapping sanity (first few) ===\n")
print(head(out[, .(dmr_id, comparison, comparison_clean, group1, group2, morph1, morph2, meth1, meth2)], 10))
cat("\n")

# -----------------------------
# Write output
# -----------------------------
fwrite(out, out_tsv, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote morph-resolved DMR table:\n  ", normalizePath(out_tsv, winslash="/", mustWork=FALSE), "\n", sep = "")

