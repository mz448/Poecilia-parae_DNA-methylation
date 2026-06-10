#!/usr/bin/env Rscript
# DATE: 20260220
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL: 
#     Join consensus regions to the set of contributing DMRs, add the meth value  
#     and compute one coverage value per morph per consensus region.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# DETAILS: 
#     (one row per consensus DMR) with morph-specific meth/cov values:
#     chr, start, end,
#     cov.score_female, cov.score_immaculata, cov.score_parae, cov.score_yellow,
#     meth_female,     meth_immaculata,     meth_parae,     meth_yellow
#
# Inputs (defaults):
#   outdir = "./dmr_consensus_real"
#   consensus BED from Script 1: "<outdir>/consensus_real_dmrs.bed"
#     columns: chr, start, end, cons_id, collapsed_dmr_ids
#   morph-resolved DMRs from Script 2: "<outdir>/real_dmrs.morph_resolved.tsv"
#
# Output (default):
#   "<outdir>/consensus_master_matrix.real.tsv"
#
# Aggregation rule (per consensus region, per morph):
#   median of all non-NA values contributed by the constituent DMRs.
#
# Optional args:
#   (none) -> use defaults
#   (outdir) -> use outdir + default filenames
#   (consensus_bed morph_resolved_tsv output_tsv) -> explicit paths; outdir inferred from output_tsv
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

default_outdir <- "."
default_consensus <- "dmr_real_consensus.bed"
default_morphs  <- "dmr_real.morph_mapped.tsv"
default_outfile   <- "dmr_real_consensus.morph_mapped.tsv"

outdir <- default_outdir
cons_bed <- file.path(outdir, default_consensus)
mr_tsv   <- file.path(outdir, default_morphs)
out_tsv  <- file.path(outdir, default_outfile)

if (length(args) == 1) {
  outdir <- args[1]
  cons_bed <- file.path(outdir, default_consensus)
  mr_tsv   <- file.path(outdir, default_morphs)
  out_tsv  <- file.path(outdir, default_outfile)
} else if (length(args) >= 3) {
  cons_bed <- args[1]
  mr_tsv   <- args[2]
  out_tsv  <- args[3]
  outdir   <- dirname(out_tsv)
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
setwd(outdir)

# Normalize paths: after setwd, use basenames if inside outdir
cons_bed_abs <- normalizePath(cons_bed, winslash="/", mustWork=FALSE)
mr_tsv_abs   <- normalizePath(mr_tsv, winslash="/", mustWork=FALSE)
out_tsv_abs  <- normalizePath(out_tsv, winslash="/", mustWork=FALSE)
outdir_abs   <- normalizePath(outdir, winslash="/", mustWork=TRUE)

if (startsWith(cons_bed_abs, outdir_abs)) cons_bed <- basename(cons_bed_abs)
if (startsWith(mr_tsv_abs, outdir_abs))   mr_tsv   <- basename(mr_tsv_abs)
if (startsWith(out_tsv_abs, outdir_abs))  out_tsv  <- basename(out_tsv_abs)

if (!file.exists(cons_bed)) stop("Consensus BED not found: ", file.path(outdir_abs, cons_bed))
if (!file.exists(mr_tsv))   stop("Morph-resolved TSV not found: ", file.path(outdir_abs, mr_tsv))

# -----------------------------
# Read inputs
# -----------------------------
cons <- fread(cons_bed, header = FALSE)
if (ncol(cons) < 5) stop("Consensus BED must have >=5 columns: chr start end cons_id collapsed_dmr_ids")

setnames(cons, 1:5, c("chr","start","end","cons_id","dmr_ids"))
cons[, start := as.integer(start)]
cons[, end   := as.integer(end)]

mr <- fread(mr_tsv)

required_mr <- c("dmr_id",
                 "meth_female","meth_parae","meth_immaculata","meth_yellow",
                 "cov.score_female","cov.score_parae","cov.score_immaculata","cov.score_yellow")
missing_mr <- setdiff(required_mr, names(mr))
if (length(missing_mr) > 0) stop("Missing required columns in morph-resolved TSV: ", paste(missing_mr, collapse=", "))

# -----------------------------
# Expand consensus mapping: cons_id -> dmr_id
# -----------------------------
map <- cons[, .(dmr_id = unlist(strsplit(dmr_ids, ",", fixed = TRUE))), by = .(cons_id)]
# Keep the coordinates attached for later
map <- merge(map, cons[, .(cons_id, chr, start, end, dmr_ids)], by = "cons_id", all.x = TRUE, sort = FALSE)

# QC: mapping size should equal sum of splits
expected_pairs <- cons[, sum(lengths(strsplit(dmr_ids, ",", fixed = TRUE)))]
if (nrow(map) != expected_pairs) {
  stop("Mapping expansion mismatch: expected ", expected_pairs, " pairs but got ", nrow(map))
}

# -----------------------------
# Join morph-resolved data onto mapping
# -----------------------------
setkey(mr, dmr_id)
setkey(map, dmr_id)
joined <- mr[map, on = "dmr_id"]

# QC: ensure no missing joins (dmr_id not found)
n_missing <- joined[is.na(dmr_id) | (is.na(meth_female) & is.na(meth_parae) & is.na(meth_immaculata) & is.na(meth_yellow)), .N]
# The check above is conservative; better: check for NA on a column that must exist post-join
n_unmatched <- joined[is.na(comparison) & is.na(condition) & is.na(chr), .N]  # fallback if those cols exist
# Safer: check which map dmr_ids didn't match
unmatched <- joined[is.na(mr$dmr_id), unique(dmr_id)]  # won't work as intended; do explicit:
unmatched_ids <- setdiff(map$dmr_id, mr$dmr_id)

if (length(unmatched_ids) > 0) {
  cat("ERROR: Some dmr_id(s) from consensus BED not found in morph-resolved TSV.\n")
  cat("First few missing IDs: ", paste(head(unmatched_ids, 20), collapse=", "), "\n", sep="")
  stop("Fix ID consistency between Script 1 and Script 2 outputs.")
}


# -----------------------------
# Aggregate to consensus x morph (median of non-NA values)
# -----------------------------
morphs <- c("female","immaculata","parae","yellow")

meth_cols <- paste0("meth_", morphs)
cov_cols  <- paste0("cov.score_", morphs)

# helper median that returns NA if all missing
med_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

agg <- joined[, {
  # counts for QC/diagnostics
  n_dmrs <- .N
  list(
    n_dmrs = n_dmrs,
    dmr_ids = unique(dmr_ids)[1],  # same for all rows in group
    meth_female      = med_na(meth_female),
    meth_immaculata  = med_na(meth_immaculata),
    meth_parae       = med_na(meth_parae),
    meth_yellow      = med_na(meth_yellow),
    cov.score_female     = med_na(cov.score_female),
    cov.score_immaculata = med_na(cov.score_immaculata),
    cov.score_parae      = med_na(cov.score_parae),
    cov.score_yellow     = med_na(cov.score_yellow)
  )
}, by = .(cons_id, chr, start, end)]

# Order columns as requested (chr/start/end first, then covs, then meths)
setcolorder(
  agg,
  c("chr","start","end",
    "cov.score_female","cov.score_immaculata","cov.score_parae","cov.score_yellow",
    "meth_female","meth_immaculata","meth_parae","meth_yellow",
    "cons_id","n_dmrs","dmr_ids")
)

# -----------------------------
# Add chromosome_type (optional but useful; derived from chr)
# -----------------------------
agg[, chromosome_type := ifelse(grepl("^Parae_12$", chr), "Sex_Ch", "Autosome")]
# If your sex chromosome has other naming variants, adjust the regex accordingly.




# -----------------------------
# QC summaries
# -----------------------------
cat("=== SETTINGS ===\n")
cat("Working directory (outdir): ", getwd(), "\n", sep = "")
cat("Consensus BED:              ", normalizePath(cons_bed, winslash="/"), "\n", sep = "")
cat("Morph-resolved TSV:         ", normalizePath(mr_tsv, winslash="/"), "\n", sep = "")
cat("Output master TSV:          ", normalizePath(out_tsv, winslash="/", mustWork=FALSE), "\n\n", sep = "")

cat("=== QC: counts ===\n")
cat("Consensus regions: ", nrow(cons), "\n", sep = "")
cat("Mapping pairs (cons_id->dmr_id): ", nrow(map), "\n", sep = "")
cat("Master rows: ", nrow(agg), "\n\n", sep = "")

cat("=== QC: n_dmrs per consensus region (how many original DMRs merged) ===\n")
print(agg[, .N, by = n_dmrs][order(n_dmrs)][1:min(30, .N)])
cat("\n")

cat("=== QC: how many morph methylation values are present per consensus region ===\n")
agg[, n_meth_nonNA := rowSums(!is.na(.SD)), .SDcols = c("meth_female","meth_immaculata","meth_parae","meth_yellow")]
print(agg[, .N, by = n_meth_nonNA][order(n_meth_nonNA)])
cat("\n")

cat("=== QC: chromosome_type breakdown ===\n")
print(agg[, .N, by = chromosome_type][order(-N)])
cat("\n")


# QC for cons_id
stopifnot(agg[, anyDuplicated(cons_id)] == 0)

# -----------------------------
# Write output
# -----------------------------
fwrite(agg, out_tsv, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote master consensus matrix:\n  ", normalizePath(out_tsv, winslash="/", mustWork=FALSE), "\n", sep = "")
