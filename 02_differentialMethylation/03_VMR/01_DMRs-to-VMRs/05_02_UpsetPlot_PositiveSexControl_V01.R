#!/usr/bin/env Rscript
# DATE:   20261026
# AUTHOR: MZF
# SCRIPT: 05_02_UpsetPlot_PositiveSexControl_V01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Build a pooled shuffled VMR x contrast membership table for UpSet plotting.
#
#   INPUTS:
#     1) vmr_shuffled_sex.bed
#        cols: chr, start, end, cons_id, dmr_ids
#
#     2) dmr_shuffled_sex.withIDs.tsv
#        contains: dmr_id, comparison, condition, ...
#
#   SELECTED CONDITIONS:
#     - shuffled_sex
#
#   CONTRASTS:
#     shuffled_sex:
#                   - f_vs_sM_01
#                   - f_vs_sM_02
#                   - f_vs_sM_03
#                   - f_vs_sM_04
#                   - f_vs_sM_05
#                   - f_vs_sM_06
#
#   OUTPUTS:
#     A) long table:
#        one row per cons_id x contributing DMR x contributing contrast,
#        retaining the shuffled condition of origin
#
#     B) binary wide table:
#        one row per cons_id, one 0/1 column per shuffled contrast
#
#     C) contrast counts:
#        number of VMRs containing each shuffled contrast
#
#   LOGIC:
#     - split collapsed dmr_ids from the pooled VMR BED
#     - join each dmr_id back to dmr_shuffled_sex.withIDs.tsv
#     - recover the condition and contrast contributing to each VMR
#     - retain all six shuffled comparisons as independent sets
#     - collapse to unique cons_id x contrast combinations
#     - cast to a wide binary matrix for UpSet plotting
#
#   NOTE:
#     Comparison suffixes are intentionally retained:
#       f_vs_sM_01
#       f_vs_sM_02
#       f_vs_sM_03
#       f_vs_sM_04
#       f_vs_sM_05
#       f_vs_sM_06
#
#     Unlike the corresponding real-data script, trailing replicate suffixes
#     MUST NOT be removed because _01, _02, and _03 define the six independent
#     shuffled comparison sets used in the UpSet analysis.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

# -----------------------------
# Defaults
# -----------------------------
default_outdir   <- "."
default_cons_bed <- "vmr_shuffled_sex.bed"
default_dmr_tsv  <- "dmr_shuffled_sex.withIDs.tsv"
default_prefix   <- "vmr_shuffled_sex.contrastSharing"

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
# Expected conditions/contrasts
# -----------------------------
expected_conditions <- c(
  "shuffled_sex" 
)

expected_contrasts <- c(
  "f_vs_sM_01",
  "f_vs_sM_02",
  "f_vs_sM_03",
  "f_vs_sM_04",
  "f_vs_sM_05",
  "f_vs_sM_06"
)

# -----------------------------
# Read inputs
# -----------------------------
cons <- fread(cons_bed, header = FALSE)

if (ncol(cons) < 5) {
  stop(
    "Pooled VMR BED must have at least 5 columns: ",
    "chr, start, end, cons_id, dmr_ids"
  )
}

setnames(
  cons,
  1:5,
  c("chr", "start", "end", "cons_id", "dmr_ids")
)

cons[, start := as.integer(start)]
cons[, end   := as.integer(end)]

dmr <- fread(dmr_tsv)

required_cols <- c(
  "dmr_id",
  "comparison",
  "condition"
)

missing_cols <- setdiff(required_cols, names(dmr))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in dmr_shuffled_sex.withIDs.tsv: ",
    paste(missing_cols, collapse = ", ")
  )
}

# -----------------------------
# Validate conditions
# -----------------------------
observed_conditions <- sort(unique(dmr$condition))

unexpected_conditions <- setdiff(
  observed_conditions,
  expected_conditions
)

missing_conditions <- setdiff(
  expected_conditions,
  observed_conditions
)

if (length(unexpected_conditions) > 0) {
  stop(
    "Unexpected conditions found in DMR table: ",
    paste(unexpected_conditions, collapse = ", "),
    "\nExpected ONLY: ",
    paste(expected_conditions, collapse = ", ")
  )
}

if (length(missing_conditions) > 0) {
  stop(
    "Expected conditions missing from DMR table: ",
    paste(missing_conditions, collapse = ", ")
  )
}

# -----------------------------
# Validate comparisons
# -----------------------------
# IMPORTANT:
# Comparison names are used exactly as present in the input.
# The _01/_02/_03 suffixes are NOT removed.
# -----------------------------
observed_contrasts <- sort(unique(dmr$comparison))

unexpected_contrasts <- setdiff(
  observed_contrasts,
  expected_contrasts
)

missing_contrasts <- setdiff(
  expected_contrasts,
  observed_contrasts
)

if (length(unexpected_contrasts) > 0) {
  stop(
    "Unexpected comparisons found in DMR table: ",
    paste(unexpected_contrasts, collapse = ", "),
    "\nExpected ONLY the six shuffled comparisons."
  )
}

if (length(missing_contrasts) > 0) {
  stop(
    "Expected shuffled comparisons missing from DMR table: ",
    paste(missing_contrasts, collapse = ", ")
  )
}

# -----------------------------
# Validate condition x contrast
# correspondence
# -----------------------------
# shuffled_all should contain ONLY sA_vs_sA_*
# shuffled_morph should contain ONLY sM_vs_sM_*
# -----------------------------
bad_condition_contrast <- dmr[
  (condition == "shuffled_sex" &
     !comparison %in% c(
       "f_vs_sM_01",
       "f_vs_sM_02",
       "f_vs_sM_03",
       "f_vs_sM_04",
       "f_vs_sM_05",
       "f_vs_sM_06"
     )) 
]

if (nrow(bad_condition_contrast) > 0) {
  cat(
    "ERROR: unexpected condition x comparison combinations detected.\n"
  )
  print(
    unique(
      bad_condition_contrast[, .(condition, comparison)]
    )
  )
  stop(
    "Check condition/comparison assignments in the shuffled DMR table."
  )
}

# -----------------------------
# Build DMR lookup table
# -----------------------------
# Keep one row per dmr_id.
#
# Unlike the real-data script, comparison names are NOT cleaned because
# _01/_02/_03 distinguish the six independent shuffled comparisons.
# -----------------------------
dmr_map <- unique(
  dmr[, .(
    dmr_id,
    condition,
    comparison
  )]
)

# Each stable dmr_id should map to exactly one condition/comparison.
dmr_id_duplicates <- dmr_map[
  ,
  .N,
  by = dmr_id
][N > 1]

if (nrow(dmr_id_duplicates) > 0) {
  cat(
    "ERROR: some dmr_id values map to more than one ",
    "condition/comparison combination.\n"
  )
  print(head(dmr_id_duplicates, 20))
  stop("DMR IDs must uniquely identify input DMR rows.")
}

# -----------------------------
# Expand pooled VMR mapping:
# one row per cons_id x dmr_id
# -----------------------------
map <- cons[, .(
  dmr_id = unlist(
    strsplit(dmr_ids, ",", fixed = TRUE)
  )
), by = .(
  chr,
  start,
  end,
  cons_id,
  dmr_ids
)]

expected_pairs <- cons[
  ,
  sum(
    lengths(
      strsplit(dmr_ids, ",", fixed = TRUE)
    )
  )
]

if (nrow(map) != expected_pairs) {
  stop(
    "Mapping expansion mismatch: expected ",
    expected_pairs,
    " cons_id-to-dmr_id pairs but got ",
    nrow(map)
  )
}

# -----------------------------
# Join dmr_id back to condition
# and contrast information
# -----------------------------
joined <- merge(
  map,
  dmr_map,
  by = "dmr_id",
  all.x = TRUE,
  sort = FALSE
)

missing_ids <- joined[
  is.na(comparison),
  unique(dmr_id)
]

if (length(missing_ids) > 0) {
  cat(
    "ERROR: some dmr_id values from pooled VMR BED were not found in ",
    "dmr_shuffled_sex.withIDs.tsv\n"
  )
  cat("First missing IDs:\n")
  print(head(missing_ids, 20))
  
  stop(
    "Fix the consistency between vmr_shuffled_sex.bed and ",
    "dmr_shuffled_sex.withIDs.tsv"
  )
}

# derive chromosome type from VMR coordinates
joined[
  ,
  chromosome_type := ifelse(
    chr == "Parae_12",
    "Sex_Ch",
    "Autosome"
  )
]

# -----------------------------
# Long table:
# one row per contributing
# cons_id x dmr_id x contrast
#
# condition is retained to preserve
# shuffled dataset provenance.
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
    condition,
    comparison
  )]
)

# -----------------------------
# Binary table for UpSet
#
# Each of the six shuffled comparisons
# is treated as an independent set.
# -----------------------------
binary_long <- unique(
  contrast_long[, .(
    chr,
    start,
    end,
    cons_id,
    dmr_ids,
    chromosome_type,
    contrast = comparison,
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
meta_cols <- c(
  "chr",
  "start",
  "end",
  "cons_id",
  "dmr_ids",
  "chromosome_type"
)

contrast_cols <- setdiff(
  names(wide),
  meta_cols
)

preferred_order <- c(
  "f_vs_sM_01",
  "f_vs_sM_02",
  "f_vs_sM_03",
  "f_vs_sM_04",
  "f_vs_sM_05",
  "f_vs_sM_06"
)

contrast_cols <- c(
  intersect(preferred_order, contrast_cols),
  setdiff(sort(contrast_cols), preferred_order)
)

setcolorder(
  wide,
  c(meta_cols, contrast_cols)
)

# -----------------------------
# Add useful summary columns
# -----------------------------
wide[
  ,
  n_contrasts := rowSums(.SD),
  .SDcols = contrast_cols
]

wide[
  ,
  contrast_list := apply(
    .SD,
    1,
    function(x) {
      hits <- contrast_cols[as.logical(x)]
      
      if (length(hits) == 0) {
        return("none")
      }
      
      paste(hits, collapse = ",")
    }
  ),
  .SDcols = contrast_cols
]

wide[
  ,
  contrast_combo := apply(
    .SD,
    1,
    function(x) {
      hits <- contrast_cols[as.logical(x)]
      
      if (length(hits) == 0) {
        return("none")
      }
      
      paste(hits, collapse = "&")
    }
  ),
  .SDcols = contrast_cols
]

# move summary columns next to metadata
setcolorder(
  wide,
  c(
    "chr",
    "start",
    "end",
    "cons_id",
    "dmr_ids",
    "chromosome_type",
    "n_contrasts",
    "contrast_list",
    "contrast_combo",
    contrast_cols
  )
)

# -----------------------------
# Simple counts per contrast
# -----------------------------
# Keep n_consensus_dmrs naming for
# backward compatibility.
# -----------------------------
contrast_counts <- data.table(
  contrast = contrast_cols,
  n_consensus_dmrs = vapply(
    contrast_cols,
    function(cc) {
      sum(wide[[cc]] == 1L)
    },
    numeric(1)
  )
)

# -----------------------------
# QC
# -----------------------------
cat("=== SETTINGS ===\n")

cat(
  "Pooled VMR BED: ",
  normalizePath(
    cons_bed,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "DMR TSV:        ",
  normalizePath(
    dmr_tsv,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Output prefix:  ",
  normalizePath(
    prefix,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n\n",
  sep = ""
)

cat("=== QC: conditions recovered ===\n")
print(
  dmr[
    ,
    .N,
    by = condition
  ][order(condition)]
)
cat("\n")

cat("=== QC: condition x contrast DMR counts ===\n")
print(
  dmr[
    ,
    .N,
    by = .(condition, comparison)
  ][order(condition, comparison)]
)
cat("\n")

cat("=== QC: counts ===\n")
cat(
  "Pooled VMR regions:                ",
  nrow(cons),
  "\n",
  sep = ""
)
cat(
  "Expanded cons_id -> dmr_id pairs:  ",
  nrow(map),
  "\n",
  sep = ""
)
cat(
  "Unique cons_id x contrast rows:    ",
  nrow(binary_long),
  "\n",
  sep = ""
)
cat(
  "Binary matrix rows:                ",
  nrow(wide),
  "\n\n",
  sep = ""
)

cat("=== QC: contrasts recovered ===\n")
print(
  contrast_counts[
    order(-n_consensus_dmrs)
  ]
)
cat("\n")

cat("=== QC: number of contrasts per VMR ===\n")
print(
  wide[
    ,
    .N,
    by = n_contrasts
  ][order(n_contrasts)]
)
cat("\n")

cat("=== QC: chromosome breakdown ===\n")
print(
  wide[
    ,
    .N,
    by = chromosome_type
  ][order(-N)]
)
cat("\n")

cat("=== QC: VMR sharing between shuffled datasets ===\n")

# Determine whether each VMR contains at least one shuffled_all
# contrast (sA) and/or at least one shuffled_morph contrast (sM).
wide[
  ,
  has_shuffled_sex :=
    f_vs_sM_01 == 1L |
    f_vs_sM_02 == 1L |
    f_vs_sM_03 == 1L |
    f_vs_sM_04 == 1L |
    f_vs_sM_05 == 1L |
    f_vs_sM_06 == 1L
]


sharing_qc <- wide[
  ,
  .N,
  by = .(
    has_shuffled_sex
  )
][
  order(
    -has_shuffled_sex
  )
]

print(sharing_qc)
cat("\n")

# Remove temporary QC columns from final binary output.
wide[
  ,
  c(
    "has_shuffled_sex"
  ) := NULL
]

# -----------------------------
# Final integrity checks
# -----------------------------
stopifnot(
  anyDuplicated(wide$cons_id) == 0
)

stopifnot(
  nrow(wide) == nrow(cons)
)

stopifnot(
  setequal(
    contrast_cols,
    expected_contrasts
  )
)

# Every VMR should contain at least one
# of the six selected shuffled contrasts.
stopifnot(
  all(wide$n_contrasts >= 1L)
)

# -----------------------------
# Write outputs
# -----------------------------
fwrite(
  contrast_long,
  out_long,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  wide,
  out_binary,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  contrast_counts,
  out_counts,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote:\n")

cat(
  "  ",
  normalizePath(
    out_long,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "  ",
  normalizePath(
    out_binary,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)

cat(
  "  ",
  normalizePath(
    out_counts,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)