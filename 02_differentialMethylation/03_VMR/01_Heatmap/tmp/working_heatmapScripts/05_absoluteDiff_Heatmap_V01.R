#!/usr/bin/env Rscript

# DATE: 20260721
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Strict unsupervised clustering of VMR methylation profiles.
#
#   Feature selection and clustering use ONLY the 12 individual-sample
#   methylation values. Morph-level means are NOT used.
#
#   Steps:
#     1. Identify the 12 individual-sample methylation columns.
#     2. For each VMR, calculate:
#
#          sample_range = max(sample methylation) - min(sample methylation)
#
#     3. Retain VMRs where sample_range >= DIFF_THRESH.
#     4. Cluster VMRs using methylation across the 12 individual samples.
#     5. Cluster samples using methylation across retained VMRs.
#     6. Add morph annotation ONLY AFTER clustering for visualization.
#
#   Therefore:
#     - Feature selection does not use morph identity.
#     - Row clustering does not use morph identity.
#     - Column/sample clustering does not use morph identity.
#
#   The morph annotation is visual only and has no effect on clustering.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(viridis)
})

# ==============================================================================
# User settings
# ==============================================================================

DIFF_THRESH <- 0.40

# Require methylation information for all samples before a VMR can be used.
# This avoids clustering being influenced by different amounts of missing data.
REQUIRE_COMPLETE <- TRUE

# Hierarchical clustering settings
DIST_METHOD   <- "euclidean"
HCLUST_METHOD <- "complete"

# ==============================================================================
# Input / output
# ==============================================================================

infile <- "vmr_methylation/vmr.morph_and_sample_observed.tsv"

outdir <- paste0(
  "../plots/strict_unsupervised_sampleRange_",
  DIFF_THRESH
)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ==============================================================================
# Read data
# ==============================================================================

dt <- fread(infile)

# Recalculate chromosome type directly from chromosome name
dt[, chromosome_type :=
     ifelse(chr == "Parae_12", "Sex_Ch", "Autosome")]

# Unique row identifier
if ("cons_id" %in% names(dt)) {
  dt[, row_id := cons_id]
} else {
  dt[, row_id := paste0(chr, ":", start, "-", end)]
}

# ==============================================================================
# Identify individual-sample methylation columns
# ==============================================================================

sample_cols <- sort(
  grep(
    "^meth_ppar.*mem[0-9]+$",
    names(dt),
    value = TRUE
  )
)

if (length(sample_cols) == 0) {
  stop("No individual-sample methylation columns were detected.")
}

message(
  "Detected ", length(sample_cols),
  " sample methylation columns:"
)

message(
  paste(sample_cols, collapse = "\n")
)

# Labels displayed in the heatmap
sample_labels <- sub("^meth_", "", sample_cols)

# Check required columns
required_cols <- c(
  "chr",
  "start",
  "end",
  "row_id",
  "chromosome_type",
  sample_cols
)

missing_cols <- setdiff(required_cols, names(dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ==============================================================================
# STRICT UNSUPERVISED FEATURE FILTER
# ==============================================================================

# IMPORTANT:
#
# No morph means are used here.
#
# For each VMR:
#
#   sample_range =
#       maximum methylation among individual samples
#       -
#       minimum methylation among individual samples
#
# Thus, selection depends only on variability among individuals.

dt[, n_meth_nonNA_samples :=
     rowSums(!is.na(.SD)),
   .SDcols = sample_cols
]

# Row-wise maximum
dt[, sample_max :=
     do.call(
       pmax,
       c(.SD, na.rm = TRUE)
     ),
   .SDcols = sample_cols
]

# Row-wise minimum
dt[, sample_min :=
     do.call(
       pmin,
       c(.SD, na.rm = TRUE)
     ),
   .SDcols = sample_cols
]

# Replace infinities generated when all observations are NA
dt[!is.finite(sample_max), sample_max := NA_real_]
dt[!is.finite(sample_min), sample_min := NA_real_]

# Calculate full range across individuals
dt[, sample_range := sample_max - sample_min]

# ==============================================================================
# Apply filter
# ==============================================================================

if (REQUIRE_COMPLETE) {
  
  # Require information from every individual sample
  dt_plot <- dt[
    n_meth_nonNA_samples == length(sample_cols) &
      !is.na(sample_range) &
      sample_range >= DIFF_THRESH
  ]
  
} else {
  
  # Alternative:
  # require at least two individuals with methylation information
  dt_plot <- dt[
    n_meth_nonNA_samples >= 2 &
      !is.na(sample_range) &
      sample_range >= DIFF_THRESH
  ]
}

message(
  "\nVMRs before filtering: ",
  nrow(dt)
)

message(
  "VMRs after sample-range filtering: ",
  nrow(dt_plot)
)

message(
  "Fraction retained: ",
  round(
    100 * nrow(dt_plot) / nrow(dt),
    2
  ),
  "%"
)

# ==============================================================================
# QC outputs
# ==============================================================================

qc_counts <- rbindlist(
  list(
    
    dt[
      ,
      .(
        stage = "before_filter",
        n = .N
      ),
      by = chromosome_type
    ],
    
    dt_plot[
      ,
      .(
        stage = "after_filter",
        n = .N
      ),
      by = chromosome_type
    ]
    
  )
)

fwrite(
  qc_counts,
  file.path(
    outdir,
    "QC_counts_before_after_unsupervised_filter.tsv"
  ),
  sep = "\t"
)

# Distribution summary
qc_range_before <- dt[
  !is.na(sample_range),
  .(
    stage  = "before_filter",
    min    = min(sample_range),
    q25    = quantile(sample_range, 0.25),
    median = median(sample_range),
    mean   = mean(sample_range),
    q75    = quantile(sample_range, 0.75),
    max    = max(sample_range)
  )
]

qc_range_after <- dt_plot[
  ,
  .(
    stage  = "after_filter",
    min    = min(sample_range),
    q25    = quantile(sample_range, 0.25),
    median = median(sample_range),
    mean   = mean(sample_range),
    q75    = quantile(sample_range, 0.75),
    max    = max(sample_range)
  )
]

fwrite(
  rbind(
    qc_range_before,
    qc_range_after
  ),
  file.path(
    outdir,
    "QC_sampleRange_summary.tsv"
  ),
  sep = "\t"
)

# Save filtered VMR table
fwrite(
  dt_plot,
  file.path(
    outdir,
    paste0(
      "vmr.strictUnsupervised.sampleRangeGE_",
      DIFF_THRESH,
      ".tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ==============================================================================
# Fixed methylation color scale: 0 - 1
# ==============================================================================

n_colors <- 101

color_palette <- inferno(n_colors)

breaks <- seq(
  0,
  1,
  length.out = n_colors + 1
)

# ==============================================================================
# Sample annotation
# ==============================================================================

# Morph identity is inferred ONLY for annotation.
#
# It is NOT supplied to:
#   - filtering
#   - distance calculations
#   - hierarchical clustering

get_morph <- function(x) {
  
  out <- rep(NA_character_, length(x))
  
  out[grepl("^pparfmem", x)] <- "Female"
  out[grepl("^pparimem", x)] <- "Immaculata"
  out[grepl("^pparpmem", x)] <- "Parae"
  out[grepl("^pparymem", x)] <- "Yellow"
  
  out
}

sample_annotation <- data.frame(
  Morph = get_morph(sample_labels)
)

rownames(sample_annotation) <- sample_labels

# Annotation colors
morph_levels <- c(
  "Female",
  "Immaculata",
  "Parae",
  "Yellow"
)

sample_annotation$Morph <- factor(
  sample_annotation$Morph,
  levels = morph_levels
)

annotation_colors <- list(
  Morph = setNames(
    viridis(length(morph_levels)),
    morph_levels
  )
)

# ==============================================================================
# Heatmap function
# ==============================================================================

plot_unsupervised_heatmap <- function(
    chrom_type,
    title,
    pdf_name,
    row_order_tsv,
    sample_order_tsv) {
  
  # --------------------------------------------------------------------------
  # Subset chromosome category
  # --------------------------------------------------------------------------
  
  sub <- dt_plot[
    chromosome_type == chrom_type
  ]
  
  setorder(
    sub,
    chr,
    start,
    end
  )
  
  if (nrow(sub) < 2) {
    
    warning(
      "Not enough VMRs for ",
      chrom_type,
      " after filtering. Skipping."
    )
    
    return(invisible(NULL))
  }
  
  # --------------------------------------------------------------------------
  # Matrix containing ONLY individual samples
  # --------------------------------------------------------------------------
  
  mat <- as.matrix(
    sub[, ..sample_cols]
  )
  
  storage.mode(mat) <- "double"
  
  rownames(mat) <- sub$row_id
  colnames(mat) <- sample_labels
  
  mat[!is.finite(mat)] <- NA_real_
  
  # --------------------------------------------------------------------------
  # Strict unsupervised ROW clustering
  #
  # Distance among VMRs based only on their methylation across the 12 samples.
  # --------------------------------------------------------------------------
  
  row_dist <- dist(
    mat,
    method = DIST_METHOD
  )
  
  hc_rows <- hclust(
    row_dist,
    method = HCLUST_METHOD
  )
  
  # --------------------------------------------------------------------------
  # Strict unsupervised SAMPLE clustering
  #
  # Transpose matrix:
  #   rows    = samples
  #   columns = VMRs
  #
  # Morph identity is NOT used here.
  # --------------------------------------------------------------------------
  
  sample_dist <- dist(
    t(mat),
    method = DIST_METHOD
  )
  
  hc_cols <- hclust(
    sample_dist,
    method = HCLUST_METHOD
  )
  
  # --------------------------------------------------------------------------
  # Save unsupervised VMR ordering
  # --------------------------------------------------------------------------
  
  ordered_row_ids <- rownames(mat)[hc_rows$order]
  
  row_order <- data.table(
    Heatmap_order = seq_along(ordered_row_ids),
    row_id = ordered_row_ids
  )
  
  row_order <- merge(
    row_order,
    sub[
      ,
      .(
        row_id,
        chr,
        start,
        end,
        chromosome_type,
        sample_range
      )
    ],
    by = "row_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  setorder(
    row_order,
    Heatmap_order
  )
  
  fwrite(
    row_order,
    file.path(
      outdir,
      row_order_tsv
    ),
    sep = "\t"
  )
  
  # --------------------------------------------------------------------------
  # Save unsupervised sample ordering
  # --------------------------------------------------------------------------
  
  ordered_samples <- colnames(mat)[hc_cols$order]
  
  sample_order <- data.table(
    Dendrogram_order = seq_along(ordered_samples),
    Sample = ordered_samples
  )
  
  # Morph added AFTER clustering only for interpretation
  sample_order[, Morph := get_morph(Sample)]
  
  fwrite(
    sample_order,
    file.path(
      outdir,
      sample_order_tsv
    ),
    sep = "\t"
  )
  
  # Print sample order to log
  message(
    "\n",
    chrom_type,
    " sample dendrogram order:\n",
    paste(
      paste0(
        sample_order$Sample,
        " [",
        sample_order$Morph,
        "]"
      ),
      collapse = "\n"
    )
  )
  
  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------
  
  pdf(
    file.path(
      outdir,
      pdf_name
    ),
    width = 7,
    height = 8
  )
  
  pheatmap(
    mat,
    
    # Use EXACT dendrograms calculated above.
    # Do not use cluster_rows = TRUE / cluster_cols = TRUE here,
    # because that would make pheatmap recalculate them.
    cluster_rows = hc_rows,
    cluster_cols = hc_cols,
    
    # Morph is annotation ONLY
    annotation_col = sample_annotation,
    annotation_colors = annotation_colors,
    
    treeheight_row = 30,
    treeheight_col = 45,
    
    main = title,
    
    fontsize = 10,
    fontsize_row = 6,
    fontsize_col = 8,
    
    angle_col = 90,
    
    color = color_palette,
    breaks = breaks,
    
    na_col = "black",
    
    show_rownames = FALSE,
    border_color = NA
  )
  
  dev.off()
}

# ==============================================================================
# Run heatmaps
# ==============================================================================

# ------------------------------------------------------------------------------
# Autosomes
# ------------------------------------------------------------------------------

plot_unsupervised_heatmap(
  
  chrom_type = "Autosome",
  
  title = paste0(
    "Strict unsupervised VMR methylation — Autosomes\n",
    "VMRs selected by individual sample range ≥ ",
    DIFF_THRESH
  ),
  
  pdf_name = paste0(
    "heatmap_AUTOSOME.strictUnsupervised.sampleRangeGE_",
    DIFF_THRESH,
    ".pdf"
  ),
  
  row_order_tsv = paste0(
    "VMR_order_AUTOSOME.strictUnsupervised.sampleRangeGE_",
    DIFF_THRESH,
    ".tsv"
  ),
  
  sample_order_tsv = paste0(
    "sample_order_AUTOSOME.strictUnsupervised.sampleRangeGE_",
    DIFF_THRESH,
    ".tsv"
  )
)

# ------------------------------------------------------------------------------
# Chromosome 12
# ------------------------------------------------------------------------------

plot_unsupervised_heatmap(
  
  chrom_type = "Sex_Ch",
  
  title = paste0(
    "Strict unsupervised VMR methylation — Chromosome 12\n",
    "VMRs selected by individual sample range ≥ ",
    DIFF_THRESH
  ),
  
  pdf_name = paste0(
    "heatmap_CH12.strictUnsupervised.sampleRangeGE_",
    DIFF_THRESH,
    ".pdf"
  ),
  
  row_order_tsv = paste0(
    "VMR_order_CH12.strictUnsupervised.sampleRangeGE_",
    DIFF_THRESH,
    ".tsv"
  ),
  
  sample_order_tsv = paste0(
    "sample_order_CH12.strictUnsupervised.sampleRangeGE_",
    DIFF_THRESH,
    ".tsv"
  )
)