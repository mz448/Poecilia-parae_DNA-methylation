#!/usr/bin/env Rscript

# DATE: 20260721
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Strict unsupervised clustering of VMR methylation profiles using
#   standard deviation (SD) across the 12 individual samples for filtering.
#
#   Steps:
#     1. Calculate SD of methylation across the 12 samples for each VMR.
#     2. Retain VMRs with SD >= SD_THRESH.
#     3. Cluster VMRs using the 12 individual methylation values.
#     4. Cluster samples using the retained VMRs.
#     5. Add morph identity ONLY as a visual annotation after clustering.
#
#   Morph means and morph identities are NOT used for:
#     - VMR filtering
#     - row clustering
#     - sample clustering
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(viridis)
})


# ==============================================================================
# User settings
# ==============================================================================

# SD_THRESH <- 0.15
SD_THRESH_MAX <- 0.5
SD_THRESH_MIN <- 0.7

DIST_METHOD   <- "euclidean"
HCLUST_METHOD <- "complete"


# ==============================================================================
# Input / output
# ==============================================================================

infile <- "vmr_methylation/vmr.morph_and_sample_observed.tsv"


outdir <- paste0("../plots/SD/", SD_THRESH_MIN,"-",SD_THRESH_MAX)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ==============================================================================
# Read data
# ==============================================================================

dt <- fread(infile)

# Define chromosome type directly from chromosome
dt[, chromosome_type :=
     ifelse(chr == "Parae_12", "Sex_Ch", "Autosome")]

# Unique VMR ID
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
  stop("No individual-sample methylation columns detected.")
}

message(
  "Detected ",
  length(sample_cols),
  " individual sample columns:"
)

message(
  paste(sample_cols, collapse = "\n")
)

# Expected: 12 samples
if (length(sample_cols) != 12) {
  warning(
    "Expected 12 sample columns, but detected ",
    length(sample_cols),
    ". Check column names."
  )
}

sample_labels <- sub(
  "^meth_",
  "",
  sample_cols
)


# ==============================================================================
# Calculate SD across individual samples
# ==============================================================================

# Matrix containing ONLY individual samples
sample_mat <- as.matrix(
  dt[, ..sample_cols]
)

storage.mode(sample_mat) <- "double"

# Number of available samples per VMR
dt[, n_meth_nonNA_samples :=
     rowSums(!is.na(sample_mat))]

# Calculate standard deviation across the 12 samples
dt[, sample_SD :=
     apply(
       sample_mat,
       1,
       sd,
       na.rm = TRUE
     )]

# SD cannot be calculated from fewer than 2 observations
dt[n_meth_nonNA_samples < 2, sample_SD := NA_real_]


# ==============================================================================
# Strict unsupervised SD filter
# ==============================================================================

# Require all 12 samples to have methylation information
dt_plot <- dt[
  n_meth_nonNA_samples == length(sample_cols) &
    !is.na(sample_SD) &
    sample_SD >= SD_THRESH_MAX &
    sample_SD <= SD_THRESH_MIN
]


message(
  "\nVMRs before filtering: ",
  nrow(dt)
)

message(
  "VMRs after SD filtering: ",
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
    "QC_counts_before_after_SD_filter.tsv"
  ),
  sep = "\t"
)


# SD distribution summary
qc_sd_before <- dt[
  !is.na(sample_SD),
  .(
    stage  = "before_filter",
    min    = min(sample_SD),
    q25    = quantile(sample_SD, 0.25),
    median = median(sample_SD),
    mean   = mean(sample_SD),
    q75    = quantile(sample_SD, 0.75),
    max    = max(sample_SD)
  )
]

qc_sd_after <- dt_plot[
  ,
  .(
    stage  = "after_filter",
    min    = min(sample_SD),
    q25    = quantile(sample_SD, 0.25),
    median = median(sample_SD),
    mean   = mean(sample_SD),
    q75    = quantile(sample_SD, 0.75),
    max    = max(sample_SD)
  )
]

fwrite(
  rbind(
    qc_sd_before,
    qc_sd_after
  ),
  file.path(
    outdir,
    "QC_sampleSD_summary.tsv"
  ),
  sep = "\t"
)


# Save filtered VMRs
fwrite(
  dt_plot,
  file.path(
    outdir,
    paste0(
      "vmr.strictUnsupervised.sampleSD_GE_",
      SD_THRESH,
      ".tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ==============================================================================
# Fixed methylation color scale: 0–1
# ==============================================================================

n_colors <- 101

color_palette <- inferno(
  n_colors
)

breaks <- seq(
  0,
  1,
  length.out = n_colors + 1
)


# ==============================================================================
# Morph annotation
# ==============================================================================

# IMPORTANT:
# Morph identity is defined only for visualization.
# It does NOT enter the filtering or clustering calculations.

get_morph <- function(x) {
  
  out <- rep(
    NA_character_,
    length(x)
  )
  
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
      " after filtering."
    )
    
    return(invisible(NULL))
  }
  
  
  # --------------------------------------------------------------------------
  # Sample methylation matrix
  # --------------------------------------------------------------------------
  
  mat <- as.matrix(
    sub[, ..sample_cols]
  )
  
  storage.mode(mat) <- "double"
  
  rownames(mat) <- sub$row_id
  colnames(mat) <- sample_labels
  
  
  # --------------------------------------------------------------------------
  # Cluster VMRs
  # --------------------------------------------------------------------------
  
  hc_rows <- hclust(
    dist(
      mat,
      method = DIST_METHOD
    ),
    method = HCLUST_METHOD
  )
  
  
  # --------------------------------------------------------------------------
  # Cluster samples
  #
  # t(mat):
  #   rows    = 12 samples
  #   columns = retained VMRs
  # --------------------------------------------------------------------------
  
  hc_cols <- hclust(
    dist(
      t(mat),
      method = DIST_METHOD
    ),
    method = HCLUST_METHOD
  )
  
  
  # --------------------------------------------------------------------------
  # Save sample dendrogram order
  # --------------------------------------------------------------------------
  
  ordered_samples <- colnames(mat)[
    hc_cols$order
  ]
  
  sample_order <- data.table(
    Dendrogram_order = seq_along(ordered_samples),
    Sample = ordered_samples
  )
  
  # Add morph AFTER clustering
  sample_order[, Morph :=
                 get_morph(Sample)]
  
  fwrite(
    sample_order,
    file.path(
      outdir,
      sample_order_tsv
    ),
    sep = "\t"
  )
  
  
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
  # Plot heatmap
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
    
    # Use explicitly calculated unsupervised dendrograms
    cluster_rows = hc_rows,
    cluster_cols = hc_cols,
    
    # Morph annotation is visual only
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
# Autosomal VMRs
# ==============================================================================

plot_unsupervised_heatmap(
  
  chrom_type = "Autosome",
  
  title = paste0(
    "Strict unsupervised VMR methylation — Autosomes\n",
    "VMRs selected by sample SD ≥ ",
    SD_THRESH
  ),
  
  pdf_name = paste0(
    "heatmap_AUTOSOME.strictUnsupervised.SD_GE_",
    SD_THRESH,
    ".pdf"
  ),
  
  sample_order_tsv = paste0(
    "sample_order_AUTOSOME.strictUnsupervised.SD_GE_",
    SD_THRESH,
    ".tsv"
  )
)


# ==============================================================================
# Chromosome 12 VMRs
# ==============================================================================

plot_unsupervised_heatmap(
  
  chrom_type = "Sex_Ch",
  
  title = paste0(
    "Strict unsupervised VMR methylation — Chromosome 12\n",
    "VMRs selected by sample SD ≥ ",
    SD_THRESH
  ),
  
  pdf_name = paste0(
    "heatmap_CH12.strictUnsupervised.SD_GE_",
    SD_THRESH,
    ".pdf"
  ),
  
  sample_order_tsv = paste0(
    "sample_order_CH12.strictUnsupervised.SD_GE_",
    SD_THRESH,
    ".tsv"
  )
)