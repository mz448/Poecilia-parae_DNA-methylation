#!/usr/bin/env Rscript

# DATE: 20260721
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Strict label-independent filtering of VMRs using a blind 3-vs-3 comparison.
#
#   For each VMR:
#     1. Sort methylation values across the 12 individual samples.
#     2. Calculate the mean of the lowest 3 values.
#     3. Calculate the mean of the highest 3 values.
#     4. Retain the VMR when:
#
#          mean(highest 3) - mean(lowest 3) >= DIFF_THRESH
#
#   No morph means, morph labels, or within-group SD values are used.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(viridis)
})

# ==============================================================================
# User settings
# ==============================================================================

DIFF_THRESH <- 0.58
# 
# DIST_METHOD   <- "euclidean"
# HCLUST_METHOD <- "complete"
# 
# # ==============================================================================
# # Input / output
# # ==============================================================================
# 
infile <- "vmr_methylation/vmr.morph_and_sample_observed.tsv"

outdir <- paste0("../plots/3vs3_extremes/", DIFF_THRESH)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)
# 
# # ==============================================================================
# # Read data
# # ==============================================================================
# 
# dt <- fread(infile)
# 
# # Define chromosome category directly from chromosome name
# dt[, chromosome_type :=
#      ifelse(chr == "Parae_12", "Sex_Ch", "Autosome")]
# 
# # Unique VMR identifier
# if ("cons_id" %in% names(dt)) {
#   dt[, row_id := cons_id]
# } else {
#   dt[, row_id := paste0(chr, ":", start, "-", end)]
# }
# 
# # ==============================================================================
# # Identify individual-sample methylation columns
# # ==============================================================================
# 
# sample_cols <- sort(
#   grep(
#     "^meth_ppar.*mem[0-9]+$",
#     names(dt),
#     value = TRUE
#   )
# )
# 
# if (length(sample_cols) != 12) {
#   stop(
#     "Expected 12 individual-sample methylation columns, but detected ",
#     length(sample_cols),
#     ":\n",
#     paste(sample_cols, collapse = "\n")
#   )
# }
# 
# sample_labels <- sub("^meth_", "", sample_cols)
# 
# message("Detected sample columns:")
# message(paste(sample_cols, collapse = "\n"))
# 
# # ==============================================================================
# # Construct individual-sample methylation matrix
# # ==============================================================================
# 
# mat_all <- as.matrix(
#   dt[, ..sample_cols]
# )
# 
# storage.mode(mat_all) <- "double"
# 
# dt[, n_meth_nonNA_samples :=
#      rowSums(!is.na(mat_all))]
# 
# # Require all 12 samples for the blind 3-vs-3 comparison
# complete_rows <- dt$n_meth_nonNA_samples == length(sample_cols)
# 
# # ==============================================================================
# # Blind 3-vs-3 for extreme filter
# # ==============================================================================
# 
# dt[, `:=`(
#   mean_low3       = NA_real_,
#   mean_high3      = NA_real_,
#   extreme_3vs3_diff = NA_real_,
#   low3_samples    = NA_character_,
#   high3_samples   = NA_character_
# )]
# 
# for (i in which(complete_rows)) {
# 
#   vals <- mat_all[i, ]
# 
#   # Sort samples from lowest to highest methylation
#   ord <- order(vals)
# 
#   low3_idx  <- ord[1:3]
#   high3_idx <- ord[10:12]
# 
#   low3_values  <- vals[low3_idx]
#   high3_values <- vals[high3_idx]
# 
#   dt$mean_low3[i]  <- mean(low3_values)
#   dt$mean_high3[i] <- mean(high3_values)
# 
#   dt$extreme_3vs3_diff[i] <-
#     dt$mean_high3[i] - dt$mean_low3[i]
# 
#   # Save sample identities for interpretation after filtering
#   dt$low3_samples[i] <-
#     paste(sample_labels[low3_idx], collapse = ",")
# 
#   dt$high3_samples[i] <-
#     paste(sample_labels[high3_idx], collapse = ",")
# }

# rm(dt_plot)
# Retain VMRs exceeding the threshold
dt_plot <- dt[
  n_meth_nonNA_samples == length(sample_cols) &
    !is.na(extreme_3vs3_diff) &
    extreme_3vs3_diff >= DIFF_THRESH
]

message("\nVMRs before filtering: ", nrow(dt))
message("VMRs after blind 3-vs-3 filtering: ", nrow(dt_plot))

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
    "QC_counts_before_after_3vs3_filter.tsv"
  ),
  sep = "\t"
)

qc_difference <- rbindlist(
  list(
    
    dt[
      !is.na(extreme_3vs3_diff),
      .(
        stage  = "before_filter",
        min    = min(extreme_3vs3_diff),
        q25    = quantile(extreme_3vs3_diff, 0.25),
        median = median(extreme_3vs3_diff),
        mean   = mean(extreme_3vs3_diff),
        q75    = quantile(extreme_3vs3_diff, 0.75),
        max    = max(extreme_3vs3_diff)
      )
    ],
    
    dt_plot[
      ,
      .(
        stage  = "after_filter",
        min    = min(extreme_3vs3_diff),
        q25    = quantile(extreme_3vs3_diff, 0.25),
        median = median(extreme_3vs3_diff),
        mean   = mean(extreme_3vs3_diff),
        q75    = quantile(extreme_3vs3_diff, 0.75),
        max    = max(extreme_3vs3_diff)
      )
    ]
    
  )
)

fwrite(
  qc_difference,
  file.path(
    outdir,
    "QC_extreme_3vs3_difference_summary.tsv"
  ),
  sep = "\t"
)

# Save retained VMRs and identities of the extreme samples
fwrite(
  dt_plot,
  file.path(
    outdir,
    paste0(
      "vmr.strictUnsupervised.extreme3vs3.diffGE_",
      DIFF_THRESH,
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

color_palette <- inferno(n_colors)

breaks <- seq(
  0,
  1,
  length.out = n_colors + 1
)

# ==============================================================================
# Morph annotation
#
# Morph identity is used only for displaying the final heatmap.
# It does not influence filtering or clustering.
# ==============================================================================

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
    sample_order_tsv,
    row_order_tsv) {
  
  sub <- dt_plot[
    chromosome_type == chrom_type
  ]
  
  setorder(sub, chr, start, end)
  
  if (nrow(sub) < 2) {
    
    warning(
      "Not enough retained VMRs for ",
      chrom_type,
      ". Skipping."
    )
    
    return(invisible(NULL))
  }
  
  # --------------------------------------------------------------------------
  # Individual-sample matrix
  # --------------------------------------------------------------------------
  
  mat <- as.matrix(
    sub[, ..sample_cols]
  )
  
  storage.mode(mat) <- "double"
  
  rownames(mat) <- sub$row_id
  colnames(mat) <- sample_labels
  
  # --------------------------------------------------------------------------
  # Unsupervised VMR clustering
  # --------------------------------------------------------------------------
  
  hc_rows <- hclust(
    dist(
      mat,
      method = DIST_METHOD
    ),
    method = HCLUST_METHOD
  )
  
  # --------------------------------------------------------------------------
  # Unsupervised sample clustering
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
  
  ordered_samples <- colnames(mat)[hc_cols$order]
  
  sample_order <- data.table(
    Dendrogram_order = seq_along(ordered_samples),
    Sample = ordered_samples
  )
  
  # Morph added only after clustering
  sample_order[, Morph := get_morph(Sample)]
  
  fwrite(
    sample_order,
    file.path(
      outdir,
      sample_order_tsv
    ),
    sep = "\t"
  )
  
  # --------------------------------------------------------------------------
  # Save VMR dendrogram order
  # --------------------------------------------------------------------------
  
  ordered_rows <- rownames(mat)[hc_rows$order]
  
  row_order <- data.table(
    Heatmap_order = seq_along(ordered_rows),
    row_id = ordered_rows
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
        mean_low3,
        mean_high3,
        extreme_3vs3_diff,
        low3_samples,
        high3_samples
      )
    ],
    by = "row_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  setorder(row_order, Heatmap_order)
  
  fwrite(
    row_order,
    file.path(
      outdir,
      row_order_tsv
    ),
    sep = "\t"
  )
  
  # --------------------------------------------------------------------------
  # Heatmap
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
    
    cluster_rows = hc_rows,
    cluster_cols = hc_cols,
    
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
# Autosomes
# ==============================================================================

plot_unsupervised_heatmap(
  
  chrom_type = "Autosome",
  
  title = paste0(
    "Blind extreme 3-vs-3 VMR methylation — Autosomes\n",
    "Mean highest 3 − mean lowest 3 ≥ ",
    DIFF_THRESH
  ),
  
  pdf_name = paste0(
    "heatmap_AUTOSOME.strictUnsupervised.extreme3vs3.diffGE_",
    DIFF_THRESH,
    ".pdf"
  ),
  
  sample_order_tsv = paste0(
    "sample_order_AUTOSOME.extreme3vs3.diffGE_",
    DIFF_THRESH,
    ".tsv"
  ),
  
  row_order_tsv = paste0(
    "VMR_order_AUTOSOME.extreme3vs3.diffGE_",
    DIFF_THRESH,
    ".tsv"
  )
)

# ==============================================================================
# Chromosome 12
# ==============================================================================

plot_unsupervised_heatmap(
  
  chrom_type = "Sex_Ch",
  
  title = paste0(
    "Blind extreme 3-vs-3 VMR methylation — Chromosome 12\n",
    "Mean highest 3 − mean lowest 3 ≥ ",
    DIFF_THRESH
  ),
  
  pdf_name = paste0(
    "heatmap_CH12.strictUnsupervised.extreme3vs3.diffGE_",
    DIFF_THRESH,
    ".pdf"
  ),
  
  sample_order_tsv = paste0(
    "sample_order_CH12.extreme3vs3.diffGE_",
    DIFF_THRESH,
    ".tsv"
  ),
  
  row_order_tsv = paste0(
    "VMR_order_CH12.extreme3vs3.diffGE_",
    DIFF_THRESH,
    ".tsv"
  )
)

