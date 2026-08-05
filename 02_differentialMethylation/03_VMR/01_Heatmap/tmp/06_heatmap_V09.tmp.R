#!/usr/bin/env Rscript

# DATE: 20260311
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Plot VMR methylation heatmaps using per-sample methylation values
#   (12 columns), while keeping the filtering and row clustering based on the
#   4 morph-mean methylation columns.
#
#   Steps:
#     - keep only consensus DMRs where at least one pairwise morph comparison
#       shows |Δ methylation| >= DIFF_THRESH (on 0–1 scale)
#     - cluster rows using the 4 morph-mean methylation columns
#     - plot the heatmap using the 12 sample methylation columns
#     - annotate clusters + insert gaps between clusters
#     - save filtered table and cluster assignments
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

library(data.table)
library(pheatmap)
library(viridis)
# install.packages("pheatmap")

# ---------- Parameters ----------
DIFF_THRESH <- 0.5   # keep regions with max pairwise abs diff >= this threshold (0..1)

##### ---------- Choose K ----------
    K_AUTO <- 8
    K_CH12 <- 7



# ---------- Inputs ----------

outdir <- paste0("../plots/",DIFF_THRESH,".real_morph")
infile <- "./dataset_Structure/real_morph_vmr_methylation/vmr.real_morph.morph_and_sample_observed.tsv"   # meth in 0..1
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)



dt <- fread(infile)


# ---------- Basic columns ----------
if (!"chromosome_type" %in% names(dt)) {
  dt[, chromosome_type := ifelse(chr == "Parae_12", "Sex_Ch", "Autosome")]
}

dt[, row_id := if ("cons_id" %in% names(dt)) cons_id else paste0(chr, ":", start, "-", end)]

morph_cols <- c("meth_female", "meth_immaculata", "meth_parae", "meth_yellow")

sample_cols <- sort(grep("^meth_ppar.*mem[0-9]+$", names(dt), value = TRUE))
sample_labels <- sub("^meth_", "", sample_cols)

required_cols <- c("chr", "start", "end", "row_id", "chromosome_type", morph_cols, sample_cols)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (length(sample_cols) == 0) {
  stop("No sample methylation columns were found.")
}

# ---------- Filter: max pairwise abs difference using morph means ----------
dt[, `:=`(
  d_F_I = abs(meth_female     - meth_immaculata),
  d_F_P = abs(meth_female     - meth_parae),
  d_F_Y = abs(meth_female     - meth_yellow),
  d_I_P = abs(meth_immaculata - meth_parae),
  d_I_Y = abs(meth_immaculata - meth_yellow),
  d_P_Y = abs(meth_parae      - meth_yellow)
)]

dt[, max_pairwise_diff := do.call(pmax, c(.SD, na.rm = TRUE)),
   .SDcols = c("d_F_I", "d_F_P", "d_F_Y", "d_I_P", "d_I_Y", "d_P_Y")]

# require at least 2 morph means present
dt[, n_meth_nonNA := rowSums(!is.na(.SD)), .SDcols = morph_cols]
dt[max_pairwise_diff == -Inf, max_pairwise_diff := NA_real_]

dt_plot <- dt[n_meth_nonNA >= 2 & !is.na(max_pairwise_diff) & max_pairwise_diff >= DIFF_THRESH]

# ---------- QC outputs ----------
qc <- rbindlist(list(
  dt[, .(stage = "before_filter", n = .N), by = chromosome_type],
  dt_plot[, .(stage = "after_filter",  n = .N), by = chromosome_type]
))
fwrite(qc,
       file.path(outdir, "QC_counts_before_after_filter.sampleHeatmap.tsv"),
       sep = "\t")

qc_diff <- dt[, .(
  stage  = "before_filter",
  min    = min(max_pairwise_diff, na.rm = TRUE),
  median = median(max_pairwise_diff, na.rm = TRUE),
  mean   = mean(max_pairwise_diff, na.rm = TRUE),
  max    = max(max_pairwise_diff, na.rm = TRUE)
)]

qc_diff2 <- dt_plot[, .(
  stage  = "after_filter",
  min    = min(max_pairwise_diff, na.rm = TRUE),
  median = median(max_pairwise_diff, na.rm = TRUE),
  mean   = mean(max_pairwise_diff, na.rm = TRUE),
  max    = max(max_pairwise_diff, na.rm = TRUE)
)]

fwrite(rbind(qc_diff, qc_diff2),
       file.path(outdir, "QC_maxPairwiseDiff_summary.sampleHeatmap.tsv"),
       sep = "\t")

fwrite(
  dt_plot,
  file.path(outdir, paste0("vmr.morph_and_sample_observed.filtered.diffGE_", DIFF_THRESH, ".tsv")),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ---------- Fixed color scale 0..1 ----------
n_colors <- 101
color_palette <- inferno(n_colors)
breaks <- seq(0, 1, length.out = length(color_palette) + 1)

# ---------- Heatmap function ----------
plot_heatmap_with_breaks <- function(chrom_type, k, title, pdf_name, clusters_tsv) {
  
  sub <- dt_plot[chromosome_type == chrom_type]
  setorder(sub, chr, start, end)
  
  if (nrow(sub) < 2) {
    warning("Not enough rows to plot for ", chrom_type, " after filtering (n < 2). Skipping.")
    return(invisible(NULL))
  }
  
  # Matrix used for row clustering (12 samples)
  mat_cluster <- as.matrix(sub[, ..sample_cols])
  
  # Matrix used for row clustering (4 morph means)
  # mat_cluster <- as.matrix(sub[, ..morph_cols])
  
  
  storage.mode(mat_cluster) <- "double"
  rownames(mat_cluster) <- sub$row_id
  mat_cluster[!is.finite(mat_cluster)] <- NA_real_

  # Matrix used for plotting (all samples)
  mat_plot <- as.matrix(sub[, ..sample_cols])
  storage.mode(mat_plot) <- "double"
  rownames(mat_plot) <- sub$row_id
  colnames(mat_plot) <- sample_labels
  mat_plot[!is.finite(mat_plot)] <- NA_real_
  
  hc <- hclust(dist(mat_cluster), method = "complete")
  cl <- cutree(hc, k = k)
  
  cl_ord <- cl[hc$order]
  gap_rows <- which(cl_ord[-length(cl_ord)] != cl_ord[-1])
  
  anno_row <- data.frame(Cluster = factor(cl, levels = sort(unique(cl))))
  rownames(anno_row) <- names(cl)
  
  anno_colors <- list(
    Cluster = setNames(viridis(length(unique(cl))), levels(anno_row$Cluster))
  )
  
  out_cl <- sub[, .(chr, start, end, row_id, max_pairwise_diff)]
  out_cl[, Cluster := cl[row_id]]
  fwrite(out_cl, file.path(outdir, clusters_tsv), sep = "\t")
  
  pdf(file.path(outdir, pdf_name), width = 6.5, height = 8)
  pheatmap(
    mat_plot,
    # cluster_rows = hc,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    gaps_row = gap_rows,
    annotation_row = anno_row,
    annotation_colors = anno_colors,
    treeheight_row = 30,
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

# ---------- Run heatmaps ----------
plot_heatmap_with_breaks(
  chrom_type = "Autosome",
  k = K_AUTO,
  title = paste0(
    "Methylation (0–1) in VMRs — Autosome\n",
    "filtered and clustered on morph means; plotted by sample (|Diff| ≥ ", DIFF_THRESH, ")"
  ),
  pdf_name = paste0("heatmap_AUTOSOME.sampleCols.clusterBreaks.diffGE_", DIFF_THRESH, ".pdf"),
  clusters_tsv = paste0("clusters_AUTOSOME.sampleHeatmap.diffGE_", DIFF_THRESH, ".tsv")
)

plot_heatmap_with_breaks(
  chrom_type = "Sex_Ch",
  k = K_CH12,
  title = paste0(
    "Methylation (0–1) in VMRs — Chromosome 12\n",
    "filtered and clustered on morph means; plotted by sample (|Diff| ≥ ", DIFF_THRESH, ")"
  ),
  pdf_name = paste0("heatmap_CH12.sampleCols.clusterBreaks.diffGE_", DIFF_THRESH, ".pdf"),
  clusters_tsv = paste0("clusters_CH12.sampleHeatmap.diffGE_", DIFF_THRESH, ".tsv")
)