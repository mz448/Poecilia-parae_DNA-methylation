#!/usr/bin/env Rscript
# DATE:       2025-09-28
# AUTHOR:     MZF & ChatGPT
# SCRIPT:     02_DEseq_02_plotting.R
# VERSION:    04
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Stage B — Load per-comparison dds objects and normalized counts saved by
#   Stage A (no DESeq() re-run). Generate plots to:
#     ../plots/<condition>/<comparison>/
#   Plots:
#     - PCA (using saved normalized_counts.tsv)
#     - heatmap3 (using VST from saved dds)
#     - pheatmap (top-50 by padj using saved normalized counts)
#     - Volcano (using results() from saved dds)
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# install.packages("heatmap3")
# install.packages("pheatmap")
suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(heatmap3)
  library(RColorBrewer)
  library(pheatmap)
  library(fs)
  library(purrr)
})

# ------------------------- 0) FIXED PARAMETERS -------------------------------
Wald_FDR <- 0.05

# ------------------------- 1) INPUT ROOT -------------------------------------
objects_root <- "deseq_objects"  # produced by Stage A

# ------------------------- 2) DISCOVER COMPARISONS ---------------------------
# Expect structure: deseq_objects/<condition>/<comparison>/dds.rds
ddsrds_paths <- dir_ls(objects_root, recurse = TRUE, type = "file", glob = "*.rds")
if (length(ddsrds_paths) == 0) quit(status = 0)

pairs_tbl <- tibble(dds_path = as.character(ddsrds_paths)) %>%
  mutate(
    condition = basename(dirname(dirname(dds_path))),
    comparison= basename(dirname(dds_path)),
    norm_path = file.path(dirname(dds_path), "normalized_counts.tsv"),
    out_dir   = file.path("../..", "plots", "DGE", condition, comparison)
  )

# ------------------------- 3) PLOTTING HELPERS -------------------------------
create_pca_plot_from_norm <- function(norm_counts_df, dds, outfile, title) {
  # norm_counts_df: gene_id + sample columns
  mat <- norm_counts_df %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  # ensure columns are the same samples (and order) as in dds
  samp <- rownames(colData(dds))
  mat  <- mat[, samp, drop = FALSE]
  
  # drop genes with zero variance (constant across samples)
  if (ncol(mat) < 2) return(invisible(NULL))
  vars <- apply(mat, 1, var, na.rm = TRUE)
  mat  <- mat[vars > 0, , drop = FALSE]
  
  # need at least 2 genes and 2 samples for PCA
  if (nrow(mat) < 2 || ncol(mat) < 2) return(invisible(NULL))
  
  # PCA on samples (variables = genes)
  pcs <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  
  pc_df <- as.data.frame(pcs$x[, 1:2, drop = FALSE])
  pc_df$Sample <- rownames(pc_df)
  
  meta <- as.data.frame(colData(dds))
  meta$Sample <- rownames(meta)  # overwrite or create a single 'Sample' column
  
  plot_df <- dplyr::left_join(pc_df, meta, by = "Sample")
  var_exp <- round((pcs$sdev[1:2]^2 / sum(pcs$sdev^2)) * 100, 2)
  
  png(filename = outfile, units = "in", width = 10, height = 8, res = 300)
  print(
    ggplot(plot_df, aes(PC1, PC2, color = Group)) +
      geom_point(size = 3) +
      labs(
        title = title,
        x = paste0("PC1 (", var_exp[1], "%)"),
        y = paste0("PC2 (", var_exp[2], "%)")
      ) +
      theme_classic(base_size = 14)
  )
  dev.off()
}

# --- HARDENED heatmap3 helper: safe for 0/1-gene and missing dimnames cases ---
create_heatmap3_from_dds <- function(dds, outfile_png, title,
                                     padj_threshold = 0.05,
                                     colors = colorRampPalette(c("blue","black","red"))(15)) {
  # Results / significant genes
  res <- results(dds, alpha = padj_threshold)
  res <- na.omit(res)
  sig_genes <- rownames(res)[res$padj < padj_threshold]
  
  message(sprintf("[heatmap3] %s | padj<=%.3f | sig genes: %d",
                  title, padj_threshold, length(sig_genes)))
  
  # Nothing to plot
  if (length(sig_genes) == 0) {
    message("[heatmap3] Skipping: no significant genes.")
    return(invisible(NULL))
  }
  
  # Keep only genes present
  sig_genes <- intersect(sig_genes, rownames(dds))
  if (length(sig_genes) == 0) {
    message("[heatmap3] Skipping: sig genes not found in dds.")
    return(invisible(NULL))
  }
  
  # VST (may return vector if only 1 gene)
  vst_obj <- varianceStabilizingTransformation(dds[sig_genes, ])
  mat <- assay(vst_obj)
  
  # Desired sample order
  samp <- rownames(colData(dds))
  
  # Coerce to proper 2D matrix with dimnames
  if (is.null(dim(mat))) {
    # 1 gene × N samples vector -> make 1xN matrix with explicit dimnames
    mat <- matrix(mat, nrow = 1, ncol = length(samp),
                  dimnames = list(sig_genes[1], samp))
  } else {
    # Ensure dimnames exist
    if (is.null(rownames(mat))) rownames(mat) <- sig_genes
    if (is.null(colnames(mat))) colnames(mat) <- colnames(vst_obj)
    # Enforce sample order to match dds
    mat <- mat[, samp, drop = FALSE]
  }
  
  # Guard for empty/degenerate matrix
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    message("[heatmap3] Skipping: degenerate matrix (0 rows or 0 cols).")
    return(invisible(NULL))
  }
  
  png(filename = outfile_png, units = "in", width = 8, height = 8, res = 300)
  heatmap3(mat,
           method = "complete", Rowv = TRUE, Colv = NA,
           col = colors, scale = "row", labRow = NA, showRowDendro = TRUE)
  title(main = title)
  dev.off()
}




create_pheatmap_from_norm <- function(norm_counts_df, dds, outfile_png, title,
                                      res_object, padj_cutoff = 0.05, top_n_genes = 50) {
  res_df <- as.data.frame(na.omit(as.data.frame(res_object)))
  sig_ids <- rownames(res_df)[res_df$padj <= padj_cutoff]
  if (length(sig_ids) == 0) return(invisible(NULL))
  
  mat <- norm_counts_df %>%
    column_to_rownames("gene_id") %>%
    as.matrix()
  
  # keep only samples present in dds and match order
  samp <- rownames(colData(dds))
  mat  <- mat[, samp, drop = FALSE]
  
  # keep only significant genes available in matrix
  sig_ids <- intersect(sig_ids, rownames(mat))
  if (length(sig_ids) == 0) return(invisible(NULL))
  
  # order by padj and select top N
  ord_ids <- sig_ids[order(res_df[sig_ids, "padj"])]
  top_ids <- head(ord_ids, top_n_genes)
  
  # drop rows with zero variance to avoid NaNs in scale()
  sds <- apply(mat[top_ids, , drop = FALSE], 1, sd, na.rm = TRUE)
  top_ids <- top_ids[sds > 0]
  if (length(top_ids) == 0) return(invisible(NULL))
  
  mat_scaled <- t(scale(t(mat[top_ids, , drop = FALSE])))
  
  ann_col <- as.data.frame(colData(dds)[, "Group", drop = FALSE])
  
  png(filename = outfile_png, units = "in", width = 10, height = 8, res = 300)
  pheatmap(mat_scaled,
           cluster_rows = TRUE, cluster_cols = TRUE,
           annotation_col = ann_col,
           show_rownames = TRUE, show_colnames = FALSE,
           color = colorRampPalette(brewer.pal(9, "RdBu"))(100),
           main = title)
  dev.off()
}

create_volcano_simple <- function(res_object, outfile_png, title,
                                  padj_cutoff = 0.05, lfc_cutoff = 0,
                                  xlim_range = c(-5, 5), ylim_range = c(0, 10),
                                  up_label = "Up (alt)", down_label = "Up (ref)") {
  res_df <- as.data.frame(na.omit(as.data.frame(res_object)))
  if (!all(c("log2FoldChange","padj") %in% colnames(res_df))) return(invisible(NULL))
  
  res_df$log10Pval <- -log10(res_df$padj)
  res_df$Significance <- "Not Significant"
  res_df$Significance[res_df$log2FoldChange >  lfc_cutoff & res_df$padj < padj_cutoff] <- up_label
  res_df$Significance[res_df$log2FoldChange < -lfc_cutoff & res_df$padj < padj_cutoff] <- down_label
  
  png(filename = outfile_png, units = "in", width = 6, height = 6, res = 300)
  print(
    ggplot(res_df, aes(x = log2FoldChange, y = log10Pval, color = Significance)) +
      geom_point(alpha = 0.8, size = 2) +
      labs(title = title,
           x = expression("Log"[2] ~ "Fold Change"),
           y = expression("-Log"[10] ~ "padj")) +
      theme_classic(base_size = 14) +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
      coord_cartesian(xlim = xlim_range, ylim = ylim_range)
  )
  dev.off()
}

# ------------------------- 4) MAIN LOOP --------------------------------------
walk2(pairs_tbl$dds_path, seq_len(nrow(pairs_tbl)), function(dds_path, idx) {
  cond_id <- pairs_tbl$condition[idx]
  comp_id <- pairs_tbl$comparison[idx]
  norm_path <- pairs_tbl$norm_path[idx]
  out_dir <- pairs_tbl$out_dir[idx]
  dir_create(out_dir, recurse = TRUE)
  
  message(sprintf("\n[loop] index=%d | condition=%s | comparison=%s", idx, cond_id, comp_id))
  
  # Load saved objects (NO DESeq() re-run)
  dds <- readRDS(dds_path)
  norm_df <- read_tsv(norm_path, show_col_types = FALSE)
  
  # Results extracted from fitted dds (cheap; not re-fitting)
  res <- results(dds, alpha = Wald_FDR)
  res <- res[order(res$padj), ]
  
  title_tag <- paste0(cond_id, " — ", comp_id)
  up_label   <- paste0("Up in ", levels(dds$Group)[2])
  down_label <- paste0("Up in ", levels(dds$Group)[1])
  
  # PCA
  message("[PCA] start")
  tryCatch({
    create_pca_plot_from_norm(
      norm_counts_df = norm_df,
      dds = dds,
      outfile = file.path(out_dir, "PCA.png"),
      title = paste0("PCA — ", title_tag)
    )
    message("[PCA] done")
  }, error = function(e) {
    message(sprintf("[PCA] ERROR: %s", conditionMessage(e)))
  })
  
  # heatmap3
  message("[heatmap3] start")
  tryCatch({
    create_heatmap3_from_dds(
      dds = dds,
      outfile_png = file.path(out_dir, "heatmap3_sig.png"),
      title = paste0("heatmap3 — ", title_tag),
      padj_threshold = Wald_FDR
    )
    message("[heatmap3] done")
  }, error = function(e) {
    message(sprintf("[heatmap3] ERROR: %s", conditionMessage(e)))
  })
  
  # pheatmap
  message("[pheatmap] start")
  tryCatch({
    create_pheatmap_from_norm(
      norm_counts_df = norm_df,
      dds = dds,
      outfile_png = file.path(out_dir, "pheatmap_top50.png"),
      title = paste0("Top-50 DEGs — ", title_tag),
      res_object = res,
      padj_cutoff = Wald_FDR,
      top_n_genes = 50
    )
    message("[pheatmap] done")
  }, error = function(e) {
    message(sprintf("[pheatmap] ERROR: %s", conditionMessage(e)))
  })
  
  # Volcano
  message("[volcano] start")
  tryCatch({
    create_volcano_simple(
      res_object = res,
      outfile_png = file.path(out_dir, "volcano.png"),
      title = paste0("Volcano — ", title_tag),
      padj_cutoff = Wald_FDR,
      lfc_cutoff = 0,
      up_label = up_label,
      down_label = down_label
    )
    message("[volcano] done")
  }, error = function(e) {
    message(sprintf("[volcano] ERROR: %s", conditionMessage(e)))
  })
})
