#!/usr/bin/env Rscript
# DATE: 20260311
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#     Decide k-means cluster number (k) for AUTOSOMES vs CH12 using the
#     12-sample methylation matrix, but only after filtering consensus DMRs
#     with the same morph-mean threshold used in the heatmap:
#       - keep only regions where at least one pairwise morph comparison
#         shows |Δ methylation| >= DIFF_THRESH
#
#     Then evaluate k using:
#       1) Elbow (total within-cluster SS)
#       2) Mean silhouette width
#       3) Gap statistic (optional; can be slow)
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

library(data.table)
library(ggplot2)
library(cluster)

# -------------------------
# Defaults
# -------------------------
DIFF_THRESH <- 0.4  # keep regions with max pairwise abs diff >= this threshold (0..1)
infile <- "vmr_methylation/vmr.morph_and_sample_observed.tsv"
outdir <- paste0("../plots/",DIFF_THRESH,"/k_choice_sampleMatrix_filteredByMeans")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# To save the messages into a text file
logfile <- file.path(outdir, "k_choice_recommendations.txt")
sink(logfile, split = TRUE)
on.exit(sink(), add = TRUE)

morph_cols <- c("meth_female", "meth_immaculata", "meth_parae", "meth_yellow")


K_MIN <- 2
K_MAX <- 25
# K_MAX <- 12

SET_SEED <- 1
NSTART   <- 30
ITER_MAX <- 200

# For very large panels, evaluate k on a subsample to keep it fast
# SUBSAMPLE_N <- 5000   # set to Inf to disable subsampling
SUBSAMPLE_N <- Inf   # set to Inf to disable subsampling

# Gap statistic settings
DO_GAP <- TRUE
GAP_B  <- 30

# -------------------------
# Load + basic prep
# -------------------------
dt <- fread(infile)

if (!"chromosome_type" %in% names(dt)) {
  dt[, chromosome_type := ifelse(chr == "Parae_12", "Sex_Ch", "Autosome")]
}

stopifnot(all(morph_cols %in% names(dt)))

sample_cols <- c(
  sort(grep("^meth_pparfmem[0-9]+$", names(dt), value = TRUE)),
  sort(grep("^meth_pparimem[0-9]+$", names(dt), value = TRUE)),
  sort(grep("^meth_pparpmem[0-9]+$", names(dt), value = TRUE)),
  sort(grep("^meth_pparymem[0-9]+$", names(dt), value = TRUE))
)

if (length(sample_cols) == 0) {
  stop("No sample methylation columns found.")
}

# -------------------------
# Apply same filter used in heatmap
# -------------------------
dt[, `:=`(
  d_F_I = abs(meth_female     - meth_immaculata),
  d_F_P = abs(meth_female     - meth_parae),
  d_F_Y = abs(meth_female     - meth_yellow),
  d_I_P = abs(meth_immaculata - meth_parae),
  d_I_Y = abs(meth_immaculata - meth_yellow),
  d_P_Y = abs(meth_parae      - meth_yellow)
)]

dt[, max_pairwise_diff := do.call(
  pmax,
  c(.SD, na.rm = TRUE)
), .SDcols = c("d_F_I", "d_F_P", "d_F_Y", "d_I_P", "d_I_Y", "d_P_Y")]

dt[, n_meth_nonNA := rowSums(!is.na(.SD)), .SDcols = morph_cols]
dt[max_pairwise_diff == -Inf, max_pairwise_diff := NA_real_]

dt_plot <- dt[
  n_meth_nonNA >= 2 &
    !is.na(max_pairwise_diff) &
    max_pairwise_diff >= DIFF_THRESH
]

# For k-means evaluation on the 12-sample matrix, rows must be complete
dt_k <- dt_plot[complete.cases(dt_plot[, ..sample_cols])]

# save filtered input used for k choice
fwrite(
  dt_k,
  file.path(outdir, paste0("vmr_real.sampleMatrix.filteredByMeans.diffGE_", DIFF_THRESH, ".tsv")),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# QC summary
qc <- data.table(
  stage = c("before_filter", "after_morph_filter", "after_complete_sample_rows"),
  n = c(nrow(dt), nrow(dt_plot), nrow(dt_k))
)
fwrite(qc, file.path(outdir, "QC_filter_counts.tsv"), sep = "\t")

cat("\n=== FILTER SUMMARY ===\n")
cat("Rows before filter:", nrow(dt), "\n")
cat("Rows after morph-mean filter:", nrow(dt_plot), "\n")
cat("Rows retained for 12-sample k evaluation:", nrow(dt_k), "\n")
cat("DIFF_THRESH:", DIFF_THRESH, "\n")

# -------------------------
# Helpers
# -------------------------
as_numeric_matrix <- function(d) {
  m <- as.matrix(d[, ..sample_cols])
  storage.mode(m) <- "double"
  m
}

subsample_rows <- function(mat, n) {
  if (!is.finite(n) || n >= nrow(mat)) return(mat)
  set.seed(SET_SEED)
  mat[sample.int(nrow(mat), n), , drop = FALSE]
}

run_kmeans <- function(mat, k) {
  set.seed(SET_SEED)
  kmeans(
    mat,
    centers = k,
    nstart = NSTART,
    iter.max = ITER_MAX
  )
}

avg_silhouette <- function(mat, k) {
  km <- run_kmeans(mat, k)
  sil <- silhouette(km$cluster, dist(mat))
  mean(sil[, "sil_width"])
}

tot_withinss <- function(mat, k) {
  run_kmeans(mat, k)$tot.withinss
}

gap_stat <- function(mat, kmax, B) {
  set.seed(SET_SEED)
  clusGap(
    mat,
    FUN = function(x, k) {
      kmeans(
        x,
        centers = k,
        nstart = NSTART,
        iter.max = ITER_MAX
      )
    },
    K.max = kmax,
    B = B
  )
}

gap_1se_k <- function(gap_obj, k_min = 2) {
  tab <- as.data.table(gap_obj$Tab)
  tab[, k := .I]
  tab <- tab[k >= k_min]
  
  if (nrow(tab) == 1) return(tab$k[1])
  
  for (i in 1:(nrow(tab) - 1)) {
    if (tab$gap[i] >= tab$gap[i + 1] - tab$SE.sim[i + 1]) {
      return(tab$k[i])
    }
  }
  
  tab$k[which.max(tab$gap)]
}

evaluate_panel <- function(panel_name, dt_panel) {
  cat("\n============================\n")
  cat("PANEL:", panel_name, "\n")
  cat("============================\n")
  
  mat_full <- as_numeric_matrix(dt_panel)
  
  if (any(!is.finite(mat_full))) {
    stop(panel_name, ": found non-finite values in the 12-sample methylation matrix.")
  }
  
  mat <- subsample_rows(mat_full, SUBSAMPLE_N)
  
  n <- nrow(mat)
  if (n < 3) stop(panel_name, ": not enough rows (n < 3) to cluster.")
  
  k_max_use <- min(K_MAX, n - 1)
  k_seq <- K_MIN:k_max_use
  
  cat("Rows (full filtered panel):", nrow(mat_full), " | Rows (used for metrics):", nrow(mat), "\n")
  cat("Columns used for clustering:", ncol(mat), "\n")
  cat("K range:", K_MIN, "to", k_max_use, "\n")
  
  # compute metrics
  res <- data.table(k = k_seq)
  res[, tot_withinss := vapply(k, function(kk) tot_withinss(mat, kk), numeric(1))]
  res[, avg_sil      := vapply(k, function(kk) avg_silhouette(mat, kk), numeric(1))]
  
  # gap statistic
  if (DO_GAP && k_max_use >= 2) {
    g <- gap_stat(mat, kmax = k_max_use, B = GAP_B)
    gtab <- as.data.table(g$Tab)
    gtab[, k := .I]
    setnames(gtab, c("gap", "SE.sim"), c("gap", "gap_se"))
    res <- merge(res, gtab[, .(k, gap, gap_se)], by = "k", all.x = TRUE, sort = FALSE)
    k_gap <- gap_1se_k(g, k_min = K_MIN)
  } else {
    res[, `:=`(gap = NA_real_, gap_se = NA_real_)]
    k_gap <- NA_integer_
  }
  
  # recommendations
  k_sil <- res[which.max(avg_sil), k]
  
  w <- res$tot_withinss
  rel_drop <- c(NA, (w[-length(w)] - w[-1]) / w[-length(w)])
  res[, rel_drop := rel_drop]
  
  k_elbow <- res[rel_drop < 0.10 & !is.na(rel_drop), k][1]
  if (is.na(k_elbow)) k_elbow <- res[which.min(tot_withinss), k]
  
  cat("Suggested K (silhouette max):", k_sil, "\n")
  if (!is.na(k_gap)) cat("Suggested K (gap 1-SE):       ", k_gap, "\n")
  cat("Suggested K (elbow heuristic):", k_elbow, "\n")
  
  fwrite(res, file.path(outdir, paste0("k_metrics_", panel_name, ".tsv")), sep = "\t")
  
  p1 <- ggplot(res, aes(x = k, y = tot_withinss)) +
    geom_line() +
    geom_point() +
    labs(
      title = paste0(panel_name, ": Elbow (12-sample matrix, filtered by morph means)"),
      x = "k",
      y = "Total within-cluster SS"
    ) +
    theme_bw()
  
  p2 <- ggplot(res, aes(x = k, y = avg_sil)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = k_sil, linetype = "dashed") +
    labs(
      title = paste0(panel_name, ": Mean silhouette (12-sample matrix, filtered by morph means)"),
      x = "k",
      y = "Mean silhouette width"
    ) +
    theme_bw()
  
  if (DO_GAP) {
    p3 <- ggplot(res, aes(x = k, y = gap)) +
      geom_line() +
      geom_point() +
      geom_errorbar(aes(ymin = gap - gap_se, ymax = gap + gap_se), width = 0.15) +
      labs(
        title = paste0(panel_name, ": Gap statistic (12-sample matrix, filtered by morph means)"),
        x = "k",
        y = "Gap"
      ) +
      theme_bw()
  }
  
  pdf(file.path(outdir, paste0("k_choice_", panel_name, ".pdf")), width = 7, height = 9)
  print(p1)
  print(p2)
  if (DO_GAP) print(p3)
  dev.off()
  
  invisible(list(metrics = res, k_sil = k_sil, k_gap = k_gap, k_elbow = k_elbow))
}

# -------------------------
# Run per panel
# -------------------------
auto_res <- evaluate_panel("AUTOSOME", dt_k[chromosome_type == "Autosome"])
ch12_res <- evaluate_panel("CH12",     dt_k[chromosome_type == "Sex_Ch"])

cat("\n=== SUMMARY ===\n")
cat(
  "AUTOSOME: K_sil = ", auto_res$k_sil,
  "; K_elbow = ", auto_res$k_elbow,
  if (!is.na(auto_res$k_gap)) paste0("; K_gap = ", auto_res$k_gap) else "",
  "\n",
  sep = ""
)
cat(
  "CH12:     K_sil = ", ch12_res$k_sil,
  "; K_elbow = ", ch12_res$k_elbow,
  if (!is.na(ch12_res$k_gap)) paste0("; K_gap = ", ch12_res$k_gap) else "",
  "\n",
  sep = ""
)
cat("Outputs in:", outdir, "\n")
sink()