#!/usr/bin/env Rscript
# DATE:       2025-10-01
# AUTHOR:     MZF
# SCRIPT:     04_plot_summary_ALL-IN-ONE_V01.R
# VERSION:    01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Run ALL plots (bubble, heatmap, scatter) for ALL datasets across ALL
#   combinations of {mode_keep ∈ [abs, signed]} × {r_to_implement ∈ [r, r2]}.
#   Each dataset writes into its own canonical output folders.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(ggplot2); library(stringr); library(ggrepel)
})

# =============================== CONFIG ======================================

pthresh <- 0.05  # significance threshold for (adjusted) p-values

# Datasets to process (summary_file + base plot dir)
datasets <- list(
  list(
    name         = "log2FC",
    summary_file = "log2FC_correlation_summary.tsv",
    base_outdir  = file.path("..","plots","log2FC","02_correlations")
  ),
  list(
    name         = "DElog2FC_x_log10pval",
    summary_file = "DElog2FCx-log10pval_correlation_summary.tsv",
    base_outdir  = file.path("..","plots","DElog2FC_x_log10pval","02_correlations")
  ),
  list(
    name         = "DEstat",
    summary_file = "DEstat_correlation_summary.tsv",
    base_outdir  = file.path("..","plots","DEstat","02_correlations")
  )
)

# Optional filters (set to NULL to keep all)
# selected_features    <- NULL
selected_features    <- c("five_prime_UTR","three_prime_UTR","geneBody_exon","geneBody_intron")

selected_conditions  <- NULL
selected_comparisons <- NULL

# ============================== HELPERS ======================================

ensure_dirs <- function(...) {
  d <- file.path(...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

ensure_r2_cols <- function(df){
  if (!"pearson_r2"    %in% names(df) && "pearson_r"    %in% names(df))    df$pearson_r2    <- df$pearson_r^2
  if (!"spearman_rho2" %in% names(df) && "spearman_rho" %in% names(df)) df$spearman_rho2 <- df$spearman_rho^2
  df
}

pick_metric_cols <- function(mode = c("r","r2")){
  mode <- match.arg(mode)
  if (mode == "r") {
    return(list(
      pear_col = "pearson_r",     rho_col  = "spearman_rho",
      pear_lab = "r",             rho_lab  = "\u03C1",
      signed   = TRUE,
      xint     = 0               # for scatter: vertical reference at 0
    ))
  } else {
    return(list(
      pear_col = "pearson_r2",    rho_col  = "spearman_rho2",
      pear_lab = expression(r^2), rho_lab  = expression(rho^2),
      signed   = FALSE,
      xint     = NA              # nonnegative; no vertical 0 line
    ))
  }
}

# ------------------------------ BUBBLES --------------------------------------

make_bubble_sig_only <- function(data, eff_col, nlogp_col, padj_col, title, outfile, eff_legend = "r", signed = TRUE,
                                 feat_order, comp_order, pthresh = 0.05) {
  base_grid <- expand.grid(
    feature = feat_order, comp_key = comp_order, stringsAsFactors = FALSE
  ) %>%
    mutate(feature = factor(feature, levels = rev(feat_order)),
           comp_key = factor(comp_key, levels = comp_order))
  
  dd <- data %>%
    transmute(
      feature, comp_key,
      effect = .data[[eff_col]],
      nlogp  = .data[[nlogp_col]],
      padj   = .data[[padj_col]]
    ) %>%
    mutate(
      feature = factor(feature, levels = rev(feat_order)),
      comp_key = factor(comp_key, levels = comp_order),
      sig = !is.na(padj) & padj < pthresh
    )
  
  dd_sig <- dd %>% filter(sig)
  
  col_scale <- {
    if (nrow(dd_sig) > 0) {
      if (signed) {
        lim <- max(abs(dd_sig$effect), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1e-8
        scale_color_gradient2(name = eff_legend, low = "#2c7bb6", mid = "white", high = "#d7191c",
                              midpoint = 0, limits = c(-lim, lim), oob = scales::squish)
      } else {
        lim <- max(dd_sig$effect, na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1e-8
        scale_color_gradient(name = eff_legend, limits = c(0, lim), oob = scales::squish)
      }
    } else {
      if (signed) {
        lim_all <- max(abs(dd$effect), na.rm = TRUE); if (!is.finite(lim_all) || lim_all == 0) lim_all <- 1e-8
        scale_color_gradient2(name = eff_legend, low = "#2c7bb6", mid = "white", high = "#d7191c",
                              midpoint = 0, limits = c(-lim_all, lim_all), oob = scales::squish)
      } else {
        lim_all <- max(dd$effect, na.rm = TRUE); if (!is.finite(lim_all) || lim_all == 0) lim_all <- 1e-8
        scale_color_gradient(name = eff_legend, limits = c(0, lim_all), oob = scales::squish)
      }
    }
  }
  
  comp_labs <- setNames(gsub("\\|", "\n", levels(base_grid$comp_key)), levels(base_grid$comp_key))
  
  p <- ggplot() +
    geom_tile(
      data = base_grid,
      aes(x = comp_key, y = feature),
      fill = NA, color = "grey90", linewidth = 0.25
    ) +
    geom_point(
      data = dd_sig,
      aes(x = comp_key, y = feature, size = nlogp, color = effect),
      alpha = 0.9, na.rm = TRUE
    ) +
    scale_x_discrete(labels = comp_labs) +
    col_scale +
    scale_size_continuous(
      name = expression(-log[10](padj)),
      range = c(1, 6),
      breaks = c(0, 1.301, 2, 3, 5),
      labels = c("0", "0.05", "0.01", "0.001", "1e-5")
    ) +
    labs(
      title = paste0(title, " — bubbles shown only if padj<", pthresh),
      x = "condition | comparison", y = "feature"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold"),
          legend.box = "vertical")
  
  ggsave(outfile, p, width = 10, height = max(3.5, length(feat_order) * 0.35), dpi = 300)
  if (nrow(dd_sig) == 0) {
    message("No significant points (padj < ", pthresh, ") — wrote empty grid: ", outfile)
  } else {
    message("Wrote: ", outfile)
  }
}

# ------------------------------ HEATMAPS -------------------------------------

make_heatmap <- function(data, eff_col, padj_col, title, outfile, legend_title="r", signed=TRUE,
                         feat_order, comp_order, pthresh = 0.05) {
  dlong <- data %>%
    transmute(
      feature, comp_key,
      effect = .data[[eff_col]],
      padj   = .data[[padj_col]],
      sig    = !is.na(padj) & padj < pthresh
    ) %>%
    complete(feature = feat_order, comp_key = comp_order) %>%
    mutate(feature = factor(feature, levels = rev(feat_order)),
           comp_key = factor(comp_key, levels = comp_order))
  
  comp_labs <- setNames(gsub("\\|", "\n", levels(dlong$comp_key)), levels(dlong$comp_key))
  
  if (signed) {
    lim <- max(abs(dlong$effect), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1e-8
    fill_scale <- scale_fill_gradient2(
      name = legend_title, low = "#2c7bb6", mid = "white", high = "#d7191c",
      midpoint = 0, na.value = "grey95", limits = c(-lim, +lim)
    )
  } else {
    lim <- max(dlong$effect, na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1e-8
    fill_scale <- scale_fill_gradient(
      name = legend_title, na.value = "grey95", limits = c(0, lim)
    )
  }
  
  p <- ggplot(dlong, aes(x = comp_key, y = feature, fill = effect)) +
    geom_tile(color = "grey90", linewidth = 0.2, na.rm = FALSE) +
    geom_point(
      data = dplyr::filter(dlong, sig),
      aes(x = comp_key, y = feature),
      shape = 21, stroke = 0.2, size = 2, fill = "black", color = "black"
    ) +
    scale_x_discrete(labels = comp_labs) +
    fill_scale +
    labs(title = title, x = "condition | comparison", y = "feature") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(face = "bold"))
  
  ggsave(outfile, p, width = 10, height = max(3.5, length(feat_order) * 0.35), dpi = 300)
  message("Wrote: ", outfile)
}

# ------------------------------ SCATTERS -------------------------------------

base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

make_scatter_panel <- function(df, x_col, y_col, title, x_lab, outfile, xint = NA, pthresh = 0.05){
  pthr <- -log10(pthresh)
  p <- ggplot(
    df, aes(x = .data[[x_col]], y = .data[[y_col]],
            color = dmr_condition, shape = feature, label = dmr_comparison)
  ) +
    geom_hline(yintercept = pthr, linetype = "dashed", linewidth = 0.4) +
    { if (!is.na(xint)) geom_vline(xintercept = xint, linetype = "dotted", linewidth = 0.3) else NULL } +
    geom_point(alpha = 0.85, size = 2) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, min.segment.length = 0, show.legend = FALSE) +
    labs(title = title, x = x_lab, y = expression(-log[10](padj)),
         color = "condition", shape = "feature") +
    base_theme
  
  ggsave(outfile, p, width = 8, height = 6, dpi = 300)
  message("Wrote: ", outfile)
}

# --------------------------- LOAD + FILTER CORE ------------------------------

load_and_prepare <- function(path_tsv, mode_keep,
                             selected_features, selected_conditions, selected_comparisons){
  readr::read_tsv(path_tsv, show_col_types = FALSE) %>%
    mutate(
      feature  = as.character(feature),
      mode     = as.character(mode),
      comp_key = paste0(dmr_condition, "|", dmr_comparison),
      nlog10_pear  = -log10(pmax(pearson_padj,  .Machine$double.xmin)),
      nlog10_spear = -log10(pmax(spearman_padj, .Machine$double.xmin))
    ) %>%
    filter(mode == mode_keep) %>%
    { if (!is.null(selected_features))    filter(., feature %in% selected_features) else . } %>%
    { if (!is.null(selected_conditions))  filter(., dmr_condition %in% selected_conditions) else . } %>%
    { if (!is.null(selected_comparisons)) filter(., dmr_comparison %in% selected_comparisons) else . } %>%
    ensure_r2_cols()
}

orders_from_df <- function(df){
  feat_order <- if (!is.null(selected_features)) selected_features else sort(unique(df$feature))
  comp_order <- df %>%
    distinct(dmr_condition, dmr_comparison, comp_key) %>%
    arrange(dmr_condition, dmr_comparison) %>%
    pull(comp_key)
  list(feat_order = feat_order, comp_order = comp_order)
}

# ================================ DRIVER =====================================

for (ds in datasets) {
  message("=== DATASET: ", ds$name, " (", ds$summary_file, ") ===")
  # Output subfolders
  out_bubble  <- ensure_dirs(ds$base_outdir, "bubbles_padj")
  out_heatmap <- ensure_dirs(ds$base_outdir, "heatmaps_padj")
  out_scatter <- ensure_dirs(ds$base_outdir, "scatter_padj")
  
  for (mode_keep in c("abs","signed")) {
    message(" -- mode_keep: ", mode_keep)
    df_all <- load_and_prepare(
      ds$summary_file, mode_keep,
      selected_features, selected_conditions, selected_comparisons
    )
    
    if (nrow(df_all) == 0) {
      message("   [SKIP] No rows after filtering for mode=", mode_keep, " in ", ds$summary_file)
      next
    }
    
    ord <- orders_from_df(df_all)
    feat_order <- ord$feat_order
    comp_order <- ord$comp_order
    
    for (metric in c("r","r2")) {
      cfg <- pick_metric_cols(metric)
      message("   .. metric: ", metric)
      
      # -------------------- BUBBLES --------------------
      make_bubble_sig_only(
        data      = df_all,
        eff_col   = cfg$pear_col,
        nlogp_col = "nlog10_pear",
        padj_col  = "pearson_padj",
        title     = paste0(ds$name, " — Pearson bubble (mode=", mode_keep, ", metric=", metric, ")"),
        outfile   = file.path(out_bubble, paste0("bubble_Pearson_metric-", metric, "_mode-", mode_keep, ".pdf")),
        eff_legend = cfg$pear_lab,
        signed    = cfg$signed,
        feat_order = feat_order,
        comp_order = comp_order,
        pthresh = pthresh
      )
      
      make_bubble_sig_only(
        data      = df_all,
        eff_col   = cfg$rho_col,
        nlogp_col = "nlog10_spear",
        padj_col  = "spearman_padj",
        title     = paste0(ds$name, " — Spearman bubble (mode=", mode_keep, ", metric=", metric, ")"),
        outfile   = file.path(out_bubble, paste0("bubble_Spearman_metric-", metric, "_mode-", mode_keep, ".pdf")),
        eff_legend = cfg$rho_lab,
        signed    = cfg$signed,
        feat_order = feat_order,
        comp_order = comp_order,
        pthresh = pthresh
      )
      
      # -------------------- HEATMAPS -------------------
      make_heatmap(
        data = df_all, eff_col = cfg$pear_col, padj_col = "pearson_padj",
        title = paste0(ds$name, " — Pearson heatmap (mode=", mode_keep, ", metric=", metric, ") — dots: padj<", pthresh),
        outfile = file.path(out_heatmap, paste0("heatmap_Pearson_metric-", metric, "_mode-", mode_keep, "_padj.pdf")),
        legend_title = cfg$pear_lab, signed = cfg$signed,
        feat_order = feat_order, comp_order = comp_order, pthresh = pthresh
      )
      
      make_heatmap(
        data = df_all, eff_col = cfg$rho_col, padj_col = "spearman_padj",
        title = paste0(ds$name, " — Spearman heatmap (mode=", mode_keep, ", metric=", metric, ") — dots: padj<", pthresh),
        outfile = file.path(out_heatmap, paste0("heatmap_Spearman_metric-", metric, "_mode-", mode_keep, "_padj.pdf")),
        legend_title = cfg$rho_lab, signed = cfg$signed,
        feat_order = feat_order, comp_order = comp_order, pthresh = pthresh
      )
      
      # -------------------- SCATTERS -------------------
      make_scatter_panel(
        df_all, x_col = cfg$pear_col, y_col = "nlog10_pear",
        title = paste0(ds$name, " — Pearson scatter (mode=", mode_keep, ", metric=", metric, ")"),
        x_lab = cfg$pear_lab,
        outfile = file.path(out_scatter, paste0("scatter_Pearson_metric-", metric, "_mode-", mode_keep, "_padj.pdf")),
        xint = cfg$xint, pthresh = pthresh
      )
      
      make_scatter_panel(
        df_all, x_col = cfg$rho_col, y_col = "nlog10_spear",
        title = paste0(ds$name, " — Spearman scatter (mode=", mode_keep, ", metric=", metric, ")"),
        x_lab = cfg$rho_lab,
        outfile = file.path(out_scatter, paste0("scatter_Spearman_metric-", metric, "_mode-", mode_keep, "_padj.pdf")),
        xint = cfg$xint, pthresh = pthresh
      )
    } # metric
  }   # mode_keep
}     # dataset

message("ALL DONE.")
