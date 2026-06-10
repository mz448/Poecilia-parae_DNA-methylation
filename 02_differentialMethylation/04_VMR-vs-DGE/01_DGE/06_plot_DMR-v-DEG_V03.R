#!/usr/bin/env Rscript
# DATE:       2026-03-18
# AUTHOR:     MZF & ChatGPT
# SCRIPT:     04_plot_DMR_DE_feature_bubbles_auto.R
# VERSION:    01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Automatically read feature-specific DMR-vs-DE result files from one folder
#   and generate three bubble plots, one per analysis.
#
# QUESTION:
#   How do the strength and significance of the three DMR-vs-DE analyses vary
#   across pairwise comparisons and annotated gene features?
#
# APPROACH:
#   Parse all feature-specific result files from a directory using:
#     - the prefix before the first "." as the feature
#     - the suffix after the first "." as the analysis type
#   Then merge all results and generate one bubble plot per analysis with:
#     - x-axis = pairwise comparison
#     - y-axis = gene feature
#     - bubble size = -log10(FDR)
#     - bubble fill = main effect statistic for that analysis
#
# INPUTS:
#   Folder containing files named as:
#     <feature>.DMR_DE_anyOverlap_wilcox.tsv
#     <feature>.DMR_DE_doseResponse_abs.tsv
#     <feature>.DMR_DE_doseResponse_signed.tsv
#
# OUTPUTS:
#   - Bubbleplot_anyOverlap_wilcox.pdf / .png
#   - Bubbleplot_doseResponse_abs.pdf / .png
#   - Bubbleplot_doseResponse_signed.pdf / .png
#   - Bubbleplot_*_plotTable.tsv
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(forcats)
})

# =========================
# 0) USER SETTINGS
# =========================
input_dir <- "./DMR-vs-DEG"

comparison_order <- c(
  "i_vs_y_01",
  "p_vs_i_01",
  "y_vs_p_01",
  "f_vs_i",
  "f_vs_p",
  "f_vs_y"
)

# Set to NULL to keep alphabetical feature order
feature_order <- c(
  # "any",
  "u2000",
  "mRNA",
  "five_prime_UTR",
  "CDS",
  "geneBody_exon",
  "geneBody_intron",
  "three_prime_UTR",
  "d2000"
)

max_p_to_plot <- 0.05
min_p_floor   <- 1e-300

out_prefix <- "../plots/Bubbleplots_DMR-vs-DEG/Bubbleplot"

# If TRUE, include files whose prefix is "any"
include_any <- FALSE

# =========================
# 1) DETECT FILES
# =========================
all_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

# Keep only files that have a feature prefix followed by ".DMR_DE_..."
feature_files <- basename(all_files) %>%
  tibble(file = ., full_path = all_files) %>%
  filter(str_detect(file, "^[^.]+\\.DMR_DE_.*\\.tsv$")) %>%
  mutate(
    feature  = str_replace(file, "^([^.]+)\\..*$", "\\1"),
    analysis = str_replace(file, "^[^.]+\\.(DMR_DE_.*)\\.tsv$", "\\1")
  )

if (!include_any) {
  feature_files <- feature_files %>% filter(feature != "any")
}

if (nrow(feature_files) == 0) {
  stop("No feature-specific files found in: ", input_dir)
}

# Supported analyses
analysis_keep <- c(
  "DMR_DE_anyOverlap_wilcox",
  "DMR_DE_doseResponse_abs",
  "DMR_DE_doseResponse_signed"
)

feature_files <- feature_files %>%
  filter(analysis %in% analysis_keep)

if (nrow(feature_files) == 0) {
  stop("No supported analysis files found in: ", input_dir)
}

# =========================
# 2) HELPERS
# =========================
safe_read_tsv <- function(path) {
  read_tsv(path, show_col_types = FALSE)
}

# add_sig_score <- function(df, p_col) {
#   df %>%
#     mutate(
#       sig_p = .data[[p_col]],
#       sig_p = ifelse(is.na(sig_p), NA_real_, pmax(sig_p, min_p_for_plot)),
#       sig_score = -log10(sig_p)
#     )
# }
# 
# prep_factor_orders <- function(df, feature_order = NULL) {
#   if (is.null(feature_order)) {
#     feature_levels <- sort(unique(df$feature))
#   } else {
#     feature_levels <- feature_order[feature_order %in% unique(df$feature)]
#     feature_levels <- c(feature_levels, sort(setdiff(unique(df$feature), feature_levels)))
#   }
#   
#   df %>%
#     mutate(
#       comparison = factor(comparison, levels = comparison_order),
#       feature    = factor(feature, levels = rev(feature_levels))
#     )
# }
add_sig_score <- function(df, p_col, max_p_to_plot = 0.05, min_p_floor = 1e-300) {
  df %>%
    mutate(
      sig_p = .data[[p_col]],
      sig_score = case_when(
        is.na(sig_p)          ~ NA_real_,
        sig_p > max_p_to_plot ~ NA_real_,   # hide non-significant points
        TRUE                  ~ -log10(pmax(sig_p, min_p_floor))
      )
    )
}

prep_factor_orders <- function(df, feature_order = NULL, keep_only_listed = TRUE) {
  if (is.null(feature_order)) {
    feature_levels <- sort(unique(df$feature))
    df_out <- df
  } else {
    if (keep_only_listed) {
      df_out <- df %>% filter(feature %in% feature_order)
      feature_levels <- feature_order[feature_order %in% unique(df_out$feature)]
    } else {
      df_out <- df
      feature_levels <- feature_order[feature_order %in% unique(df_out$feature)]
      feature_levels <- c(feature_levels, sort(setdiff(unique(df_out$feature), feature_levels)))
    }
  }
  
  df_out %>%
    mutate(
      comparison = factor(comparison, levels = comparison_order),
      feature    = factor(feature, levels = rev(feature_levels))
    )
}


# plot_bubbles <- function(df, title_text, fill_label, size_label, out_stub,
#                          width = 8.8, height = 5.2) {
#   
#   p <- ggplot(df, aes(x = comparison, y = feature)) +
#     geom_point(
#       aes(size = sig_score, fill = effect),
#       shape = 21,
#       color = "black",
#       alpha = 0.9,
#       stroke = 0.3
#     ) +
#     scale_size_continuous(
#       name = size_label,
#       range = c(2.5, 14)
#     ) +
#     scale_fill_gradient2(
#       name = fill_label,
#       low = "#2166AC",
#       mid = "white",
#       high = "#B2182B",
#       midpoint = 0
#     ) +
#     labs(
#       title = title_text,
#       x = "Pairwise comparison",
#       y = "Gene feature"
#     ) +
#     theme_bw(base_size = 11) +
#     theme(
#       panel.grid.major = element_line(linewidth = 0.2, color = "grey85"),
#       panel.grid.minor = element_blank(),
#       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#       plot.title = element_text(face = "bold"),
#       legend.box = "vertical"
#     )
#   
#   ggsave(paste0(out_stub, ".pdf"), p, width = width, height = height)
#   ggsave(paste0(out_stub, ".png"), p, width = width, height = height, dpi = 300)
#   
#   invisible(p)
# }

# plot_bubbles <- function(df, title_text, fill_label, size_label, out_stub,
#                          width = 8.8, height = 5.2) {
#   
#   df_plot <- df %>% filter(!is.na(sig_score))
#   
#   p <- ggplot(df_plot, aes(x = comparison, y = feature)) +
#     geom_point(
#       aes(size = sig_score, fill = effect),
#       shape = 21,
#       color = "black",
#       alpha = 0.9,
#       stroke = 0.3
#     ) +
#     scale_size_continuous(
#       name = size_label,
#       range = c(2.5, 14)
#     ) +
#     scale_fill_gradient2(
#       name = fill_label,
#       low = "#2166AC",
#       mid = "white",
#       high = "#B2182B",
#       midpoint = 0
#     ) +
#     labs(
#       title = title_text,
#       x = "Pairwise comparison",
#       y = "Gene feature"
#     ) +
#     theme_bw(base_size = 11) +
#     theme(
#       panel.grid.major = element_line(linewidth = 0.2, color = "grey85"),
#       panel.grid.minor = element_blank(),
#       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#       plot.title = element_text(face = "bold"),
#       legend.box = "vertical"
#     )
#   
#   ggsave(paste0(out_stub, ".pdf"), p, width = width, height = height)
#   ggsave(paste0(out_stub, ".png"), p, width = width, height = height, dpi = 300)
#   
#   invisible(p)
# }

plot_bubbles <- function(df, title_text, fill_label, size_label, out_stub,
                         width = 8.8, height = 5.2) {
  
  p <- ggplot(df, aes(x = comparison, y = feature)) +
    geom_tile(fill = "white", color = "grey90", linewidth = 0.2) +
    geom_point(
      aes(size = sig_score, fill = effect),
      shape = 21,
      color = "black",
      alpha = 0.9,
      stroke = 0.3,
      na.rm = TRUE
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    scale_size_continuous(
      name = size_label,
      range = c(2.5, 14)
    ) +
    scale_fill_gradient2(
      name = fill_label,
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      na.value = NA
    ) +
    labs(
      title = title_text,
      x = "Pairwise comparison",
      y = "Gene feature"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(face = "bold"),
      legend.box = "vertical"
    )
  
  ggsave(paste0(out_stub, ".pdf"), p, width = width, height = height)
  # ggsave(paste0(out_stub, ".png"), p, width = width, height = height, dpi = 300)
  
  invisible(p)
}

read_feature_analysis <- function(feature_files_tbl, analysis_name) {
  feature_files_tbl %>%
    filter(analysis == analysis_name) %>%
    mutate(data = map(full_path, safe_read_tsv)) %>%
    select(feature, analysis, data) %>%
    unnest(data)
}

# =========================
# 3) READ EACH ANALYSIS SET
# =========================
all_any <- read_feature_analysis(feature_files, "DMR_DE_anyOverlap_wilcox")
all_abs <- read_feature_analysis(feature_files, "DMR_DE_doseResponse_abs")
all_signed <- read_feature_analysis(feature_files, "DMR_DE_doseResponse_signed")

# =========================
# 4) PREP PLOT TABLES
# =========================

# # ---- Analysis 1: any-overlap Wilcoxon
# # Fill = median |LFC| shift (DMR - noDMR)
# # Size = -log10(Wilcoxon FDR)
# plot_any <- all_any %>%
#   mutate(
#     effect = median_abs_l2fc_DMR - median_abs_l2fc_noDMR
#   ) %>%
#   add_sig_score("p_wilcox_fdr") %>%
#   prep_factor_orders(feature_order = feature_order)
# 
# # ---- Analysis 2: absolute dose-response
# # Fill = LM beta for weighted abs DMR difference
# # Size = -log10(LM FDR)
# plot_abs <- all_abs %>%
#   mutate(
#     effect = beta_weighted_abs_dmr_diff
#   ) %>%
#   add_sig_score("p_lm_fdr") %>%
#   prep_factor_orders(feature_order = feature_order)
# 
# # ---- Analysis 3: signed dose-response
# # Fill = LM beta for aligned signed DMR difference
# # Size = -log10(LM FDR)
# plot_signed <- all_signed %>%
#   mutate(
#     effect = beta_weighted_signed_dmr_diff_aligned
#   ) %>%
#   add_sig_score("p_lm_fdr") %>%
#   prep_factor_orders(feature_order = feature_order)

# Make the full grid for plotting
complete_plot_grid <- function(df, feature_order, comparison_order) {
  tidyr::expand_grid(
    comparison = comparison_order,
    feature = feature_order
  ) %>%
    left_join(df, by = c("comparison", "feature")) %>%
    mutate(
      comparison = factor(comparison, levels = comparison_order),
      feature    = factor(feature, levels = rev(feature_order))
    )
}


# ---- Analysis 1: any-overlap Wilcoxon
# Fill = median |LFC| shift (DMR - noDMR)
# Size = -log10(Wilcoxon FDR)
plot_any <- all_any %>%
  mutate(
    effect = median_abs_l2fc_DMR - median_abs_l2fc_noDMR
  ) %>%
  add_sig_score("p_wilcox_fdr", max_p_to_plot = max_p_to_plot, min_p_floor = min_p_floor) %>%
  filter(feature %in% feature_order) %>%
  complete_plot_grid(feature_order = feature_order, comparison_order = comparison_order)

# ---- Analysis 2: absolute dose-response
# Fill = LM beta for weighted abs DMR difference
# Size = -log10(LM FDR)

plot_abs <- all_abs %>%
  mutate(
    effect = beta_weighted_abs_dmr_diff
  ) %>%
  add_sig_score("p_lm_fdr", max_p_to_plot = max_p_to_plot, min_p_floor = min_p_floor) %>%
  filter(feature %in% feature_order) %>%
  complete_plot_grid(feature_order = feature_order, comparison_order = comparison_order)

# ---- Analysis 3: signed dose-response
# Fill = LM beta for aligned signed DMR difference
# Size = -log10(LM FDR)

plot_abs <- all_abs %>%
  mutate(
    effect = beta_weighted_abs_dmr_diff
  ) %>%
  add_sig_score("p_lm_fdr", max_p_to_plot = max_p_to_plot, min_p_floor = min_p_floor) %>%
  filter(feature %in% feature_order) %>%
  complete_plot_grid(feature_order = feature_order, comparison_order = comparison_order)




# =========================
# 5) WRITE PLOT TABLES
# =========================
write_tsv(plot_any,    paste0(out_prefix, "_anyOverlap_plotTable.tsv"))
write_tsv(plot_abs,    paste0(out_prefix, "_doseResponse_abs_plotTable.tsv"))
write_tsv(plot_signed, paste0(out_prefix, "_doseResponse_signed_plotTable.tsv"))

# =========================
# 6) DRAW PLOTS
# =========================
plot_bubbles(
  df = plot_any,
  title_text = "Any-overlap analysis",
  fill_label = "Median |LFC| shift\n(DMR - no DMR)",
  size_label = expression(-log[10]("Wilcoxon FDR")),
  out_stub = paste0(out_prefix, "_anyOverlap_wilcox")
)

plot_bubbles(
  df = plot_abs,
  title_text = "Absolute dose-response analysis",
  fill_label = "LM beta\n(weighted abs DMR diff)",
  size_label = expression(-log[10]("LM FDR")),
  out_stub = paste0(out_prefix, "_doseResponse_abs")
)

plot_bubbles(
  df = plot_signed,
  title_text = "Signed dose-response analysis",
  fill_label = "LM beta\n(aligned signed DMR diff)",
  size_label = expression(-log[10]("LM FDR")),
  out_stub = paste0(out_prefix, "_doseResponse_signed")
)

message("Done.")