#!/usr/bin/env Rscript
# DATE:       2025-09-30
# AUTHOR:     MZF & ChatGPT
# SCRIPT:     03_plot_02_DMRdiff_vs_DElog2FCx-log10pval_V04.R
# VERSION:    04
# GOAL: Correlate DiffMethylation vs (log2FC * -log10(padj)), signed & absolute.
#       V04 — adds r^2 exports + robust saving/plotting guards (safe_save, cor_safe)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2); library(stringr); library(tidyr)
})

# ---------------------- PARAMETER BLOCK (edit here) ---------------------------
compiled_file <- "compiled_intersect_DEGs_AUTO_CPMfilter.tsv"
dat <- readr::read_tsv(compiled_file, na = c("", "NA"), show_col_types = FALSE) %>%
  dplyr::mutate(feature = stringr::str_trim(as.character(feature)))

# All features present in the compiled table
feature_list <- unique(dat$feature)

# Each entry: c(dmr_condition_sel, dmr_comparison_sel, deg_l2fc_col, deg_padj_col, deg_Wald-stat)
column_sets <- list(
  c("real_sex",       "f_vs_p",     "f_vs_p_l2fc",     "f_vs_p_padj",     "f_vs_p_stat"),
  c("real_sex",       "f_vs_y",     "f_vs_y_l2fc",     "f_vs_y_padj",     "f_vs_y_stat"),
  c("real_sex",       "f_vs_i",     "f_vs_i_l2fc",     "f_vs_i_padj",     "f_vs_i_stat"),
  c("real_morph",     "p_vs_i_01",  "p_vs_i_01_l2fc",  "p_vs_i_01_padj",  "p_vs_i_01_stat"),
  c("real_morph",     "y_vs_p_01",  "y_vs_p_01_l2fc",  "y_vs_p_01_padj",  "y_vs_p_01_stat"),
  c("real_morph",     "i_vs_y_01",  "i_vs_y_01_l2fc",  "i_vs_y_01_padj",  "i_vs_y_01_stat"),
  c("shuffled_sex",   "f_vs_sM_01", "f_vs_sM_01_l2fc", "f_vs_sM_01_padj", "f_vs_sM_01_stat"),
  c("shuffled_sex",   "f_vs_sM_02", "f_vs_sM_02_l2fc", "f_vs_sM_02_padj", "f_vs_sM_02_stat"),
  c("shuffled_sex",   "f_vs_sM_03", "f_vs_sM_03_l2fc", "f_vs_sM_03_padj", "f_vs_sM_03_stat"),
  c("shuffled_sex",   "f_vs_sM_04", "f_vs_sM_04_l2fc", "f_vs_sM_04_padj", "f_vs_sM_04_stat"),
  c("shuffled_sex",   "f_vs_sM_05", "f_vs_sM_05_l2fc", "f_vs_sM_05_padj", "f_vs_sM_05_stat"),
  c("shuffled_sex",   "f_vs_sM_06", "f_vs_sM_06_l2fc", "f_vs_sM_06_padj", "f_vs_sM_06_stat"),
  c("shuffled_morph", "sM_vs_sM_01","sM_vs_sM_01_l2fc","sM_vs_sM_01_padj","sM_vs_sM_01_stat"),
  c("shuffled_morph", "sM_vs_sM_02","sM_vs_sM_02_l2fc","sM_vs_sM_02_padj","sM_vs_sM_02_stat"),
  c("shuffled_morph", "sM_vs_sM_03","sM_vs_sM_03_l2fc","sM_vs_sM_03_padj","sM_vs_sM_03_stat"),
  c("shuffled_all",   "sA_vs_sA_01","sA_vs_sA_01_l2fc","sA_vs_sA_01_padj","sA_vs_sA_01_stat"),
  c("shuffled_all",   "sA_vs_sA_02","sA_vs_sA_02_l2fc","sA_vs_sA_02_padj","sA_vs_sA_02_stat"),
  c("shuffled_all",   "sA_vs_sA_03","sA_vs_sA_03_l2fc","sA_vs_sA_03_padj","sA_vs_sA_03_stat")
)

# Output summary
summary_dir  <- file.path("..", "data")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(summary_dir, "DElog2FCx-log10pval_correlation_summary.tsv")

# Plot params
fdr_alpha   <- 0.05
point_alpha <- 0.30
point_size  <- 1.0

# Collector
all_results <- list()
# -----------------------------------------------------------------------------


# ---------------------------- Helpers ----------------------------------------

# Safe ggsave: never pass a NULL path
safe_save <- function(outfile, plot, width, height, dpi = 300) {
  d <- dirname(outfile); if (!nzchar(d)) d <- "."
  f <- basename(outfile)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename = f, plot = plot, path = d, width = width, height = height, dpi = dpi)
}

# Correlation wrapper that won’t error on tiny/degenerate slices
# Correlation wrapper that won’t error on tiny/degenerate slices
cor_safe <- function(x, y, method = c("pearson","spearman")) {
  method <- match.arg(method)
  out <- tryCatch(
    suppressWarnings(stats::cor.test(x, y, method = method)),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(list(estimate = NA_real_, p.value = NA_real_))
  }
  est <- tryCatch(unname(out$estimate[[1]]), error = function(e) NA_real_)
  p   <- tryCatch(out$p.value,              error = function(e) NA_real_)
  list(estimate = est, p.value = p)
}


# Plot + stats with guards; uses de_comb for y
make_plot <- function(df, xlab, ylab, title, outfile) {
  pear  <- cor_safe(df$dmr_diff, df$de_comb, method = "pearson")
  spear <- cor_safe(df$dmr_diff, df$de_comb, method = "spearman")
  
  message(title)
  message(sprintf("  Pearson:  r = %.3f,  p = %.3g", pear$estimate, pear$p.value))
  message(sprintf("  Spearman: rho = %.3f, p = %.3g", spear$estimate, spear$p.value))
  
  n_ok   <- nrow(df) >= 3
  var_ok <- isTRUE(stats::sd(df$dmr_diff) > 0) && isTRUE(stats::sd(df$de_comb) > 0)
  add_smooth <- n_ok && var_ok
  
  xr <- range(df$dmr_diff, na.rm = TRUE)
  yr <- range(df$de_comb,  na.rm = TRUE)
  x_annot <- xr[1] + 0.02 * diff(xr)
  y_annot <- yr[2] - 0.02 * diff(yr)
  annot_text <- sprintf("Pearson r = %.3f (p=%.3g)\nSpearman rho = %.3f (p=%.3g)",
                        pear$estimate, pear$p.value, spear$estimate, spear$p.value)
  
  p <- ggplot(df, aes(x = dmr_diff, y = de_comb, color = sig)) +
    geom_hline(yintercept = 0, linewidth = 0.3, linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.3, linetype = "dashed") +
    geom_point(alpha = point_alpha, size = point_size) +
    { if (add_smooth) geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, color = "blue") } +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red"),
                       labels = c(`FALSE` = paste0("padj >= ", fdr_alpha),
                                  `TRUE`  = paste0("padj < ",  fdr_alpha)),
                       name = "DE significance") +
    annotate("text", x = x_annot, y = y_annot, hjust = 0, vjust = 1,
             label = annot_text, size = 3.3) +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  
  print(p)
  if (!is.null(outfile) && nzchar(outfile)) {
    safe_save(outfile, p, width = 6.5, height = 5.0, dpi = 300)
    message("Saved plot: ", outfile)
  }
  
  list(pear = pear, spear = spear)
}
# -----------------------------------------------------------------------------


for (set in column_sets) {
  # Parameters for this comparison
  dmr_condition_sel  <- set[[1]]
  dmr_comparison_sel <- set[[2]]
  deg_l2fc_col       <- set[[3]]
  deg_padj_col       <- set[[4]]
  
  for (feat in feature_list) {
    
    # Output folders & filenames
    out_dir <- file.path("..", "plots", "DElog2FC_x_log10pval", dmr_condition_sel, feat)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_signed_pdf <- file.path(out_dir, paste0("DMRdiff_vs_DElog2FCx-log10pval_SIGNED_", dmr_comparison_sel, ".pdf"))
    out_abs_pdf    <- file.path(out_dir, paste0("DMRdiff_vs_DElog2FCx-log10pval_ABS_"   , dmr_comparison_sel, ".pdf"))
    
    # Slice & compute combined DE metric
    df_base <- dat %>%
      filter(
        .data[["dmr_condition"]]  == dmr_condition_sel,
        .data[["dmr_comparison"]] == dmr_comparison_sel,
        .data[["feature"]]        == feat
      ) %>%
      transmute(
        gene,
        feature  = .data[["feature"]],
        dmr_diff = .data[["dmr_diff"]],
        de_comb  = .data[[deg_l2fc_col]] * (-log10(pmax(.data[[deg_padj_col]], .Machine$double.xmin))),
        de_padj  = .data[[deg_padj_col]],
        sig      = !is.na(.data[[deg_padj_col]]) & .data[[deg_padj_col]] < fdr_alpha
      ) %>%
      filter(!is.na(dmr_diff), !is.na(de_comb))
    
    # Sanity checks
    if (nrow(df_base) == 0) {
      message("Skipping: no data for feature=", feat, " set=", paste(set, collapse = ", "))
      next
    }
    if (length(unique(df_base$feature)) != 1 || unique(df_base$feature) != feat) {
      stop("Feature mismatch: requested ", feat,
           " but found: ", paste(unique(df_base$feature), collapse = ", "))
    }
    
    message("N points after filtering: ", nrow(df_base), " | feature=", feat,
            " | comp=", dmr_comparison_sel)
    
    # SIGNED
    df_signed <- df_base
    res_signed <- make_plot(
      df_signed,
      xlab   = "DMR difference (meth1 - meth2)",
      ylab   = paste0("log2FC × -log10(padj) (", deg_l2fc_col, ", ", deg_padj_col, ")"),
      title  = paste0("Signed: DMR diff. vs log2FC×-log10(padj) — ", dmr_comparison_sel, " — feature=", feat),
      outfile = out_signed_pdf
    )
    all_results[[length(all_results) + 1]] <- data.frame(
      dmr_condition   = dmr_condition_sel,
      dmr_comparison  = dmr_comparison_sel,
      feature         = feat,
      deg_metric      = paste0(deg_l2fc_col, " × -log10(", deg_padj_col, ")"),
      mode            = "signed",
      pearson_r       = unname(res_signed$pear$estimate),
      pearson_p       = res_signed$pear$p.value,
      spearman_rho    = unname(res_signed$spear$estimate),
      spearman_p      = res_signed$spear$p.value,
      pearson_r2         = unname(res_signed$pear$estimate)^2,
      spearman_rho2      = unname(res_signed$spear$estimate)^2,
      signed_pearson_r2  = sign(unname(res_signed$pear$estimate))  * (unname(res_signed$pear$estimate)^2),
      signed_spearman_r2 = sign(unname(res_signed$spear$estimate)) * (unname(res_signed$spear$estimate)^2),
      n_points        = nrow(df_signed),
      stringsAsFactors = FALSE
    )
    
    # ABS
    df_abs <- df_base %>% mutate(dmr_diff = abs(.data$dmr_diff), de_comb = abs(.data$de_comb))
    res_abs <- make_plot(
      df_abs,
      xlab   = "|DMR difference|",
      ylab   = paste0("|log2FC × -log10(padj)| (", deg_l2fc_col, ", ", deg_padj_col, ")"),
      title  = paste0("Absolute: |DMR diff.| vs |log2FC×-log10(padj)| — ", dmr_comparison_sel, " — feature=", feat),
      outfile = out_abs_pdf
    )
    all_results[[length(all_results) + 1]] <- data.frame(
      dmr_condition   = dmr_condition_sel,
      dmr_comparison  = dmr_comparison_sel,
      feature         = feat,
      deg_metric      = paste0(deg_l2fc_col, " × -log10(", deg_padj_col, ")"),
      mode            = "abs",
      pearson_r       = unname(res_abs$pear$estimate),
      pearson_p       = res_abs$pear$p.value,
      spearman_rho    = unname(res_abs$spear$estimate),
      spearman_p      = res_abs$spear$p.value,
      pearson_r2         = unname(res_abs$pear$estimate)^2,
      spearman_rho2      = unname(res_abs$spear$estimate)^2,
      signed_pearson_r2  = sign(unname(res_abs$pear$estimate))  * (unname(res_abs$pear$estimate)^2),
      signed_spearman_r2 = sign(unname(res_abs$spear$estimate)) * (unname(res_abs$spear$estimate)^2),
      n_points        = nrow(df_abs),
      stringsAsFactors = FALSE
    )
  }
}

# ------------------------ Save correlation summary ---------------------------
message("Total correlation rows collected: ", length(all_results))
if (length(all_results) > 0) {
  results_df <- dplyr::bind_rows(all_results) %>%
    dplyr::mutate(
      pearson_padj  = p.adjust(pearson_p,  method = "BH"),
      spearman_padj = p.adjust(spearman_p, method = "BH")
    )
  readr::write_tsv(results_df, summary_path)
  message("Wrote correlation summary: ", summary_path)
} else {
  message("No correlation rows collected; not writing summary.")
}
