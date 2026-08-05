#!/usr/bin/env Rscript
# DATE: 20260317
# AUTHOR: MZF
# SCRIPT: 05_scatterPlots_meth-vs-expression_V01.R
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Generate feature-specific methylation vs expression plots and panel-level
#   statistics from the joined long table produced by:
#       04_Join_GeneCounts_V04.R
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# DESCRIPTION:
#   This script reads the reusable joined long table:
#       gene x feature x morph
#
#   It then:
#     A) filters the joined table to one chosen feature
#     B) collapses to one row per gene x morph for plotting
#     C) writes feature-specific tables
#     D) computes panel-level statistics
#     E) generates plots of methylation vs expression
#
# IMPORTANT DESIGN CHOICE:
#   The feature filter is applied near the top of the script so the entire
#   plotting workflow is explicitly tied to one biological feature at a time.
#
# FEATURE EXAMPLES:
#   - "gene"
#   - "geneBody"
#   - "d2000"
#   - "intron"
#
# INPUT FILE:
#   - ./scDMR_geneSubsets_cpm/methPerFeature_AND_cpm.long.tsv
#
# OUTPUTS:
#   1) Feature-specific joined table
#   2) Collapsed plotting table (one row per gene x morph)
#   3) Linear / correlation panel statistics table
#   4) Quadratic panel statistics table
#   5) GAM panel statistics table
#   6) Linear scatter plot with panel annotations
#   7) Quadratic plot with binned boxplots + panel annotations
#
# COLLAPSING RULE:
#   If a gene has multiple rows for the chosen feature within a morph:
#     - use weighted mean methylation when coverage is available
#     - otherwise use simple mean methylation
#
# STATISTICAL NOTES:
#   - Linear plot annotations use:
#       linear model (R²) + Spearman correlation
#   - Quadratic plot annotations use:
#       y ~ x + x^2
#     and report whether the pattern is consistent with an inverted-U
#   - GAM statistics are also written as a flexible non-linear diagnostic
#
# EXPECTED RUN ORDER:
#   - Run 04_Join_GeneCounts_V04.R first
#   - Then run this plotting script
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
 
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(mgcv)
})

# Clear grid# Clear environment and graphics
rm(list = ls())
graphics.off()

# ------------------------------------------------------------
# User-defined feature filter
# Move feature filtering to the top of the plotting workflow
# ------------------------------------------------------------
# feature_to_plot <- "u2000"
# feature_to_plot <- "gene"
# feature_to_plot <- "geneBody"
feature_to_plot <- "d2000"

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
# joined_file   <- "./scDMR_geneSubsets_cpm/methPerFeature_AND_cpm.long.tsv"
joined_file   <- "./methPerFeature_AND_cpm.long.tsv"
plots_folder  <- paste0("../../plots/GM-vs-GE",feature_to_plot)
tables_folder <- "./gene_meth-vs-cpm"

dir.create(plots_folder,  showWarnings = FALSE, recursive = TRUE)
dir.create(tables_folder, showWarnings = FALSE, recursive = TRUE)



# ------------------------------------------------------------
# Plot settings
# ------------------------------------------------------------
morph_order <- c("female", "immaculata", "parae", "yellow")

morph_colors <- c(
  "female"      = "#149954",
  "immaculata"  = "#888888",
  "parae"       = "#9e9ac8",
  "yellow"      = "#fd8d3c"
)

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------
safe_feature <- gsub("[^A-Za-z0-9_\\-]", "_", feature_to_plot)

out_feature_table <- file.path(
  tables_folder,
  paste0("methPerFeature_AND_cpm.feature_", safe_feature, ".tsv")
)

out_plot_table <- file.path(
  tables_folder,
  paste0("methPerFeature_AND_cpm.feature_", safe_feature, ".collapsed.tsv")
)

out_linear_stats <- file.path(
  plots_folder,
  paste0("scatter_meth_vs_mean_logCPM1_feature_", safe_feature, ".panel_stats_linear_spearman.tsv")
)

out_quadratic_stats <- file.path(
  plots_folder,
  paste0("scatter_meth_vs_mean_logCPM1_feature_", safe_feature, ".panel_stats_quadratic.tsv")
)

out_gam_stats <- file.path(
  plots_folder,
  paste0("scatter_meth_vs_mean_logCPM1_feature_", safe_feature, ".panel_stats_gam.tsv")
)

out_linear_plot <- file.path(
  plots_folder,
  paste0("scatter_meth_vs_mean_logCPM1_feature_", safe_feature, ".spearman.pdf")
)

out_quadratic_plot <- file.path(
  plots_folder,
  paste0("boxplots_meth_vs_mean_logCPM1_feature_", safe_feature, ".quadratic.pdf")
)

# ------------------------------------------------------------
# Read joined table
# ------------------------------------------------------------
joined_dt <- fread(joined_file, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))

cat("Read joined rows: ", nrow(joined_dt), "\n", sep = "")
cat("Feature selected for plotting: ", feature_to_plot, "\n", sep = "")

# ------------------------------------------------------------
# Required columns check
# ------------------------------------------------------------
required_cols <- c(
  "gene_id", "Morph_join", "Morph_meth", "Morph_cpm",
  "feature", "chromosome_type",
  "cov_morph", "meth_morph",
  "mean_CPM", "mean_logCPM1", "InSubset", "n_rep"
)

missing_required <- setdiff(required_cols, names(joined_dt))
if (length(missing_required) > 0) {
  stop(
    "Missing expected columns in joined table: ",
    paste(missing_required, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Feature filter
# ------------------------------------------------------------
feature_dt <- joined_dt[
  feature == feature_to_plot &
    !is.na(meth_morph) &
    !is.na(mean_logCPM1)
]

if (nrow(feature_dt) == 0) {
  stop("No rows remained after filtering for feature == '", feature_to_plot, "'.")
}

fwrite(
  feature_dt,
  out_feature_table,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote feature-filtered joined table: ", out_feature_table, "\n", sep = "")
cat("Rows after feature filter: ", nrow(feature_dt), "\n", sep = "")

# ------------------------------------------------------------
# Collapse to one row per gene x morph for plotting
# ------------------------------------------------------------
gene_plot_dt <- feature_dt[, {
  valid_w <- !is.na(meth_morph) & !is.na(cov_morph) & cov_morph > 0
  valid_x <- !is.na(meth_morph)
  
  meth_weighted <- if (any(valid_w)) {
    weighted.mean(meth_morph[valid_w], cov_morph[valid_w])
  } else if (any(valid_x)) {
    mean(meth_morph[valid_x])
  } else {
    NA_real_
  }
  
  list(
    n_feature_rows  = .N,
    cov_morph_sum   = if (all(is.na(cov_morph))) NA_real_ else sum(cov_morph, na.rm = TRUE),
    meth_morph_plot = meth_weighted,
    mean_CPM        = mean_CPM[1],
    mean_logCPM1    = mean_logCPM1[1],
    InSubset        = InSubset[1],
    n_rep           = n_rep[1],
    chromosome_type = chromosome_type[1]
  )
}, by = .(gene_id, Morph_join, Morph_meth, Morph_cpm)]

fwrite(
  gene_plot_dt,
  out_plot_table,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote collapsed plotting table: ", out_plot_table, "\n", sep = "")
cat("Rows in collapsed plotting table: ", nrow(gene_plot_dt), "\n", sep = "")

# ------------------------------------------------------------
# Final plotting table
# ------------------------------------------------------------
plot_dt <- gene_plot_dt[
  !is.na(meth_morph_plot) &
    !is.na(mean_logCPM1)
]

if (nrow(plot_dt) == 0) {
  stop("No rows remained in plot_dt after removing missing methylation/expression values.")
}

plot_dt[, Morph_join := factor(Morph_join, levels = morph_order)]

# ------------------------------------------------------------
# Linear model + Spearman correlation stats for panel annotation
# ------------------------------------------------------------
linear_stats_dt <- plot_dt[, {
  ok <- is.finite(meth_morph_plot) & is.finite(mean_logCPM1)
  x <- meth_morph_plot[ok]
  y <- mean_logCPM1[ok]
  n <- length(x)
  
  x_rng <- range(x, na.rm = TRUE)
  y_rng <- range(y, na.rm = TRUE)
  x_span <- diff(x_rng)
  y_span <- diff(y_rng)
  
  x_pos <- if (is.finite(x_span) && x_span > 0) x_rng[1] + 0.03 * x_span else x_rng[1]
  y_pos <- if (is.finite(y_span) && y_span > 0) y_rng[2] - 0.03 * y_span else y_rng[2]
  
  if (n >= 3 && length(unique(x)) >= 2) {
    fit <- lm(y ~ x)
    fit_sum <- summary(fit)
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    
    intercept  <- unname(coef(fit)[1])
    slope      <- unname(coef(fit)[2])
    r2         <- fit_sum$r.squared
    adj_r2     <- fit_sum$adj.r.squared
    spearman_r <- unname(ct$estimate)
    p_spearman <- ct$p.value
    
    label <- paste0(
      "n = ", n, "\n",
      "y = ", sprintf("%.3f", intercept),
      ifelse(slope >= 0, " + ", " - "),
      sprintf("%.3f", abs(slope)), "x\n",
      "R² = ", sprintf("%.3f", r2), "\n",
      "rho = ", sprintf("%.3f", spearman_r),
      ", p = ", format.pval(p_spearman, digits = 2, eps = 1e-3)
    )
  } else {
    intercept  <- NA_real_
    slope      <- NA_real_
    r2         <- NA_real_
    adj_r2     <- NA_real_
    spearman_r <- NA_real_
    p_spearman <- NA_real_
    
    label <- paste0("n = ", n, "\n", "Insufficient variation")
  }
  
  .(
    n = n,
    intercept = intercept,
    slope = slope,
    r2 = r2,
    adj_r2 = adj_r2,
    spearman_r = spearman_r,
    p_spearman = p_spearman,
    x_pos = x_pos,
    y_pos = y_pos,
    label = label
  )
}, by = .(chromosome_type, Morph_join)]

fwrite(
  linear_stats_dt,
  out_linear_stats,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote linear/spearman panel statistics table: ", out_linear_stats, "\n", sep = "")

# ------------------------------------------------------------
# Quadratic model stats for panel annotation
# Model: y ~ x + x^2
# ------------------------------------------------------------
quadratic_stats_dt <- plot_dt[, {
  ok <- is.finite(meth_morph_plot) & is.finite(mean_logCPM1)
  x <- meth_morph_plot[ok]
  y <- mean_logCPM1[ok]
  n <- length(x)
  
  x_rng <- range(x, na.rm = TRUE)
  y_rng <- range(y, na.rm = TRUE)
  x_span <- diff(x_rng)
  y_span <- diff(y_rng)
  
  x_pos <- if (is.finite(x_span) && x_span > 0) x_rng[1] + 0.03 * x_span else x_rng[1]
  y_pos <- if (is.finite(y_span) && y_span > 0) y_rng[2] - 0.03 * y_span else y_rng[2]
  
  if (n >= 5 && length(unique(x)) >= 3) {
    fit_lin  <- lm(y ~ x)
    fit_quad <- lm(y ~ x + I(x^2))
    
    s_quad <- summary(fit_quad)
    coefs  <- coef(s_quad)
    
    b0 <- unname(coef(fit_quad)[1])
    b1 <- unname(coef(fit_quad)[2])
    b2 <- unname(coef(fit_quad)[3])
    
    p_x  <- coefs["x", "Pr(>|t|)"]
    p_x2 <- coefs["I(x^2)", "Pr(>|t|)"]
    
    r2     <- s_quad$r.squared
    adj_r2 <- s_quad$adj.r.squared
    
    cmp <- anova(fit_lin, fit_quad)
    p_quad_vs_linear <- cmp$`Pr(>F)`[2]
    
    vertex_x <- if (!is.na(b2) && b2 != 0) -b1 / (2 * b2) else NA_real_
    
    shape <- if (!is.na(b2) && !is.na(p_x2) && p_x2 < 0.05) {
      if (b2 < 0) "inverted-U" else if (b2 > 0) "U-shaped" else "no clear quadratic"
    } else {
      "no clear quadratic"
    }
    
    label <- paste0(
      "n = ", n, "\n",
      "beta_x = ", sprintf("%.3f", b1), "\n",
      "beta_x2 = ", sprintf("%.3f", b2), "\n",
      "p_x2 = ", format.pval(p_x2, digits = 2, eps = 1e-3), "\n",
      "p_quad_vs_linear = ", format.pval(p_quad_vs_linear, digits = 2, eps = 1e-3), "\n",
      "vertex_x = ", ifelse(is.na(vertex_x), "NA", sprintf("%.3f", vertex_x)), "\n",
      "R² = ", sprintf("%.3f", r2), "\n",
      "shape = ", shape
    )
  } else {
    b0 <- b1 <- b2 <- p_x <- p_x2 <- r2 <- adj_r2 <- p_quad_vs_linear <- vertex_x <- NA_real_
    shape <- "insufficient data"
    label <- paste0("n = ", n, "\n", "insufficient data")
  }
  
  .(
    n = n,
    intercept = b0,
    beta_x = b1,
    beta_x2 = b2,
    p_x = p_x,
    p_x2 = p_x2,
    r2 = r2,
    adj_r2 = adj_r2,
    p_quad_vs_linear = p_quad_vs_linear,
    vertex_x = vertex_x,
    shape = shape,
    x_pos = x_pos,
    y_pos = y_pos,
    label = label
  )
}, by = .(chromosome_type, Morph_join)]

fwrite(
  quadratic_stats_dt,
  out_quadratic_stats,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote quadratic panel statistics table: ", out_quadratic_stats, "\n", sep = "")

# ------------------------------------------------------------
# GAM stats
# ------------------------------------------------------------
gam_stats_dt <- plot_dt[, {
  ok <- is.finite(meth_morph_plot) & is.finite(mean_logCPM1)
  x <- meth_morph_plot[ok]
  y <- mean_logCPM1[ok]
  n <- length(x)
  
  if (n >= 10 && length(unique(x)) >= 5) {
    d <- data.frame(x = x, y = y)
    
    fit_lin <- lm(y ~ x, data = d)
    fit_gam <- gam(y ~ s(x, k = 5), data = d, method = "REML")
    
    s_gam <- summary(fit_gam)
    
    edf      <- s_gam$s.table[1, "edf"]
    fval     <- s_gam$s.table[1, "F"]
    p_smooth <- s_gam$s.table[1, "p-value"]
    dev_expl <- s_gam$dev.expl
    aic_lin  <- AIC(fit_lin)
    aic_gam  <- AIC(fit_gam)
  } else {
    edf <- fval <- p_smooth <- dev_expl <- aic_lin <- aic_gam <- NA_real_
  }
  
  .(
    n = n,
    edf = edf,
    F_smooth = fval,
    p_smooth = p_smooth,
    deviance_explained = dev_expl,
    AIC_linear = aic_lin,
    AIC_gam = aic_gam
  )
}, by = .(chromosome_type, Morph_join)]

fwrite(
  gam_stats_dt,
  out_gam_stats,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote GAM panel statistics table: ", out_gam_stats, "\n", sep = "")

# ------------------------------------------------------------
# Linear scatter plot with linear fit + Spearman annotation
# ------------------------------------------------------------
p_linear <- ggplot(
  plot_dt,
  aes(x = meth_morph_plot, y = mean_logCPM1, color = Morph_join)
) +
  geom_point(alpha = 0.10, size = 1) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_text(
    data = linear_stats_dt,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  scale_color_manual(values = morph_colors, drop = FALSE) +
  facet_grid(
    rows = vars(chromosome_type),
    cols = vars(Morph_join),
    scales = "free"
  ) +
  labs(
    title = "Gene methylation vs expression",
    subtitle = paste0("feature == '", feature_to_plot, "'; one point per gene x morph"),
    x = "Morph-level methylation",
    y = "mean_logCPM1",
    color = "Morph"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95")
  )

print(p_linear)

ggsave(out_linear_plot, p_linear, width = 10, height = 8)

cat("Wrote linear scatter plot: ", out_linear_plot, "\n", sep = "")

# ------------------------------------------------------------
# Quadratic plot with binned boxplots + quadratic annotation
# ------------------------------------------------------------
p_quadratic <- ggplot(
  plot_dt,
  aes(x = meth_morph_plot, y = mean_logCPM1, color = Morph_join)
) +
  # geom_point(alpha = 0.05, size = 0.8) +
  geom_boxplot(
    aes(group = cut_width(meth_morph_plot, width = 0.05)),
    outlier.alpha = 0.08,
    width = 0.04,
    linewidth = 0.35
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x + I(x^2),
    se = TRUE,
    linewidth = 0.9,
    color = "black"
  ) +
  geom_text(
    data = quadratic_stats_dt,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.8,
    lineheight = 0.95
  ) +
  scale_color_manual(values = morph_colors, drop = FALSE) +
  facet_grid(
    rows = vars(chromosome_type),
    cols = vars(Morph_join),
    scales = "free"
  ) +
  labs(
    title = "Gene methylation vs expression",
    subtitle = paste0(
      "feature == '", feature_to_plot,
      "'; one point per gene x morph; quadratic fit per panel"
    ),
    x = "% methylation",
    y = "gene expression (mean_logCPM1)",
    color = "Morph"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95")
  )

print(p_quadratic)

ggsave(out_quadratic_plot, p_quadratic, width = 12, height = 8)

cat("Wrote quadratic plot: ", out_quadratic_plot, "\n", sep = "")
cat("Done.\n")

