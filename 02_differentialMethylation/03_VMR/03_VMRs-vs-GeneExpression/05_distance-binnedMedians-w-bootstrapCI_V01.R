# DATE: 20260317
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Generate a metadistance plot showing how differential expression effect
#   size varies with the genomic distance between cDMRs and
#   the closest gene.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# BIOLOGICAL CONTEXT:
#   Regulatory DNA methylation effects are often strongest near genes and
#   decay with genomic distance. This script summarizes the relationship
#   between cDMR proximity to genes and the magnitude of
#   differential gene expression.
#
#   Instead of fitting a smoothing model (e.g. LOESS), which can be biased by
#   uneven data density, this approach summarizes the data in biologically
#   interpretable distance bins and computes robust statistics per bin.
#
#   This type of "metadistance" plot is commonly used in genomics to visualize
#   regulatory decay patterns relative to genomic landmarks such as TSS,
#   enhancers, or CpG islands.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# INPUT:
#   A table generated from the cDMR–DEG join step containing at least:
#
#     gene_id
#     comparison
#     abs_l2fc        (absolute log2 fold change)
#     abs_dist_bp     (absolute distance between scDMR and closest gene)
#
#   Typically:
#
#     cDMRs_diff<DIFF_THRESH>.overlap_to_mRNA.joined_DEGs.tsv
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# METHOD OVERVIEW:
#
#   1) Filter dataset to scDMRs within a defined maximum distance from genes.
#
#   2) Assign each observation to a biologically meaningful distance bin:
#
#        0 bp
#        1–500 bp
#        501–1 kb
#        1–2 kb
#        2–5 kb
#        5–10 kb
#        10–20 kb
#        20–50 kb
#
#   3) For each comparison and distance bin:
#
#        - compute the median |log2FC|
#        - estimate uncertainty using bootstrap resampling
#
#   4) Plot:
#
#        X-axis  : genomic distance bin
#        Y-axis  : median |log2FC|
#
#        Ribbon  : 95% bootstrap confidence interval
#        Points  : median per bin
#        Point size proportional to number of scDMRs in the bin
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# OUTPUT:
#
#   Summary statistics table:
#
#     *.absDmrDistance_vs_absL2FC.metadistance.tsv
#
#   Metadistance figure:
#
#     *.absDmrDistance_vs_absL2FC.metadistance.pdf
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# ADVANTAGES OF THIS APPROACH:
#
#   • Robust to uneven sampling density across genomic distance
#   • Avoids LOESS overfitting artifacts
#   • Directly interpretable biological bins
#   • Provides uncertainty estimates via bootstrap confidence intervals
#
#   This methodology is widely used in analyses of:
#
#     - DNA methylation vs TSS distance
#     - enhancer–gene regulatory proximity
#     - ChIP-seq signal decay
#     - chromatin interaction distance effects
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%


library(dplyr)
library(ggplot2)

#. NOTE!!!
# Run right after script 04


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set Bin sizes ----------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Bin_size_bp <- 100
XMAX <- 100000
N_BOOT <- 1000


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare distance bins ----------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_df_boot <- plot_df3 %>%
  filter(!is.na(abs_dist_bp), abs_dist_bp <= XMAX) %>%
  mutate(
    dist_bin = case_when(
      abs_dist_bp == 0 ~ 0,
      abs_dist_bp > 0  ~ ((abs_dist_bp - 1) %/% Bin_size_bp + 1) * Bin_size_bp
    )
  )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bootstrap function --------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boot_median <- function(x, nboot = N_BOOT) {
  if(length(x) < 2) return(c(med = median(x), lo = NA, hi = NA))
  
  boots <- replicate(nboot, median(sample(x, replace = TRUE)))
  
  c(
    med = median(x),
    lo  = quantile(boots, 0.025),
    hi  = quantile(boots, 0.975)
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute bin statistics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boot_df <- plot_df_boot %>%
  group_by(comparison, dist_bin) %>%
  summarise(
    stats = list(boot_median(abs_l2fc)),
    n = n(),
    .groups = "drop"
  ) %>%
  tidyr::unnest_wider(stats)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot =======================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p_boot <- ggplot(boot_df, aes(x = dist_bin, y = med)) +
  
  geom_ribbon(
    aes(ymin = lo, ymax = hi),
    # aes(ymin = lo, ymax = hi),
    alpha = 0.2,
    fill = "blue"
  ) +
  
  geom_line(
    color = "steelblue",
    linewidth = 1
  ) +
  
  geom_point(
    size = 1
  ) +
  
  facet_wrap(~comparison, scales = "free_y") +
  
  labs(
    title = paste0("Effect size vs distance to closest gene, bootstrap = ",N_BOOT),
    x = paste0("Distance to closest gene (bp)- intervals of ",Bin_size_bp),
    y = "|log2 fold change| (median per bin)"
  ) +

  # scale_x_log10()+
  theme_classic(base_size = 12)

print(p_boot)

ggsave(
  paste0(
    plots_folder,
    "/closest_mRNA.cDMRs_diff",
    DIFF_THRESH,
    ".absDmrDistance_vs_absL2FC.bootstrapTrend.pdf"
  ),
  p_boot,
  width = 10,
  height = 7
)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Metadistance plot: =========================================
# median |log2FC| across distance bins to closest mRNA
# with bootstrap 95% CI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(dplyr)
library(tidyr)
library(ggplot2)

title_base <- "Metadistance plot: cDMR distance to closest gene vs differential expression"
y_lab <- "|log2 fold change|"
N_BOOT <- 1000

# distance breaks (customizable)
dist_breaks <- c(0, 1, 500, 1000, 2000, 5000, 10000, 20000, 50000)
dist_labels <- c("0", "1-500", "501-1k", "1-2k", "2-5k", "5-10k", "10-20k", "20-50k")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_meta <- plot_df3 %>%
  filter(!is.na(abs_dist_bp), !is.na(abs_l2fc), abs_dist_bp <= 50000) %>%
  mutate(
    dist_bin = cut(
      abs_dist_bp,
      breaks = dist_breaks,
      labels = dist_labels,
      include.lowest = TRUE,
      right = TRUE
    ),
    dist_bin = factor(dist_bin, levels = dist_labels)
  ) %>%
  filter(!is.na(dist_bin))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bootstrap helper
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
boot_median <- function(x, nboot = N_BOOT) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(c(med = NA, lo = NA, hi = NA))
  if (length(x) == 1) return(c(med = x, lo = x, hi = x))
  
  boots <- replicate(nboot, median(sample(x, replace = TRUE)))
  c(
    med = median(x),
    lo  = as.numeric(quantile(boots, 0.025)),
    hi  = as.numeric(quantile(boots, 0.975))
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summarize per comparison x distance bin
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meta_df <- plot_meta %>%
  group_by(comparison, dist_bin) %>%
  summarise(
    n = n(),
    stats = list(boot_median(abs_l2fc)),
    .groups = "drop"
  ) %>%
  unnest_wider(stats) %>%
  mutate(
    bin_index = as.numeric(dist_bin)
  )

# optional export
fwrite(
  as.data.table(meta_df),
  paste0(
    plots_folder,
    "/closest_mRNA.cDMRs_diff",
    DIFF_THRESH,
    ".absDmrDistance_vs_absL2FC.metadistance.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p_meta <- ggplot(meta_df, aes(x = bin_index, y = med, group = 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey70", alpha = 0.35) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(size = n), alpha = 0.9) +
  facet_wrap(~ comparison, scales = "free_y") +
  scale_x_continuous(
    breaks = seq_along(dist_labels),
    labels = dist_labels
  ) +
  scale_size_continuous(name = "n scDMRs") +
  labs(
    title = title_base,
    x = "Absolute distance to closest gene (bp)",
    y = paste0(y_lab, " (median ± 95% bootstrap CI)")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_meta)

ggsave(
  paste0(
    plots_folder,
    "/closest_mRNA.cDMRs_diff",
    DIFF_THRESH,
    ".absDmrDistance_vs_absL2FC.metadistance.pdf"
  ),
  plot = p_meta,
  width = 10,
  height = 7
)
