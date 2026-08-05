#!/usr/bin/env Rscript
# DATE:         20250908
# AUTHOR:       MZF & ChatGPT
# VERSION:      01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   1) Generates a subtrate table nesscesary for subsequent analysis 'per_morph_gene_means.tsv'  
#   2) Runs wilcoxon tests and plots Gene counts within each morph. So this can 
#      be used to compare 2 groups of genes within each morph.
# DETAILS:
#   Using the tidy CPM table (from Step 1), compute per-morph gene-level
#   mean CPM, compare subset vs. non-subset genes via Wilcoxon rank-sum
#   within each Morph, export summary stats, and generate boxplots.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# INPUT:
#   - tidy_cpm_long.tsv (from 01_cpm_and_flag.R)
# OUTPUT (to --out_dir):
#   - per_morph_gene_means.tsv
#   - wilcoxon_per_morph.tsv        (raw p-values + BH FDR across morphs)
#   - boxplots_per_morph.pdf        (log2(CPM+1), subset vs others)
# USAGE:
#   Rscript 02_stats_and_boxplots.R \
#     --tidy ./out_cpm_flag/tidy_cpm_long.tsv \
#     --out_dir ./out_stats_plots
# NOTE:
#   Wilcoxon runs on per-gene mean CPM (by morph). Log transform is used
#   for visualization only; the test uses raw CPM means (rank-based).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

opt <- OptionParser()
opt <- add_option(opt, c("--tidy"),    type="character", help="Path to tidy_cpm_long.tsv", default = "./cpm_genes/tidy_cpm_long.tsv")
opt <- add_option(opt, c("--out_dir"), type="character", help="Output directory", default = "../../plots/GeneCounts")
opt <- add_option(opt, c("--feature"), type="character", help="Feature For labeling the plots", default = "allGenes")
args <- parse_args(opt)

feat <- args$feature
stopifnot(!is.null(args$tidy))
out_dir <- args$out_dir
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("[i] Reading tidy CPM table: ", args$tidy)
long <- read_tsv(args$tidy, show_col_types = FALSE)

req_cols <- c("ID","Sample","CPM","logCPM1","Morph","Tissue","Treatment","InSubset")
missing_cols <- setdiff(req_cols, names(long))
if (length(missing_cols) > 0) stop("tidy_cpm_long.tsv missing columns: ", paste(missing_cols, collapse=", "))

# ---------- Per-morph gene means ----------
per_gene_morph <- long %>%
  group_by(Morph, ID, InSubset) %>%
  summarize(mean_CPM = mean(CPM, na.rm = TRUE),
            mean_logCPM1 = mean(logCPM1, na.rm = TRUE),
            n_rep = n(), .groups = "drop")

write_tsv(per_gene_morph, file.path(out_dir, "per_morph_gene_means.tsv"))

# ---------- Wilcoxon tests within each Morph ----------
wilc <- per_gene_morph %>%
  group_by(Morph) %>%
  summarize(
    n_subset = sum(InSubset),
    n_others = sum(!InSubset),
    median_subset = median(mean_CPM[InSubset], na.rm = TRUE),
    median_others = median(mean_CPM[!InSubset], na.rm = TRUE),
    mean_subset   = mean(mean_CPM[InSubset], na.rm = TRUE),
    mean_others   = mean(mean_CPM[!InSubset], na.rm = TRUE),
    # Wilcoxon rank-sum on mean CPM
    p_value = {
      x <- mean_CPM[InSubset]
      y <- mean_CPM[!InSubset]
      if (length(x) > 1 && length(y) > 1) {
        suppressWarnings(wilcox.test(x, y, alternative = "two.sided")$p.value)
      } else NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

write_tsv(wilc, file.path(out_dir, "wilcoxon_per_morph.tsv"))

# Build labels of significance for annotation
wilc_labels <- wilc %>%
  mutate(label = paste0("FDR = ", signif(p_adj_BH, 3))) %>%
  select(Morph, label)


# ---------- Boxplots (log scale for viz only) ----------
# Use per-gene mean logCPM1 for visualization (less sample-to-sample noise)
plot_df <- per_gene_morph %>%
  mutate(Group = ifelse(InSubset, paste0("Genes in ",feat), "Other genes"))

# p <- ggplot(plot_df, aes(x = Group, y = mean_logCPM1)) +
#   geom_boxplot(outlier.shape = 19, alpha = 0.9) +
#   facet_wrap(~ Morph, scales = "free_y") +
#   labs(
#     title = "Mean expression per gene (log2(CPM+1)) — Subset vs. Others, by Morph",
#     x = NULL,
#     y = "Mean log2(CPM+1) per gene"
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     panel.grid.minor = element_blank(),
#     axis.text.x = element_text(angle = 20, hjust = 1)
#   )

p <- ggplot(plot_df, aes(x = Group, y = mean_logCPM1)) +
  # geom_boxplot(outlier.shape = 19, alpha = 0.9) +
  geom_point(position = position_jitter(seed = 1, width = 0.1))+
  geom_violin(alpha=0.7) +
  stat_summary(fun = mean, geom = "point", size = 2.8, shape = 23, fill = "white") +
  stat_summary(fun = mean, fun.min = mean, fun.max = mean,
               geom = "errorbar", width = 0.25, linewidth = 0.4) +
  # facet_wrap(~ Morph, scales = "free_y") +
  facet_wrap(~ Morph, nrow=1, scales = "free_y") +
  labs(
    title = paste0("Mean expression per gene (log2(CPM+1)) — overlap with ",feat," -vs- non overlap, by Morph"),
    x = NULL,
    y = "Mean log2(CPM+1) per gene"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  # ← annotate FDR on each facet
  geom_text(
    data = wilc_labels,
    aes(x = 1.5, y = max(plot_df$mean_logCPM1, na.rm = TRUE), label = label),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 3.5
  )
# print(p)


ggsave(file.path(out_dir, paste0("violin_per_morph_",feat,".pdf")), p, width = 9, height = 7)
ggsave(file.path(out_dir, paste0("violin_per_morph_",feat,".png")), p, width = 9, height = 7, dpi = 300)

message("[✓] Done. Outputs written to: ", out_dir)
