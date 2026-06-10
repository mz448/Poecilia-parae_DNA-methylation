#!/usr/bin/env Rscript
# DATE: 20260316
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Plot cDMR score against absolute differential expression effect size
#   using the joined cDMR–DEG table.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# INPUT:
#   cDMRs_diff<DIFF_THRESH>.overlap_to_mRNA.joined_DEGs.tsv
#
# X:
#   dmr_score
#
# Y:
#   abs(log2FoldChange)
#
# DETAILS:
#   - Select only the requested DEG comparisons
#   - Reshape *_l2fc columns into long format
#   - Compute abs_l2fc = abs(l2fc)
#   - Keep only rows with non-NA dmr_score and abs_l2fc
#   - Plot one scatter plot faceted by comparison
#
# IMPORTANT:
#   - Because the join keeps all DEG rows, genes without an associated cDMR
#     will have NA in dmr_score and will be excluded from plotting.
#   - If one gene is linked to multiple cDMRs, it can appear multiple times.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# ~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Inputs
# ~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
DIFF_THRESH <- "0"
plots_folder <- paste0("../plots_", DIFF_THRESH)


# INTRAGENIC (Overlap with mRNA)
infile <- paste0(plots_folder,"/cDMRs_diff",DIFF_THRESH,".overlap_to_mRNA.joined_DEGs.tsv")
out_plot <- paste0(plots_folder,"/overlap_to_mRNA.cDMRs_diff",DIFF_THRESH,".dmrScore_vs_absL2FC.scatter.pdf")
out_plot_png <- paste0(plots_folder,"/overlap_to_mRNA.cDMRs_diff",DIFF_THRESH,".dmrScore_vs_absL2FC.scatter.png")
out_long <- paste0(plots_folder,"/overlap_to_mRNA.cDMRs_diff",DIFF_THRESH,".dmrScore_vs_absL2FC.plotTable.tsv")


# closest_mRNA
# infile <- paste0(plots_folder,"/cDMRs_diff",DIFF_THRESH,".closest_mRNA.joined_DEGs.tsv")
# out_plot <- paste0(plots_folder,"/closest_mRNA.cDMRs_diff",DIFF_THRESH,".dmrScore_vs_absL2FC.scatter.pdf")
# out_plot_png <- paste0(plots_folder,"/closest_mRNA.cDMRs_diff",DIFF_THRESH,".dmrScore_vs_absL2FC.scatter.png")
# out_long <- paste0(plots_folder,"/closest_mRNA.cDMRs_diff",DIFF_THRESH,".dmrScore_vs_absL2FC.plotTable.tsv")




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read data ----------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merged2 <- fread(infile, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic checks ----------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
required_cols <- c("gene_id", "gene_name", "dmr_score")
missing_cols <- setdiff(required_cols, names(merged2))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Select desired comparisons + *_l2fc columns ----------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("[1/4] Selecting requested *_l2fc columns")

wanted_prefix <- c(
  "i_vs_y_01",
  "p_vs_i_01",
  "y_vs_p_01",
  "f_vs_i",
  "f_vs_p",
  "f_vs_y"
)

l2fc_cols <- names(merged2)[
  str_detect(names(merged2), paste0("(", paste(wanted_prefix, collapse = "|"), ")")) &
    str_detect(names(merged2), "_l2fc$")
]

if (length(l2fc_cols) == 0) {
  stop("No *_l2fc columns found matching requested comparisons. Check column names/prefixes.")
}
message("  Found ", length(l2fc_cols), " l2fc columns.")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INTRAGENIC #################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Build long table for plotting ----------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("[2/4] Building long plotting table")


# Intragenic Dataframe ---------------------------------------
plot_df <- merged2 %>%
  filter(
    abs_dist_bp==0
  ) %>%
  select(
    gene_id,
    gene_name,
    dmr_chr,
    dmr_start,
    dmr_end,
    cons_id,
    cluster,
    chromosome_type,
    dist_signed,
    abs_dist_bp,
    dmr_score, # Adding the DMR score == Max meth difference observed across groups (means)
    all_of(l2fc_cols)
  ) %>%
  pivot_longer(
    cols = all_of(l2fc_cols),
    names_to = "comparison",
    values_to = "l2fc"
  ) %>%
  mutate(
    l2fc = suppressWarnings(as.numeric(l2fc)),
    abs_l2fc = abs(l2fc),
    comparison = str_remove(comparison, "_l2fc$")
  ) %>%
  filter(!is.na(dmr_score), !is.na(abs_l2fc))

if (nrow(plot_df) == 0) {
  stop("No rows left after filtering for non-NA dmr_score and abs_l2fc values.")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write plotting table ----------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("[3/4] Writing plotting table")
fwrite(plot_df, out_long, sep = "\t", quote = FALSE, na = "NA")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Scatter: cDMR MethDiff vs differential expression------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("[4/4] Plotting intragenic")

title_base <- "cDMR MethDiff vs differential expression"
x_lab <- "cDMR score"
y_lab <- "|log2 fold change|"

p <- ggplot(plot_df, aes(x = dmr_score, y = abs_l2fc)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
  facet_wrap(~ comparison, scales = "free_y") +
  labs(
    title = title_base,
    x = x_lab,
    y = y_lab
  ) +
  theme_classic(base_size = 11)
print (p)
ggsave(out_plot, plot = p, width = 10, height = 7)
# ggsave(out_plot_png, plot = p, width = 10, height = 7, dpi = 300)

cat("Wrote: ", out_long, "\n", sep = "")
cat("Wrote: ", out_plot, "\n", sep = "")
cat("Wrote: ", out_plot_png, "\n", sep = "")
cat("Rows in plotting table: ", nrow(plot_df), "\n", sep = "")




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# INTERGENIC #################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# intergenic
infile <- paste0(plots_folder,"/cDMRs_diff",DIFF_THRESH,".intergenic.joined_DEGs.tsv")
out_plot <- paste0(plots_folder,"/intergenic.cDMRs_diff",DIFF_THRESH,".dmrScore_vs_absL2FC.scatter.pdf")
out_plot_png <- paste0(plots_folder,"/intergenic.cDMRs_diff",DIFF_THRESH,".dmrScore_vs_absL2FC.scatter.png")
out_long <- paste0(plots_folder,"/intergenic.cDMRs_diff",DIFF_THRESH,".dmrScore_vs_absL2FC.plotTable.tsv")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merged2 <- fread(infile, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic checks
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
required_cols <- c("gene_id", "gene_name", "dmr_score")
missing_cols <- setdiff(required_cols, names(merged2))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Select desired comparisons + *_l2fc columns
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("[1/4] Selecting requested *_l2fc columns")

wanted_prefix <- c(
  "i_vs_y_01",
  "p_vs_i_01",
  "y_vs_p_01",
  "f_vs_i",
  "f_vs_p",
  "f_vs_y"
)

l2fc_cols <- names(merged2)[
  str_detect(names(merged2), paste0("(", paste(wanted_prefix, collapse = "|"), ")")) &
    str_detect(names(merged2), "_l2fc$")
]

if (length(l2fc_cols) == 0) {
  stop("No *_l2fc columns found matching requested comparisons. Check column names/prefixes.")
}
message("  Found ", length(l2fc_cols), " l2fc columns.")



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot Scatter vs distance (intergenic cDMRs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Intergenic Dataframe
plot_df2 <- merged2 %>%
  filter(
    abs_dist_bp!=0
  ) %>%
  select(
    gene_id,
    gene_name,
    dmr_chr,
    dmr_start,
    dmr_end,
    cons_id,
    cluster,
    chromosome_type,
    dmr_score,
    dist_signed, # Distance
    abs_dist_bp, # ABS distance
    all_of(l2fc_cols)
  ) %>%
  pivot_longer(
    cols = all_of(l2fc_cols),
    names_to = "comparison",
    values_to = "l2fc"
  ) %>%
  mutate(
    l2fc = suppressWarnings(as.numeric(l2fc)),
    abs_l2fc = abs(l2fc),
    comparison = str_remove(comparison, "_l2fc$")
  ) %>%
  filter(!is.na(dist_signed), !is.na(abs_l2fc))

# plot
title_base <- "cDMR distance (to a Gene) vs differential expression"
x_lab <- "cDMR distance to gene"
y_lab <- "|log2 fold change|"


p2 <- ggplot(plot_df2, aes(x = dist_signed, y = abs_l2fc, color = dmr_score)) +
  geom_point(alpha = 1, size = .5) +
  facet_wrap(~ comparison, scales = "free_y") +
  labs(
    title = title_base,
    x = x_lab,
    y = y_lab,
    color = "cDMR score"
  ) +
  scale_color_gradient(
    low = "grey",
    high = "black"
  ) +
  theme_classic(base_size = 11)

print(p2)

ggsave(paste0(plots_folder,"/intergenic.cDMRs_diff",DIFF_THRESH,".dmrDistance_vs_absL2FC.scatter.pdf"), plot = p2, width = 10, height = 7)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot boxplots per bin of distance (intergenic cDMRs)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Plot details
title_base <- "cDMR distance (to a Gene) vs differential expression"
x_lab <- "cDMR distance to gene"
y_lab <- "|log2 fold change|"

Bin_size_bp <-  5000
p3 <- ggplot(plot_df2, aes(x = dist_signed, y = abs_l2fc)) +
  geom_boxplot(aes(group = cut_width(dist_signed, Bin_size_bp)), outlier.alpha = 0.1)+
  facet_wrap(~ comparison, scales = "free_y") +
  labs(
    title = title_base,
    x = x_lab,
    y = y_lab
  ) +
  scale_color_gradient(
    low = "grey",
    high = "black"
  ) +
  coord_cartesian(xlim = c(-50000,50000))+
  theme_classic(base_size = 11)

print(p3)

ggsave(paste0(plots_folder,"/intergenic.cDMRs_diff",DIFF_THRESH,".dmrDistance_vs_absL2FC.boxplots.pdf"), 
       plot = p3, width = 10, height = 7)

# Load all the distance variables 
# All variables Dataframe
plot_df3 <- merged2 %>%
  select(
    gene_id,
    gene_name,
    dmr_chr,
    dmr_start,
    dmr_end,
    cons_id,
    cluster,
    chromosome_type,
    dmr_score,
    dist_signed,
    abs_dist_bp,
    all_of(l2fc_cols)
  ) %>%
  pivot_longer(
    cols = all_of(l2fc_cols),
    names_to = "comparison",
    values_to = "l2fc"
  ) %>%
  mutate(
    l2fc = suppressWarnings(as.numeric(l2fc)),
    abs_l2fc = abs(l2fc),
    comparison = str_remove(comparison, "_l2fc$")
  ) %>%
  filter(!is.na(dist_signed), !is.na(abs_l2fc))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# == Plot boxplots per bin of distance (intergenic cDMRs) ====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot details
title_base <- "cDMR distance (to a Gene) vs differential expression"
Bin_size_bp <- 5000
x_lab <- paste0("cDMR distance to gene (0 as single bin; ±", Bin_size_bp, "-bp bins)")
y_lab <- "|log2 fold change|"

plot_df4 <- plot_df3 %>%
  filter(dist_signed >= -50000, dist_signed <= 50000) %>%
  mutate(
    dist_bin_start = case_when(
      dist_signed == 0 ~ 0,
      dist_signed > 0  ~ ((dist_signed - 1) %/% Bin_size_bp) * Bin_size_bp + 1,
      dist_signed < 0  ~ -((((abs(dist_signed) - 1) %/% Bin_size_bp) + 1) * Bin_size_bp)
    ),
    dist_bin_end = case_when(
      dist_signed == 0 ~ 0,
      dist_signed > 0  ~ dist_bin_start + Bin_size_bp - 1,
      dist_signed < 0  ~ dist_bin_start + Bin_size_bp - 1
    ),
    dist_bin_tick = case_when(
      dist_signed == 0 ~ "0",
      dist_signed > 0  ~ as.character(dist_bin_end),
      dist_signed < 0  ~ as.character(dist_bin_start)
    )
  )

# order bins numerically
bin_levels <- plot_df4 %>%
  distinct(dist_bin_start, dist_bin_tick) %>%
  arrange(dist_bin_start) %>%
  pull(dist_bin_tick)

plot_df4 <- plot_df4 %>%
  mutate(dist_bin_tick = factor(dist_bin_tick, levels = bin_levels))

p4 <- ggplot(plot_df4, aes(x = dist_bin_tick, y = abs_l2fc)) +
  geom_boxplot(outlier.alpha = 0.1, width = 0.7) +
  facet_wrap(~ comparison, scales = "free_y") +
  labs(
    title = title_base,
    x = x_lab,
    y = y_lab
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p4)

ggsave(
  paste0(
    plots_folder,
    "/closest_mRNA.cDMRs_diff",
    DIFF_THRESH,
    ".dmrDistance_vs_absL2FC.boxplots.pdf"
  ),
  plot = p4,
  width = 10,
  height = 7
)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot boxplots per bin of ABSOLUTE distance (intergenic cDMRs) =====
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Plot details
title_base <- "cDMR absolute distance to closest gene vs differential expression"
Bin_size_bp <- 1000
x_lab <- paste0("Absolute cDMR distance to gene (bp; 0 as single bin, then ", Bin_size_bp, "-bp bins)")
y_lab <- "|log2 fold change|"

plot_df4 <- plot_df3 %>%
  filter(!is.na(abs_dist_bp), abs_dist_bp <= 20000) %>%
  mutate(
    dist_bin_start = case_when(
      abs_dist_bp == 0 ~ 0,
      abs_dist_bp > 0  ~ ((abs_dist_bp - 1) %/% Bin_size_bp) * Bin_size_bp + 1
    ),
    dist_bin_end = case_when(
      abs_dist_bp == 0 ~ 0,
      abs_dist_bp > 0  ~ dist_bin_start + Bin_size_bp - 1
    ),
    dist_bin_tick = case_when(
      abs_dist_bp == 0 ~ "0",
      abs_dist_bp > 0  ~ as.character(dist_bin_end)
    )
  )

# order bins numerically
bin_levels <- plot_df4 %>%
  distinct(dist_bin_start, dist_bin_tick) %>%
  arrange(dist_bin_start) %>%
  pull(dist_bin_tick)

plot_df4 <- plot_df4 %>%
  mutate(dist_bin_tick = factor(dist_bin_tick, levels = bin_levels))

p4 <- ggplot(plot_df4, aes(x = dist_bin_tick, y = abs_l2fc)) +
  geom_boxplot(outlier.alpha = 0.1, width = 0.7) +
  facet_wrap(~ comparison, scales = "free_y") +
  labs(
    title = title_base,
    x = x_lab,
    y = y_lab
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )



print(p4)

ggsave(
  paste0(
    plots_folder,
    "/closest_mRNA.cDMRs_diff",
    DIFF_THRESH,
    ".absDmrDistance_vs_absL2FC.boxplots.pdf"
  ),
  plot = p4,
  width = 10,
  height = 7
)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ===== Plot boxplots per bin of ABSOLUTE distance (intergenic cDMRs) ===
# ========================== with fit line ==============================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot details
title_base <- "cDMR absolute distance to closest gene vs differential expression"
Bin_size_bp <- 1000
x_lab <- paste0("Absolute cDMR distance to gene (bp; 0 as single bin, then ", Bin_size_bp, "-bp bins)")
y_lab <- "|log2 fold change|"

XMAX <- 15000

plot_df4 <- plot_df3 %>%
  filter(!is.na(abs_dist_bp), abs_dist_bp <= XMAX) %>%
  mutate(
    dist_bin_start = case_when(
      abs_dist_bp == 0 ~ 0,
      abs_dist_bp > 0  ~ ((abs_dist_bp - 1) %/% Bin_size_bp) * Bin_size_bp + 1
    ),
    dist_bin_end = case_when(
      abs_dist_bp == 0 ~ 0,
      abs_dist_bp > 0  ~ dist_bin_start + Bin_size_bp - 1
    ),
    dist_bin_tick = case_when(
      abs_dist_bp == 0 ~ "0",
      abs_dist_bp > 0  ~ as.character(dist_bin_end)
    ),
    dist_bin_x = case_when(
      abs_dist_bp == 0 ~ 0,
      abs_dist_bp > 0  ~ dist_bin_end
    )
  )

# order bins numerically
bin_levels <- plot_df4 %>%
  distinct(dist_bin_start, dist_bin_tick) %>%
  arrange(dist_bin_start) %>%
  pull(dist_bin_tick)

plot_df4 <- plot_df4 %>%
  mutate(dist_bin_tick = factor(dist_bin_tick, levels = bin_levels))

# numeric x-position for smoothing
plot_df4 <- plot_df4 %>%
  mutate(
    dist_bin_x = case_when(
      abs_dist_bp == 0 ~ 0,
      abs_dist_bp > 0  ~ dist_bin_end
    )
  )

# one summarized value per bin per comparison for the smooth curve
smooth_df <- plot_df4 %>%
  group_by(comparison, dist_bin_x) %>%
  summarise(
    abs_l2fc_median = median(abs_l2fc, na.rm = TRUE),
    # abs_l2fc_median = mean(abs_l2fc, na.rm = TRUE),
    .groups = "drop"
  )
p4 <- ggplot(plot_df4, aes(x = dist_bin_tick, y = abs_l2fc)) +
  geom_boxplot(outlier.alpha = 0.1, width = 0.7) +
  # geom_smooth(
  #   data = smooth_df,
  #   aes(x = dist_bin_x, y = abs_l2fc_median, group = 1),
  #   inherit.aes = FALSE,
  #   method = "loess",
  #   se = FALSE,
  #   linewidth = 0.8,
  #   color = "red",
  #   span = 0.15
  # ) +
  facet_wrap(~ comparison, scales = "free_y") +
  # facet_wrap(~ comparison) +
  labs(
    title = title_base,
    x = x_lab,
    y = y_lab
  ) +
  # scale_x_log10()+
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p4)

ggsave(
  paste0(
    plots_folder,
    "/closest_mRNA.cDMRs_diff",
    DIFF_THRESH,
    ".absDmrDistance_vs_absL2FC.boxplots.withSmooth.pdf"
  ),
  plot = p4,
  width = 10,
  height = 7
)
