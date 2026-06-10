# DATE:       20250730
# AUTHOR:     Maximiliano Zuluaga Forero
# SCRIPT:     04_descriptivePlots_withinCondition.R
# VERSION:    03

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DESCRIPTION:
# This script filters DMRs by a user-specified condition and generates:
#   - scatter plots (meth1 vs meth2), colored by annotation
#   - histograms of methylation difference (meth1 - meth2), with density
#   - basic summary statistics
# Adds identity of each group based on comparison name
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dmr_file <- "/local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_differential_methylation/03_methylasso_DMRs/data/dmr_output/q0.05_cov3_sym_run/dmr_all.tsv"
selected_condition <- "real_sex"
# selected_condition <- "real_morph"
# selected_condition <- "shuffled_sex"
# selected_condition <- "shuffled_morph"
# selected_condition <- "shuffled_all"
output_dir <- paste0("/local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_differential_methylation/03_methylasso_DMRs/plots/dmr_output/q0.05_cov3_sym_run/", selected_condition)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# LOAD & FILTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dmr_all <- read_tsv(dmr_file, show_col_types = FALSE)
dmr_subset <- dmr_all %>%
  filter(condition == selected_condition) %>%
  mutate(length = end - start)

# EXTRACT GROUP NAMES FROM COMPARISON ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dmr_subset <- dmr_subset %>%
  separate(comparison, into = c("group1", "group2"), sep = "_vs_", remove = FALSE) %>%
  mutate(
    group1_label = case_when(
      group1 %in% c("f", "rf")        ~ "Female",
      group1 %in% c("i", "i_01")      ~ "Immaculata-morph",
      group1 %in% c("y", "y_01")      ~ "Yellow-morph",
      group1 %in% c("p", "p_01")      ~ "Parae-morph",
      str_detect(group1, "^sM")       ~ "Shuffled males",
      str_detect(group1, "^sA")       ~ "Shuffled all",
      TRUE                            ~ "Unknown"
    ),
    group2_label = case_when(
      group2 %in% c("f", "rf")        ~ "Female",
      group2 %in% c("i", "i_01")      ~ "Immaculata-morph",
      group2 %in% c("y", "y_01")      ~ "Yellow-morph",
      group2 %in% c("p", "p_01")      ~ "Parae-morph",
      str_detect(group2, "^sM")       ~ "Shuffled males",
      str_detect(group2, "^sA")       ~ "Shuffled all",
      TRUE                            ~ "Unknown"
    ),
    comparison_label = paste0(group1_label, " vs ", group2_label)
  )


# colors <- c("Female" = "#149954", 
#             "Immaculata-morph" = "#888888", 
#             "Parae-morph" = "#9E9AC6", 
#             "Yellow-morph" = "#F68C41", 
#             "Shuffled males" = "#991459", 
#             "Shuffled all" = "#E4312b")

# # Set individual colors for each plot ##########################################
# ## REAL_sex
# colors <- c(
#   f_vs_i = "#888888", 
#   f_vs_p = "#9E9AC6", 
#   f_vs_y = "#F68C41")
# plotting_function("REAL_sex",colors)
# 
# ## REAL_morph
# colors <- c(
#   p_vs_i_01 = "#000000", 
#   y_vs_p_01 = "#000000", 
#   i_vs_y_01 = "#000000", 
#   i_vs_p_02 = "#000000", 
#   p_vs_y_02 = "#000000", 
#   y_vs_i_02 = "#000000")
# plotting_function("REAL_morph",colors)
# 
# ## SHUFFLED_sex
# colors <- c(
#   rf_vs_sM_01 = "#149954", 
#   rf_vs_sM_02 = "#149954", 
#   rf_vs_sM_03 = "#149954", 
#   rf_vs_sM_04 = "#149954", 
#   rf_vs_sM_05 = "#149954", 
#   rf_vs_sM_06 = "#149954")
# plotting_function("SHUFFLED_sex",colors)
# 
# ## SHUFFLED_morph
# colors <- c(
#   sM_vs_sM_01 = "#991459", 
#   sM_vs_sM_02 = "#991459", 
#   sM_vs_sM_03 = "#991459"
# )
# plotting_function("SHUFFLED_morph",colors)
# 
# ## SHUFFLED_all
# colors <- c(
#   sA_vs_sA_01 = "#E4312b", 
#   sA_vs_sA_02 = "#E4312b", 
#   sA_vs_sA_03 = "#E4312b")
# plotting_function("SHUFFLED_all",colors)


# scale_colour_viridis_d(option = "inferno") +

  
  
# Define chromosome to highlight
highlight_chr <- "Parae_12"

# Add a new logical column to identify highlight
dmr_subset <- dmr_subset %>%
  mutate(highlight_chr = chr == highlight_chr)

# PLOT 1: SCATTER PLOT COLORED BY ANNOTATION ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 <- ggplot(dmr_subset, aes(x = meth1, y = meth2, color = annotation)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~comparison, scales = "fixed") +
  labs(
    title = paste("Methylation1 vs Methylation2 -", selected_condition),
    x = "Methylation Level Group 1 (%)",
    y = "Methylation Level Group 2 (%)",
    color = "Annotation"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, paste0("scatter_meth1_vs_meth2_", selected_condition,"_annot", ".pdf")),
       plot = p1, width = 10, height = 6)

# PLOT 3: SCATTER PLOT COLORED BY CHROMOSOME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p3 <- ggplot(dmr_subset, aes(x = meth1, y = meth2, color = chr)) +
  geom_point(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color ="black") +
  # scale_colour_viridis_d(option = "inferno") +
  scale_colour_viridis_d("plasma") +
  facet_wrap(~comparison, scales = "fixed") +
  labs(
    title = paste("Methylation1 vs Methylation2 -", selected_condition),
    x = "Methylation Level Group 1 (%)",
    y = "Methylation Level Group 2 (%)",
    color = "Chromosome"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, paste0("scatter_meth1_vs_meth2_", selected_condition,"_chr", ".pdf")),
       plot = p3, width = 10, height = 6)


# PLOT 4: SCATER PLOT HIGHLIGHTING 12 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plot with fill as chromosome color, outline based on highlight
p4 <- ggplot(dmr_subset, aes(x = meth1, y = meth2)) +
  geom_point(
    aes(fill = chr, color = highlight_chr),
    size = 1.8,
    shape = 21,
    stroke = 0.8
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~comparison, scales = "fixed") +
  # scale_color_manual(
  #   values = c(`TRUE` = "black", `FALSE` = "transparent"),
  #   guide = "none"
  # ) +
  scale_color_manual(
    values = c(`TRUE` = "black", `FALSE` = "transparent"),
    name = paste("Highlighted\nChr:", highlight_chr)
  )+
  labs(
    title = paste("DMRs -", selected_condition),
    x = "Methylation Level Group 1 (%)",
    y = "Methylation Level Group 2 (%)",
    fill = "Chromosome"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, paste0("scatter_highlight_chr_", highlight_chr, "_", selected_condition, ".pdf")),
       plot = p4, width = 10, height = 6)


# PLOT 35 SCATTER PLOT COLORED BY CHROMOSOME ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p5 <- ggplot(dmr_subset, aes(x = meth1, y = meth2, color = chrosomome_type)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color ="grey") +
  # scale_colour_viridis_d(option = "inferno") +
  # scale_colour_viridis_d("plasma") +
  facet_wrap(~comparison, scales = "fixed") +
  labs(
    title = paste("Chromosome type -", selected_condition),
    x = "Methylation Level Group 1 (%)",
    y = "Methylation Level Group 2 (%)",
    color = "Chromosome"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, paste0("scatter_chrom_type", selected_condition,"_chr", ".pdf")),
       plot = p5, width = 10, height = 6)

# # Plot with fill as chromosome color, outline based on highlight
# p4 <- ggplot(dmr_subset, aes(x = meth1, y = meth2)) +
#   geom_point(
#     aes(fill = chr, color = highlight_chr, alpha = 0.8),
#     size = 1,
#     shape = 21,
#     stroke = 0.5,
#   ) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
#   facet_wrap(~comparison_label, scales = "fixed") +
#   scale_color_manual(
#     values = c(`TRUE` = "black", `FALSE` = "transparent"),
#     guide = "none"
#   ) +
#   labs(
#     # title = paste("Highlight:", highlight_chr, "-", selected_condition),
#     title = paste("DMRS -", selected_condition),
#     x = "Methylation Level Group 1 (%)",
#     y = "Methylation Level Group 2 (%)",
#     fill = "Chromosome"
#   ) +
#   theme_minimal() +
#   theme(legend.position = "bottom")
# 
# ggsave(filename = file.path(output_dir, paste0("scatter_highlight_chr_", highlight_chr, "_", selected_condition, ".pdf")),
#        plot = p4, width = 10, height = 6)

# # PLOT 2: HISTOGRAM OF DIFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# p2 <- ggplot(dmr_subset, aes(x = diff)) +
#   geom_histogram(binwidth = 2, fill = "#FDAE61", color = "black") +
#   facet_wrap(~comparison, scales = "fixed") +
#   labs(
#     title = paste("Distribution of Methylation Difference -", selected_condition),
#     # title = paste0("Distribution of Methylation Difference - ", group1_label, "vs", group2_label, " - ", selected_condition),
#     x = "Methylation Difference (meth1 - meth2)",
#     y = "Count"
#   ) +
#   theme_minimal()
# 
# ggsave(filename = file.path(output_dir, paste0("hist_diff_", selected_condition, ".pdf")),
#        plot = p2, width = 10, height = 6)


# PLOT 2: HISTOGRAM WITH DENSITY ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# p2 <- ggplot(dmr_subset, aes(x = diff, fill = annotation)) +
p2 <- ggplot(dmr_subset, aes(x = diff)) +  
  geom_histogram(binwidth = 2, fill = "#FDAE61", color = "black", alpha = 0.3, position = "identity") +
  geom_density(aes(y = ..count..), color = "black", fill = NA, linetype = "dotted") +
  facet_wrap(~comparison, scales = "fixed") +
  labs(
    title = paste("Distribution of Methylation Difference -", selected_condition),
    x = "Methylation Difference (meth1 - meth2)",
    y = "Count",
    fill = "Annotation"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(filename = file.path(output_dir, paste0("hist_diff_", selected_condition, ".pdf")),
       plot = p2, width = 10, height = 6)

# STATS TABLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stats_summary <- dmr_subset %>%
  group_by(comparison, group1_label, group2_label) %>%
  summarise(
    n_DMRs = n(),
    mean_diff = mean(diff, na.rm = TRUE),
    median_diff = median(diff, na.rm = TRUE),
    mean_meth1 = mean(meth1, na.rm = TRUE),
    mean_meth2 = mean(meth2, na.rm = TRUE),
    mean_length = mean(length, na.rm = TRUE),
    mean_cpgs = mean((num.cpgs1 + num.cpgs2) / 2, na.rm = TRUE),
    .groups = "drop"
  )

# Save outputs
write_tsv(stats_summary, file.path(output_dir, paste0("summary_stats-", selected_condition, ".tsv")))
write_tsv(dmr_subset, file.path(output_dir, paste0("dmr_list-", selected_condition, ".tsv")))
cat("✅ Enhanced plots and summaries saved for condition:", selected_condition, "\n")
