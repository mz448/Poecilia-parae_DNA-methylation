#!/usr/bin/env Rscript
# DATE:       20250730
# AUTHOR:     Maximiliano Zuluaga Forero
# SCRIPT:     04_descriptivePlots_withinCondition.R
# VERSION:    02

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DESCRIPTION:
# This script filters DMRs by a user-specified condition and generates:
#   - scatter plots (meth1 vs meth2)
#   - histograms of methylation difference (meth1 - meth2)
#   - basic summary statistics
# Adds identity of each group (e.g. Female, Parae-morph) based on comparison name
# All plots are faceted by comparison within the condition.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load libraries
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(tidyr)

# PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dmr_file <- "/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/02_DMR/data/01_methylasso/methylasso_DMRs/dmr_all.tsv"
# selected_condition <- "real_morph"
selected_condition <- "real_sex"
output_dir <- paste0("/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/02_DMR/plots/01_methylasso", selected_condition,"V01")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# LOAD & FILTER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# dmr_all <- read_tsv(dmr_file, show_col_types = FALSE)
# dmr_subset <- dmr_all %>%
#   filter(condition == selected_condition) %>%
#   mutate(length = end - start) %>%
#   rename(chromosome_type = chromosome_type)

dmr_all <- read_tsv(dmr_file, show_col_types = FALSE)

dmr_subset <- dmr_all %>%
  filter(condition == selected_condition) %>%
  mutate(
    DMR_lenght = abs(end - start)
  ) %>%
  rename(chromosome_type = chromosome_type)

# Change the order of the Chromosome type, so the scatter plot is more clear.
dmr_subset <- dmr_subset %>%
  mutate(
    chromosome_type = factor(chromosome_type,
                             levels = c("Sex_Ch", "Autosome"))
  ) %>%
  # arrange(desc(pvalue)) #
  # arrange(desc(chromosome_type)) # plots Sex_Ch first in the scatter
  arrange(chromosome_type) # plots autosome first in the scatter


dmr_subset <- dmr_all %>%
  filter(condition == selected_condition) %>%
  mutate(DMR_lenght = abs(end - start)) %>%
  rename(chromosome_type = chromosome_type)


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
    )
  )

# PLOT 1: SCATTER PLOT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 <- ggplot(dmr_subset, aes(x = meth1, y = meth2, color = chromosome_type)) +
  # geom_point(alpha = 0.3, size = 1, color = "#2C7BB6") +
  # geom_point(alpha = 1/10, size = 2, stroke = "black") +
  
  # geom_point(
  #   aes(fill = chromosome_type),
  #   shape = 21,
  #   color = "black",
  #   stroke = 0.1,
  #   alpha = 1,
  #   size = 0.75
  # ) +
geom_point(
  aes(fill = chromosome_type),
  shape = 21,
  color = "black" ,
  stroke = 0,
  alpha = 1,
  size = 2
) +
  scale_fill_manual(values = c(
    "Sex_Ch" = "red",
    "Autosome" = "blue"
  ))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~comparison, scales = "fixed") +
  coord_flip() +
  labs(
    title = paste("Methylation1 vs Methylation2 -", selected_condition),
    # title = paste0(group1_label, "vs", group2_label, "-", selected_condition),
    x = "Methylation Level Group 1 (%)",
    y = "Methylation Level Group 2 (%)"
  ) +
  theme_minimal()+
  theme(
    legend.position = "top",
    legend.direction = "horizontal"
  )
print (p1)
ggsave(filename = file.path(output_dir, paste0("scatter_meth1_vs_meth2_", selected_condition, ".pdf")),
       plot = p1, width = 10, height = 5)

# PLOT 2: HISTOGRAM OF DIFF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p2 <- ggplot(dmr_subset, aes(x = diff)) +
  geom_histogram(binwidth = 3, fill = "#00a651", linewidth = 0.1, color = "black") +
  geom_vline(xintercept = 0 ,linetype = "dashed", color = "black") +
  facet_wrap(~comparison, scales = "fixed") +
  # scale_x_continuous(breaks = seq(-100, 100, by = 25)) +
  scale_x_reverse(
    limits = c(100, -100),
    breaks = seq(100, -100, by = -25)
  )+
  scale_y_continuous(breaks = seq(0, 900, by = 100)) +
  coord_cartesian(xlim = c(-100, 100)) +
  labs(
    title = paste("Distribution of Methylation Difference -", selected_condition),
    # title = paste0("Distribution of Methylation Difference - ", group1_label, "vs", group2_label, " - ", selected_condition),
    x = "% of Methylation difference",
    y = "DMR count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print (p2)
ggsave(filename = file.path(output_dir, paste0("hist_diff_", selected_condition, ".pdf")),
       plot = p2, width = 10, height = 5)


# PLOT 3: DMR LENGTH BY METHYLATION DIFFERENCE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# X axis = signed methylation difference, fixed from -100 to 100
# Y axis = DMR length in bp




p3 <- ggplot(dmr_subset, aes(x = diff, y = DMR_lenght)) +
  geom_point(
    alpha = 0.3,
    size = 1,
    color = "#00a651"
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~comparison, scales = "fixed") +
  scale_y_continuous(
    breaks = seq(
      0,
      ceiling(max(dmr_subset$DMR_lenght, na.rm = TRUE) / 1000) * 1000,
      by = 1000
    ),
    minor_breaks = seq(
      0,
      ceiling(max(dmr_subset$DMR_lenght, na.rm = TRUE) / 500) * 500,
      by = 500
    )
  )+
  scale_x_reverse(
    limits = c(100, -100),
    breaks = seq(100, -100, by = -10)
  ) +
  coord_cartesian(xlim = c(-100, 100)) +
  labs(
    title = paste("DMR length by methylation difference -", selected_condition),
    x = "% of methylation difference",
    y = "DMR length (bp)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

print(p3)

ggsave(
  filename = file.path(output_dir, paste0("DMR_length_by_diff_", selected_condition, ".pdf")),
  plot = p3,
  width = 10,
  height = 5
)

# STATS TABLE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
stats_summary <- dmr_subset %>%
  group_by(comparison, group1_label, group2_label) %>%
  summarise(
    n_DMRs = n(),
    mean_diff = mean(diff, na.rm = TRUE),
    median_diff = median(diff, na.rm = TRUE),
    mean_meth1 = mean(meth1, na.rm = TRUE),
    mean_meth2 = mean(meth2, na.rm = TRUE),
    mean_DMR_lenght = mean(DMR_lenght, na.rm = TRUE),
    mean_cpgs = mean((num.cpgs1 + num.cpgs2) / 2, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(stats_summary, file.path(output_dir, paste0("summary_stats-", selected_condition, ".tsv")))
write_tsv(dmr_subset, file.path(output_dir, paste0("dmr_list-", selected_condition, ".tsv")))
cat("✅ All plots and summaries saved for condition:", selected_condition, "\n")
