#!/usr/bin/env Rscript
# DATE:       2026-03-22
# AUTHOR:     MZF
# SCRIPT:     01_plot_bismark_summary_contexts.R
# VERSION:    01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Read a Bismark summary table and generate four ggplot figures:
#
#   1) One stacked bar per sample showing the cumulative percentage of cytosines
#      in CpG, CHG, CHH, and Other contexts.
#   2) One single stacked bar showing the same context composition averaged across
#      all samples.
#   3) One bar plot per sample showing methylation percentage in each context
#      (CpG, CHG, CHH).
#   4) One bar plot showing mean methylation percentage in each context across
#      all samples.
#   5) One bar plot with mean methylation percentag in each contex across 
#      all samples. Use a single axis and report mean and error (SEM)
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# INPUT:
#   A tab-delimited Bismark summary file with at least these columns:
#     File
#     Total Cs
#     Methylated CpGs
#     Unmethylated CpGs
#     Methylated chgs
#     Unmethylated chgs
#     Methylated CHHs
#     Unmethylated CHHs
#
# MAIN CALCULATIONS:
#   - Total context-specific Cs:
#       CpG = Methylated CpGs + Unmethylated CpGs
#       CHG = Methylated chgs + Unmethylated chgs
#       CHH = Methylated CHHs + Unmethylated CHHs
#       Other = Total Cs - (CpG + CHG + CHH)
#
#   - Context composition (% of all Cs):
#       100 * context_total / Total Cs
#
#   - Context-specific methylation (%):
#       100 * methylated_context / (methylated_context + unmethylated_context)
#
# NOTES:
#   - "Other" represents cytosines not accounted for in CpG, CHG, or CHH classes
#     in the summary table.
#   - Sample names are cleaned by removing common Bismark/Bowtie2 suffixes.
#   - The script assumes one row per sample.
#
# OUTPUT:
#   Four PDF files saved in the working directory:
#     01_context_composition_per_sample.pdf
#     02_context_composition_mean_all_samples.pdf
#     03_methylation_percent_per_sample.pdf
#     04_methylation_percent_mean_all_samples.pdf
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggbreak) 
})
#install.packages("ggbreak")
# To start fresh ---------------------------------------------------------------

graphics.off()
rm(list = ls())

#===============================================================================
# BLOCK 1: LOAD TABLE ----------------------------------------------------------
#===============================================================================

input_file <-   "/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bismarkAlignment/bismark/bismark/summary/bismark_summary_report.txt"
plots_folder <- "/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/plots/"

dat <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "File",
  "Total Cs",
  "Methylated CpGs",
  "Unmethylated CpGs",
  "Methylated chgs",
  "Unmethylated chgs",
  "Methylated CHHs",
  "Unmethylated CHHs"
)

missing_cols <- setdiff(required_cols, colnames(dat))
if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing from the input table:\n",
    paste(missing_cols, collapse = ", ")
  )
}

# Ensure numeric columns are truly numeric
numeric_cols <- setdiff(required_cols, "File")
dat[numeric_cols] <- lapply(dat[numeric_cols], function(x) as.numeric(gsub(",", "", x)))

#===============================================================================
# BLOCK 2: CLEAN SAMPLE NAMES --------------------------------------------------
#===============================================================================
sample_order <- c(
  "pparfmem001",
  "pparfmem002",
  "pparfmem003",
  "pparfmem004",
  "pparfmem005",
  "pparfmem006",
  "pparfmem007",
  "pparfmem008",
  "pparfmem009",
  "pparfmem010",
  "pparfmem011",
  "pparfmem012"
)


# dat <- dat %>%
#   mutate(
#     Sample = File,
#     Sample = str_remove(Sample, "_1_val_1_bismark_bt2_pe\\.bam$"),
#     Sample = str_trim(Sample)
#   )

dat2 <- dat %>%
  mutate(
    Sample = File,
    Sample = str_remove(Sample, "_1_val_1_bismark_bt2_pe\\.bam$"),
    Sample = str_trim(Sample),
    Sample = factor(Sample, levels = sample_order)
  )

# Preserve file order in plots
# dat$Sample <- factor(dat$Sample, levels = dat$Sample)

#===============================================================================
# BLOCK 3: CALCULATIONS --------------------------------------------------------
#===============================================================================

dat_calc <- dat2 %>%
  mutate(
    CpG_total = `Methylated CpGs` + `Unmethylated CpGs`,
    CHG_total = `Methylated chgs` + `Unmethylated chgs`,
    CHH_total = `Methylated CHHs` + `Unmethylated CHHs`,
    Other_total = `Total Cs` - CpG_total - CHG_total - CHH_total,
    Other_total = ifelse(Other_total < 0, 0, Other_total),
    
    CpG_pct_of_allCs = 100 * CpG_total / `Total Cs`,
    CHG_pct_of_allCs = 100 * CHG_total / `Total Cs`,
    CHH_pct_of_allCs = 100 * CHH_total / `Total Cs`,
    Other_pct_of_allCs = 100 * Other_total / `Total Cs`,
    
    CpG_meth_pct = ifelse(CpG_total > 0, 100 * `Methylated CpGs` / CpG_total, NA),
    CHG_meth_pct = ifelse(CHG_total > 0, 100 * `Methylated chgs` / CHG_total, NA),
    CHH_meth_pct = ifelse(CHH_total > 0, 100 * `Methylated CHHs` / CHH_total, NA)
  )

# Long table for Figure 1
context_comp_long <- dat_calc %>%
  select(
    Sample,
    CpG_pct_of_allCs,
    CHG_pct_of_allCs,
    CHH_pct_of_allCs,
    Other_pct_of_allCs
  ) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Context",
    values_to = "Percent"
  ) %>%
  mutate(
    Context = recode(
      Context,
      CpG_pct_of_allCs = "CpG",
      CHG_pct_of_allCs = "CHG",
      CHH_pct_of_allCs = "CHH",
      Other_pct_of_allCs = "Other"
    ),
    Context = factor(Context, levels = c("CpG", "CHG", "CHH", "Other"))
  )

# Summary table for Figure 2
context_comp_mean <- dat_calc %>%
  summarise(
    CpG = mean(CpG_pct_of_allCs, na.rm = TRUE),
    CHG = mean(CHG_pct_of_allCs, na.rm = TRUE),
    CHH = mean(CHH_pct_of_allCs, na.rm = TRUE),
    Other = mean(Other_pct_of_allCs, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Context",
    values_to = "Percent"
  ) %>%
  mutate(
    SampleGroup = "Mean across samples",
    Context = factor(Context, levels = c("CpG", "CHG", "CHH", "Other"))
  )

context_comp_mean <- context_comp_mean %>%
  filter (Context != "Other")

# Long table for Figure 3
meth_context_long <- dat_calc %>%
  select(Sample, CpG_meth_pct, CHG_meth_pct, CHH_meth_pct) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Context",
    values_to = "MethylationPercent"
  ) %>%
  mutate(
    Context = recode(
      Context,
      CpG_meth_pct = "CpG",
      CHG_meth_pct = "CHG",
      CHH_meth_pct = "CHH"
    ),
    Context = factor(Context, levels = c("CpG", "CHG", "CHH"))
  )

# Summary table for Figure 4
meth_context_mean <- dat_calc %>%
  summarise(
    CpG = mean(CpG_meth_pct, na.rm = TRUE),
    CHG = mean(CHG_meth_pct, na.rm = TRUE),
    CHH = mean(CHH_meth_pct, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Context",
    values_to = "MethylationPercent"
  )

#===============================================================================
# BLOCK 4: PLOTS (ggplot)
#===============================================================================

# context_colors <- c(
#   "CpG"   = "#1b9e77",
#   "CHG"   = "#d95f02",
#   "CHH"   = "#7570b3",
#   "Other" = "#666666"
# )

context_colors <- c(
  "CpG"   = "#E4312b",
  "CHG"   = "#149954",
  "CHH"   = "#0096FF",
  "Other" = "#666666"
)

breaks = c("CpG", "CHG", "CHH", "Other")


base_theme <- theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

#--------------------------------------------------------------------------
# FIGURE 1: stacked context composition per sample ------------------------
#--------------------------------------------------------------------------

p1 <- ggplot(context_comp_long, aes(x = Sample, y = Percent, fill = Context)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = context_colors,
    breaks = c("CpG", "CHG", "CHH", "Other"),
    drop = FALSE
  )+
  labs(
    title = "Cytosine context composition per sample",
    x = "Sample",
    y = "% of all cytosines",
    fill = "Context"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  base_theme
print (p1)
ggsave(
  filename = paste0(plots_folder,"01_context_composition_per_sample.pdf"),
  plot = p1,
  width = 11,
  height = 5
)

#--------------------------------------------------------------------------
# FIGURE 2: stacked mean context composition across all samples -----------
#--------------------------------------------------------------------------


p2 <- ggplot(context_comp_mean, aes(x = SampleGroup, y = Percent, fill = Context)) +
  geom_col(width = 0.6) +
  scale_fill_manual(
    values = context_colors,
    breaks = c("CpG", "CHG", "CHH", "Other"),
    drop = FALSE
  )+
  geom_text(
    aes(label = sprintf("%.1f%%", Percent)),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white"
  ) +
  labs(
    title = "Cumulative % of all cytosines",
    x = NULL,
    y = "% of all cytosines",
    fill = "Context"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )
print (p2)
# ggsave(
#   filename = paste0(plots_folder,"02_context_composition_mean_all_samples.pdf"),
#   plot = p2,
#   width = 5,
#   height = 6
# )
ggsave(
  filename = paste0(plots_folder,"02_context_composition_mean_all_samples.pdf"),
  plot = p2,
  width = 3,
  height = 4
)

#--------------------------------------------------------------------------
# FIGURE 3: methylation percentage by context for each sample -------------
#--------------------------------------------------------------------------
library(ggbreak) 

p3 <- ggplot(meth_context_long, aes(x = Sample, y = MethylationPercent, fill = Context)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    values = context_colors,
    breaks = c("CpG", "CHG", "CHH", "Other"),
    drop = FALSE
  )+
  labs(
    title = "Cytosine methylation percentage by context per sample",
    x = "Sample",
    y = "% mC/C",
    fill = "Context"
  ) +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0)) +
  scale_y_break(c(1, 70)) +
  # scale_y_break(c(75, 100)) +
  base_theme
print (p3)
ggsave(
  filename = paste0(plots_folder,"03_methylation_percent_per_sample.pdf"),
  plot = p3,
  width = 11,
  height = 6
)

#--------------------------------------------------------------------------
# FIGURE 4: mean methylation percentage by context across all samples -----
#--------------------------------------------------------------------------

p4 <- ggplot(meth_context_mean, aes(x = Context, y = MethylationPercent, fill = Context)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = context_colors,
    breaks = c("CpG", "CHG", "CHH", "Other"),
    drop = FALSE
  )+
  geom_text(
    aes(label = sprintf("%.2f%%", MethylationPercent)),
    position = position_stack(vjust = 0.5),
    size = 4,
    # color = "white"
  ) +
  labs(
    title = "Mean cytosine methylation % by context across all samples",
    x = "Context",
    y = "% mC/C",
    fill = "Context"
  ) +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0)) +
  theme_classic(base_size = 12) +
  scale_y_break(c(1.0, 70))
# scale_y_break(c(75, 100)) +
# theme(
#   panel.grid.major.x = element_blank(),
#   panel.grid.minor = element_blank()
# )
print (p4)
ggsave(
  filename = paste0(plots_folder,"04_methylation_percent_mean_all_samples.pdf"),
  plot = p4,
  width = 6,
  height = 6
)

#--------------------------------------------------------------------------
# Fasceted mean methylation percentage by context across all samples
#--------------------------------------------------------------------------
meth_context_mean2 <- meth_context_mean %>%
  mutate(
    AxisGroup = ifelse(Context == "CpG", "CpG", "CHG + CHH")
  )

p5 <- ggplot(meth_context_mean2, aes(x = Context, y = MethylationPercent, fill = Context)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = context_colors,
    breaks = c("CpG", "CHG", "CHH", "Other"),
    drop = FALSE
  )+
  geom_text(
    aes(label = sprintf("%.2f%%", MethylationPercent)),
    vjust = 2,
    size = 4,
    color = "white"
  ) +
  # facet_wrap(~ AxisGroup, nrow = 1)+
  facet_wrap(~ AxisGroup, scales = "free_x", nrow = 1)+
  labs(
    title = "Mean cytosine methylation % by context across all samples",
    x = "Context",
    y = "% mC/C",
    fill = "Context"
  ) +
  # scale_y_continuous(limits = c(0, 80), expand = c(0, 0)) +
  coord_flip()+
  theme_classic(base_size = 12) 
# theme(
#   panel.grid.major.x = element_blank(),
#   panel.grid.minor = element_blank(),
#   strip.background = element_rect(fill = "grey90"),
#   strip.text = element_text(face = "bold")
# )
# base_theme

print(p5)
ggsave(
  filename = paste0(plots_folder,"04_methylation_percent_mean_all_samples_fascet.pdf"),
  plot = p5,
  width = 6,
  height = 6
)

#--------------------------------------------------------------------------
# FIGURE 5: mean methylation percentage by context with SEM ----------------
#--------------------------------------------------------------------------

meth_context_summary <- meth_context_long %>%
  group_by(Context) %>%
  summarise(
    mean_meth_pct = mean(MethylationPercent, na.rm = TRUE),
    sd_meth_pct   = sd(MethylationPercent, na.rm = TRUE),
    n             = sum(!is.na(MethylationPercent)),
    sem_meth_pct  = sd_meth_pct / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    Context = factor(Context, levels = c("CpG", "CHG", "CHH"))
  )

p6 <- ggplot(
  meth_context_summary,
  aes(x = Context, y = mean_meth_pct, fill = Context)
) +
  geom_col(width = 0.7) +
  geom_errorbar(
    aes(
      ymin = mean_meth_pct - sem_meth_pct,
      ymax = mean_meth_pct + sem_meth_pct
    ),
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      y = mean_meth_pct + sem_meth_pct,
      label = sprintf("%.2f%% ± %.2f", mean_meth_pct, sem_meth_pct)
    ),
    vjust = -0.6,
    size = 4
  ) +
  scale_fill_manual(
    values = context_colors,
    breaks = c("CpG", "CHG", "CHH"),
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, max(meth_context_summary$mean_meth_pct + meth_context_summary$sem_meth_pct, na.rm = TRUE) * 1.12),
    expand = c(0, 0)
  ) +
  labs(
    title = "Mean cytosine methylation % by context across all samples",
    x = "Context",
    y = "% mC/C",
    fill = "Context"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

print(p6)

ggsave(
  filename = paste0(plots_folder, "05_methylation_percent_mean_SEM_all_samples.pdf"),
  plot = p6,
  width = 4,
  height = 4
)

