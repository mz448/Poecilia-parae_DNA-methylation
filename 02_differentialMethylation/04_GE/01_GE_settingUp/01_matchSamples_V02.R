#!/usr/bin/env Rscript
# DATE:       2025-09-27
# AUTHOR:     MZF & ChatGPT
# SCRIPT:     01_matchSamples_V02.R
# VERSION:    02
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:       Reproduce methylation real/shuffled group assignments in RNA-seq colData for DESeq2
#             Add one column per methylation 'comparison' to sample_info.tsv.
#             Each new column has the group label for that sample in that comparison,
#             or "FALSE" if the sample is not included.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%


library(readr)
library(dplyr)
library(tidyr)
library(stringr)



# 0) input ------------------------------------------------------------------------
# Edit paths as needed
methyl_path  <- "methyl.tsv"                # columns: condition, comparison, group, Rep1, Rep2, Rep3
sample_path  <- "Sample_Info.tsv"           # columns include: Sample, Rep, Species, Morph, Sex, Tissue, Treatment, ...
map_path     <- "rnaseq_to_methyl_map.tsv"  # two columns: RNASample, methylationSample

# Output path (non-destructive)
out_path     <- "sample_info_with_comparisons.tsv"

# 1) Read inputs ---------------------------------------------------------------
methyl <- read_tsv(methyl_path, show_col_types = FALSE)
sample_info <- read_tsv(sample_path, show_col_types = FALSE)
rna2meth <- read_tsv(map_path, show_col_types = FALSE,
                     col_names = c("RNASample","methylationSample"))

# 2) Expand methyl to per-sample rows ------------------------------------------
# From Rep1/Rep2/Rep3 -> long 'methylationSample'
methyl_long <- methyl %>%
  pivot_longer(cols = starts_with("Rep"),
               names_to = "rep_label",
               values_to = "methylationSample") %>%
  filter(!is.na(methylationSample), methylationSample != "") %>%
  select(comparison, group, methylationSample)

# 3) Map methyl IDs → RNA IDs --------------------------------------------------
# Keep only mapped samples; any unmapped methyl IDs are dropped.
assignments_rna <- methyl_long %>%
  left_join(rna2meth, by = "methylationSample") %>%
  filter(!is.na(RNASample)) %>%
  transmute(Sample = RNASample, comparison, group)

# 4) Pivot to wide: one column per comparison ----------------------------------
# For each Sample, for each comparison, value = group
wide_assign <- assignments_rna %>%
  distinct(Sample, comparison, group) %>%
  pivot_wider(names_from = comparison, values_from = group)

# 5) Join to sample_info and fill non-participants -----------------------------
# Left-join preserves the exact rows/order of sample_info.
dkey <- sample_info %>%
  left_join(wide_assign, by = c("Sample"))

# Identify the comparison columns we just added
comp_cols <- setdiff(names(dkey), names(sample_info))

# Replace NA in comparison columns with "FALSE"
if (length(comp_cols) > 0) {
  dkey <- dkey %>%
    mutate(across(all_of(comp_cols),
                  ~ ifelse(is.na(.x) | .x == "", "FALSE", .x)))
}

# 6) Write out -----------------------------------------------------------------
write_tsv(dkey, out_path)
