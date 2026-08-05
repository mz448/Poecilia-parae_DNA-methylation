#!/usr/bin/env Rscript
# DATE:       2025-09-29
# AUTHOR:     MZF & ChatGPT
# SCRIPT:     02_DEseq_02_geneNamesLabeling_V01.R
# VERSION:    01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Read the combined DE table for Muscle and label/overwrite the gene_name
#   column using the mapping in gene_names.tsv. Write a labeled TSV.
#   Inputs: - Muscle_real-and-shuffled_combined_DEGs.tsv, 
#           - gene_names.tsv
#   Output: Muscle_real-and-shuffled_combined_DEGs_labeled.tsv
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# ------------------------- 1) FILE PATHS --------------------------------------
in_de_path   <- "Muscle_real-and-shuffled_combined_DEGs_CPMfilter.tsv"
names_path   <- "gene_names.tsv"
out_de_path  <- "Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv"

# ------------------------- 2) LOAD -------------------------------------------
de  <- read_tsv(in_de_path, show_col_types = FALSE)
gn  <- read_tsv(names_path, show_col_types = FALSE) %>%
  distinct(gene_id, .keep_all = TRUE)

# ------------------------- 3) LABEL / OVERWRITE -------------------------------
# If input already has gene_name, overwrite missing or placeholder entries with mapping.
# If not present, create it.
de_labeled <- de %>%
  left_join(gn, by = "gene_id", suffix = c("", "_map")) %>%
  mutate(
    gene_name = dplyr::coalesce(.data$gene_name_map, .data$gene_name)  # prefer mapped name
  ) %>%
  select(gene_id, gene_name, everything(), -gene_name_map)

# ------------------------- 4) WRITE -------------------------------------------
write_tsv(de_labeled, out_de_path)
