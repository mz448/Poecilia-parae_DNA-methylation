#!/usr/bin/env Rscript
# DATE:       2026-05-06
# AUTHOR:     MZF
# SCRIPT:     03_Summary_DML_levels.R
# VERSION:    01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#             Generate DML summary at the whole chromosome and window scale
# DETAILS: 
#              - The weighted percentage is the chromosome/partition-level DML 
#                rate and should be stable.  
#              - The mean of window percentages is a description of local 
#                spatial variation and is expected to change with window size.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# User settings
# ==============================================================================
output_dir <- "DNMTools_DML_summary"

input_file <- file.path(output_dir,"12_DNMTools_window_level_DML_summary_by_partition.tsv")

 dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Load window-level table
# ==============================================================================

dt <- fread(input_file)

required_cols <- c(
  "contrast",
  "partition",
  "chr",
  "window_start",
  "window_end",
  "covered_CpGs",
  "DMLs_FDR_le_0_01",
  "percentage_DMLs"
)

missing_cols <- setdiff(required_cols, names(dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# ==============================================================================
# 1. Descriptive statistics per contrast, partition, and chromosome
# ==============================================================================

chromosome_window_stats <- dt[
  ,
  .(
    n_windows = .N,
    total_covered_CpGs = sum(covered_CpGs, na.rm = TRUE),
    total_DMLs = sum(DMLs_FDR_le_0_01, na.rm = TRUE),
    
    # Weighted chromosome-level DML percentage
    weighted_percentage_DMLs = 100 *
      sum(DMLs_FDR_le_0_01, na.rm = TRUE) /
      sum(covered_CpGs, na.rm = TRUE),
    
    # Distribution of window-level percentages
    mean_window_percentage_DMLs = mean(percentage_DMLs, na.rm = TRUE),
    sd_window_percentage_DMLs = sd(percentage_DMLs, na.rm = TRUE),
    sem_window_percentage_DMLs = sd(percentage_DMLs, na.rm = TRUE) / sqrt(.N),
    median_window_percentage_DMLs = median(percentage_DMLs, na.rm = TRUE),
    min_window_percentage_DMLs = min(percentage_DMLs, na.rm = TRUE),
    q25_window_percentage_DMLs = quantile(percentage_DMLs, 0.25, na.rm = TRUE),
    q75_window_percentage_DMLs = quantile(percentage_DMLs, 0.75, na.rm = TRUE),
    max_window_percentage_DMLs = max(percentage_DMLs, na.rm = TRUE)
  ),
  by = .(contrast, partition, chr)
]

chromosome_window_stats[, chr_number := suppressWarnings(as.integer(sub("^Parae_", "", chr)))]

chromosome_window_stats[, partition_order := fifelse(
  partition == "WG", 1L,
  fifelse(partition == "AUTO", 2L,
          fifelse(partition == "CH12", 3L, 99L))
)]

setorder(
  chromosome_window_stats,
  contrast,
  partition_order,
  chr_number
)

chromosome_window_stats[, c("partition_order", "chr_number") := NULL]

chromosome_stats_file <- file.path(
  output_dir,
  "18_DNMTools_window_percentage_DML_stats_by_chromosome_partition.tsv"
)

fwrite(
  chromosome_window_stats,
  chromosome_stats_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Saved chromosome-level window statistics to: ", chromosome_stats_file)

# ==============================================================================
# 2. Descriptive statistics per contrast and partition
# ==============================================================================

partition_window_stats <- dt[
  ,
  .(
    n_windows = .N,
    n_chromosomes = uniqueN(chr),
    total_covered_CpGs = sum(covered_CpGs, na.rm = TRUE),
    total_DMLs = sum(DMLs_FDR_le_0_01, na.rm = TRUE),
    
    # Weighted partition-level DML percentage
    weighted_percentage_DMLs = 100 *
      sum(DMLs_FDR_le_0_01, na.rm = TRUE) /
      sum(covered_CpGs, na.rm = TRUE),
    
    # Distribution of window-level percentages
    mean_window_percentage_DMLs = mean(percentage_DMLs, na.rm = TRUE),
    sd_window_percentage_DMLs = sd(percentage_DMLs, na.rm = TRUE),
    sem_window_percentage_DMLs = sd(percentage_DMLs, na.rm = TRUE) / sqrt(.N),
    median_window_percentage_DMLs = median(percentage_DMLs, na.rm = TRUE),
    min_window_percentage_DMLs = min(percentage_DMLs, na.rm = TRUE),
    q25_window_percentage_DMLs = quantile(percentage_DMLs, 0.25, na.rm = TRUE),
    q75_window_percentage_DMLs = quantile(percentage_DMLs, 0.75, na.rm = TRUE),
    max_window_percentage_DMLs = max(percentage_DMLs, na.rm = TRUE)
  ),
  by = .(contrast, partition)
]

partition_window_stats[, partition_order := fifelse(
  partition == "WG", 1L,
  fifelse(partition == "AUTO", 2L,
          fifelse(partition == "CH12", 3L, 99L))
)]

setorder(
  partition_window_stats,
  contrast,
  partition_order
)

partition_window_stats[, partition_order := NULL]

partition_stats_file <- file.path(
  output_dir,
  "19_DNMTools_window_percentage_DML_stats_by_partition.tsv"
)

fwrite(
  partition_window_stats,
  partition_stats_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Saved partition-level window statistics to: ", partition_stats_file)

# ==============================================================================
# 3. Print outputs
# ==============================================================================

print(chromosome_window_stats)
print(partition_window_stats)
