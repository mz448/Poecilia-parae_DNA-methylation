#!/usr/bin/env Rscript
# DATE:       2026-05-04
# AUTHOR:     MZF
# SCRIPT:     01_summarise_DNMTools_DML_counts_by_partition.R
# VERSION:    02
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Summarise DNMTools RADMeth/RADAdjust output by contrast and genomic partition.
#
#   This script reports:
#
#   1) Percentage of covered CpGs identified as DMLs per contrast and partition.
#   2) Mean ± SD/SEM of DML percentages across contrasts.
#   3) Mean ± SD/SEM of DML percentages by comparison class and partition.
#   4) Chromosome-scale distribution tests for WG and AUTO.
#   5) AUTO-vs-CH12 partition-level enrichment/depletion tests.
#   6) Window-level binomial/quasibinomial GLMs for WG and AUTO.
#   7) Within-CH12 window-level chi-square tests.
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# DETAILS:
#   For each contrast and partition, report:
#
#   1) Covered CpGs:
#        CpG loci with non-zero read coverage in both contrast groups.
#
#   2) DMLs:
#        Covered CpGs with RADAdjust FDR <= 0.01.
#
#   3) Percentage:
#        100 * DMLs / Covered CpGs.
#
#   The script also adds a classification column that separates CpGs into
#   count-based RADMeth status categories analogous to the information obtained
#   with radmeth -na-info:
#
#      - NA_LOW_COV_NO_COVERAGE_BOTH_GROUPS
#      - NA_LOW_COV_ONE_GROUP_ONLY
#      - NA_EXTREME_FULL_CONVERSION_ALL_UNMETHYLATED
#      - NA_EXTREME_NO_CONVERSION_ALL_METHYLATED
#      - NA_OTHER_OR_NUMERICAL
#      - TESTED_OR_NON_NA_PVALUE
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# INPUT:
#   One or more RADAdjust BED files generated from DNMTools, with no header.
#
#   Expected format:
#     Column 1   chr
#     Column 2   pos
#     Column 3   strand
#     Column 4   context
#     Column 5   radmeth_pval
#     Column 6   adjusted_pval
#     Column 7   fdr
#     Column 8   case_total
#     Column 9   case_methylated
#     Column 10  control_total
#     Column 11  control_methylated
#
#   Example input file names:
#     radadjust_f-vs-p.bed
#     radadjust_f-vs-y.bed
#     radadjust_f-vs-i.bed
#     radadjust_p-vs-i.bed
#     radadjust_p-vs-y.bed
#     radadjust_y-vs-i.bed
#
# MAIN CALCULATIONS:
#   - Covered CpG:
#       case_total > 0 AND control_total > 0
#
#   - DML:
#       covered CpG AND fdr <= 0.01
#
#   - Percentage of covered CpGs identified as DMLs:
#       100 * DMLs / Covered CpGs
#
#   - Partition definitions:
#       WG   = all nuclear chromosomes matching NUCLEAR_CHR_PATTERN
#       AUTO = nuclear chromosomes excluding CH12_NAME
#       CH12 = chromosome matching CH12_NAME
#
#   - Chromosome-scale null:
#       DMLs are distributed across chromosomes in proportion to the number of
#       covered CpGs per chromosome.
#
#   - AUTO-vs-CH12 null:
#       DMLs are distributed between AUTO and CH12 in proportion to the number of
#       covered CpGs in each partition.
#
#   - CH12 window null:
#       DMLs are distributed along CH12 windows in proportion to the number of
#       covered CpGs per window.
#
# NOTES:
#   - Fully methylated and fully unmethylated covered CpGs are retained in the
#     denominator because they are covered loci, even when RADMeth returns an NA
#     p-value due to extreme counts.
#   - CpGs with no coverage in one or both groups are excluded from the covered
#     CpG denominator because they are not informative for a two-group DML test.
#   - WG and AUTO can be tested for chromosome-scale heterogeneity.
#   - CH12 cannot be tested for among-chromosome heterogeneity because it contains
#     only one chromosome. For CH12, this script instead provides a window-level
#     goodness-of-fit test along the chromosome.
#   - Modify CH12_NAME and NUCLEAR_CHR_PATTERN below if chromosome names differ.
#
# OUTPUT:
#   01_DNMTools_DML_summary_by_contrast_partition.tsv
#   02_DNMTools_radadjust_with_locus_classification.tsv
#   03_DNMTools_locus_classification_summary.tsv
#   04_DNMTools_mean_DML_percentage_by_partition.tsv
#   05_DNMTools_mean_DML_percentage_by_comparison_class_partition.tsv
#   06_DNMTools_formatted_DML_percentage_report.tsv
#   07_DNMTools_chromosome_level_DML_summary_by_partition.tsv
#   08_DNMTools_chisq_chromosome_distribution_by_partition_contrast.tsv
#   09_DNMTools_chisq_chromosome_residuals_by_partition_contrast.tsv
#   10_DNMTools_AUTO_vs_CH12_partition_chisq_by_contrast.tsv
#   11_DNMTools_AUTO_vs_CH12_partition_residuals_by_contrast.tsv
#   12_DNMTools_window_level_DML_summary_by_partition.tsv
#   13_DNMTools_window_binomial_GLM_chromosome_effect_by_partition_contrast.tsv
#   14_DNMTools_window_quasibinomial_GLM_chromosome_effect_by_partition_contrast.tsv
#   15_DNMTools_CH12_window_chisq_by_contrast.tsv
#   16_DNMTools_CH12_window_chisq_residuals_by_contrast.tsv
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# 0. User-defined settings
# ==============================================================================

input_dir <- "."
output_dir <- "DNMTools_DML_summary_V03"

input_pattern <- "^radadjust_.*\\.bed$"

fdr_cutoff <- 0.01

# Chromosome 12 name in your genome
CH12_NAME <- "Parae_12"

# Nuclear chromosome pattern.
# This keeps chromosomes such as Parae_01, Parae_02, ..., Parae_12, etc.
# It excludes LambdaNEB, mitochondria, pUC19, and other non-nuclear controls.
NUCLEAR_CHR_PATTERN <- "^Parae_[0-9]+$"

# Window settings for window-level models
WINDOW_SIZE_BP <- 100000
MIN_COVERED_CPGS_PER_WINDOW <- 10

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 0A. Helper functions
# ==============================================================================

safe_p_adjust <- function(p_values, method = "BH") {
  out <- rep(NA_real_, length(p_values))
  keep <- !is.na(p_values)
  out[keep] <- p.adjust(p_values[keep], method = method)
  return(out)
}

add_partition_order <- function(dt) {
  dt[, partition_order := fifelse(
    partition == "WG", 1L,
    fifelse(partition == "AUTO", 2L,
            fifelse(partition == "CH12", 3L, 99L))
  )]
  return(dt)
}

add_comparison_class <- function(dt) {
  dt[, comparison_class := fifelse(
    contrast %in% c("f-vs-p", "f-vs-y", "f-vs-i"),
    "sex",
    fifelse(
      contrast %in% c("p-vs-i", "p-vs-y", "y-vs-i"),
      "morph",
      "unknown"
    )
  )]
  return(dt)
}

add_chromosome_number <- function(dt) {
  dt[, chr_number := suppressWarnings(as.integer(sub("^Parae_", "", chr)))]
  return(dt)
}

# ==============================================================================
# 1. Locate input files
# ==============================================================================

radadjust_files <- list.files(
  path = input_dir,
  pattern = input_pattern,
  full.names = TRUE
)

if (length(radadjust_files) == 0) {
  stop(
    "No RADAdjust files found. Expected files matching pattern: ",
    input_pattern,
    " in directory: ",
    input_dir
  )
}

message("Found ", length(radadjust_files), " RADAdjust file(s).")

# ==============================================================================
# 2. Function to read one RADAdjust file
# ==============================================================================

read_radadjust_file <- function(file_path) {
  
  message("Reading: ", basename(file_path))
  
  dt <- fread(
    file_path,
    header = FALSE,
    sep = "\t",
    na.strings = c("NA", "NaN", "nan", ".")
  )
  
  if (ncol(dt) < 11) {
    stop(
      "File has fewer than 11 columns and does not match expected RADAdjust format: ",
      file_path
    )
  }
  
  dt <- dt[, 1:11]
  
  setnames(
    dt,
    c(
      "chr",
      "pos",
      "strand",
      "context",
      "radmeth_pval",
      "adjusted_pval",
      "fdr",
      "case_total",
      "case_methylated",
      "control_total",
      "control_methylated"
    )
  )
  
  contrast_name <- basename(file_path)
  contrast_name <- sub("^radadjust_", "", contrast_name)
  contrast_name <- sub("\\.bed$", "", contrast_name)
  
  dt[, contrast := contrast_name]
  
  numeric_cols <- c(
    "pos",
    "radmeth_pval",
    "adjusted_pval",
    "fdr",
    "case_total",
    "case_methylated",
    "control_total",
    "control_methylated"
  )
  
  for (col in numeric_cols) {
    dt[, (col) := as.numeric(get(col))]
  }
  
  return(dt)
}

# ==============================================================================
# 3. Load all RADAdjust files into one dataframe
# ==============================================================================

dt_all <- rbindlist(
  lapply(radadjust_files, read_radadjust_file),
  use.names = TRUE,
  fill = TRUE
)

message("Loaded ", nrow(dt_all), " total CpG rows.")

# ==============================================================================
# 4. Add locus-level classification columns
# ==============================================================================

dt_all[, covered_both_groups := case_total > 0 & control_total > 0]

dt_all[, fully_unmethylated_both_groups :=
         covered_both_groups &
         case_methylated == 0 &
         control_methylated == 0]

dt_all[, fully_methylated_both_groups :=
         covered_both_groups &
         case_methylated == case_total &
         control_methylated == control_total]

dt_all[, locus_class := fifelse(
  case_total == 0 & control_total == 0,
  "NA_LOW_COV_NO_COVERAGE_BOTH_GROUPS",
  fifelse(
    case_total == 0 | control_total == 0,
    "NA_LOW_COV_ONE_GROUP_ONLY",
    fifelse(
      fully_unmethylated_both_groups,
      "NA_EXTREME_FULL_CONVERSION_ALL_UNMETHYLATED",
      fifelse(
        fully_methylated_both_groups,
        "NA_EXTREME_NO_CONVERSION_ALL_METHYLATED",
        fifelse(
          is.na(radmeth_pval),
          "NA_OTHER_OR_NUMERICAL",
          "TESTED_OR_NON_NA_PVALUE"
        )
      )
    )
  )
)]

dt_all[, is_DML := covered_both_groups & !is.na(fdr) & fdr <= fdr_cutoff]

# ==============================================================================
# 5. Define chromosome classes
# ==============================================================================

dt_all[, is_nuclear := grepl(NUCLEAR_CHR_PATTERN, chr)]

dt_all[, chromosome_class := fifelse(
  chr == CH12_NAME,
  "CH12",
  fifelse(
    is_nuclear & chr != CH12_NAME,
    "AUTO",
    "NON_NUCLEAR_OR_CONTROL"
  )
)]

# ==============================================================================
# 6. Create partition-expanded table
# ==============================================================================

dt_wg <- dt_all[is_nuclear == TRUE]
dt_wg[, partition := "WG"]

dt_auto <- dt_all[chromosome_class == "AUTO"]
dt_auto[, partition := "AUTO"]

dt_ch12 <- dt_all[chromosome_class == "CH12"]
dt_ch12[, partition := "CH12"]

dt_partitioned <- rbindlist(
  list(dt_wg, dt_auto, dt_ch12),
  use.names = TRUE,
  fill = TRUE
)

if (nrow(dt_partitioned) == 0) {
  stop(
    "No rows were assigned to WG, AUTO, or CH12. Check CH12_NAME and NUCLEAR_CHR_PATTERN."
  )
}

# ==============================================================================
# 7. Summarise covered CpGs, DMLs, and percentage per contrast and partition
# ==============================================================================

summary_dt <- dt_partitioned[
  ,
  .(
    total_CpG_rows_in_partition = .N,
    covered_CpGs = sum(covered_both_groups, na.rm = TRUE),
    DMLs_FDR_le_0_01 = sum(is_DML, na.rm = TRUE)
  ),
  by = .(contrast, partition)
]

summary_dt[, percentage_covered_CpGs_identified_as_DMLs :=
             fifelse(
               covered_CpGs > 0,
               100 * DMLs_FDR_le_0_01 / covered_CpGs,
               NA_real_
             )]

summary_dt <- add_partition_order(summary_dt)
setorder(summary_dt, contrast, partition_order)
summary_dt[, partition_order := NULL]

# ==============================================================================
# 8. Classification summary
# ==============================================================================

classification_summary <- dt_all[
  ,
  .N,
  by = .(contrast, chromosome_class, locus_class)
]

setorder(classification_summary, contrast, chromosome_class, locus_class)

# ==============================================================================
# 8A. Average DML percentage across contrasts per partition
# ==============================================================================

partition_percentage_summary <- summary_dt[
  !is.na(percentage_covered_CpGs_identified_as_DMLs),
  .(
    n_contrasts = .N,
    mean_percentage_DMLs = mean(percentage_covered_CpGs_identified_as_DMLs, na.rm = TRUE),
    sd_percentage_DMLs = sd(percentage_covered_CpGs_identified_as_DMLs, na.rm = TRUE),
    sem_percentage_DMLs = sd(percentage_covered_CpGs_identified_as_DMLs, na.rm = TRUE) / sqrt(.N)
  ),
  by = partition
]

partition_percentage_summary <- add_partition_order(partition_percentage_summary)
setorder(partition_percentage_summary, partition_order)
partition_percentage_summary[, partition_order := NULL]

# ==============================================================================
# 8B. Mean percentage of DMLs by comparison class and partition
# ==============================================================================

summary_dt <- add_comparison_class(summary_dt)

if (any(summary_dt$comparison_class == "unknown")) {
  warning(
    "Some contrasts were classified as 'unknown'. Check contrast names: ",
    paste(unique(summary_dt[comparison_class == "unknown", contrast]), collapse = ", ")
  )
}

comparison_partition_summary <- summary_dt[
  !is.na(percentage_covered_CpGs_identified_as_DMLs) &
    comparison_class != "unknown",
  .(
    n_contrasts = .N,
    mean_percentage_DMLs = mean(
      percentage_covered_CpGs_identified_as_DMLs,
      na.rm = TRUE
    ),
    sd_percentage_DMLs = sd(
      percentage_covered_CpGs_identified_as_DMLs,
      na.rm = TRUE
    ),
    sem_percentage_DMLs = sd(
      percentage_covered_CpGs_identified_as_DMLs,
      na.rm = TRUE
    ) / sqrt(.N)
  ),
  by = .(comparison_class, partition)
]

comparison_partition_summary <- add_partition_order(comparison_partition_summary)

comparison_partition_summary[, class_order := fifelse(
  comparison_class == "sex", 1L,
  fifelse(comparison_class == "morph", 2L, 99L)
)]

setorder(comparison_partition_summary, class_order, partition_order)

comparison_partition_summary[, c("class_order", "partition_order") := NULL]

comparison_partition_report <- copy(comparison_partition_summary)

comparison_partition_report[, mean_SD_SEM := sprintf(
  "%.4f ± %.4f / %.4f",
  mean_percentage_DMLs,
  sd_percentage_DMLs,
  sem_percentage_DMLs
)]

comparison_partition_report <- comparison_partition_report[
  ,
  .(
    comparison_class,
    partition,
    n_contrasts,
    `mean % DMLs ± SD / SEM` = mean_SD_SEM
  )
]

# ==============================================================================
# 9. Save core summary outputs
# ==============================================================================

main_summary_file <- file.path(
  output_dir,
  "01_DNMTools_DML_summary_by_contrast_partition.tsv"
)

classified_table_file <- file.path(
  output_dir,
  "02_DNMTools_radadjust_with_locus_classification.tsv"
)

classification_summary_file <- file.path(
  output_dir,
  "03_DNMTools_locus_classification_summary.tsv"
)

partition_percentage_summary_file <- file.path(
  output_dir,
  "04_DNMTools_mean_DML_percentage_by_partition.tsv"
)

comparison_partition_summary_file <- file.path(
  output_dir,
  "05_DNMTools_mean_DML_percentage_by_comparison_class_partition.tsv"
)

comparison_partition_report_file <- file.path(
  output_dir,
  "06_DNMTools_formatted_DML_percentage_report.tsv"
)

fwrite(summary_dt, main_summary_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(dt_all, classified_table_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(classification_summary, classification_summary_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(partition_percentage_summary, partition_percentage_summary_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(comparison_partition_summary, comparison_partition_summary_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(comparison_partition_report, comparison_partition_report_file, sep = "\t", quote = FALSE, na = "NA")

message("Saved main summary to: ", main_summary_file)
message("Saved classified RADAdjust table to: ", classified_table_file)
message("Saved classification summary to: ", classification_summary_file)
message("Saved mean DML percentage summary to: ", partition_percentage_summary_file)
message("Saved numeric comparison-class summary to: ", comparison_partition_summary_file)
message("Saved formatted comparison-class report to: ", comparison_partition_report_file)

# ==============================================================================
# 10. Chromosome-level DML summary by contrast and partition
# ==============================================================================

chromosome_summary_by_partition <- dt_partitioned[
  ,
  .(
    covered_CpGs = sum(covered_both_groups, na.rm = TRUE),
    DMLs_FDR_le_0_01 = sum(is_DML, na.rm = TRUE)
  ),
  by = .(contrast, partition, chr)
]

chromosome_summary_by_partition[, non_DMLs := covered_CpGs - DMLs_FDR_le_0_01]

chromosome_summary_by_partition[, percentage_DMLs := fifelse(
  covered_CpGs > 0,
  100 * DMLs_FDR_le_0_01 / covered_CpGs,
  NA_real_
)]

chromosome_summary_by_partition <- chromosome_summary_by_partition[
  covered_CpGs > 0
]

chromosome_summary_by_partition <- add_partition_order(chromosome_summary_by_partition)
chromosome_summary_by_partition <- add_chromosome_number(chromosome_summary_by_partition)

setorder(
  chromosome_summary_by_partition,
  contrast,
  partition_order,
  chr_number
)

chromosome_summary_by_partition[, c("partition_order", "chr_number") := NULL]

chromosome_summary_file <- file.path(
  output_dir,
  "07_DNMTools_chromosome_level_DML_summary_by_partition.tsv"
)

fwrite(
  chromosome_summary_by_partition,
  chromosome_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Saved chromosome-level DML summary to: ", chromosome_summary_file)

# ==============================================================================
# 11. Chi-square goodness-of-fit across chromosomes by contrast and partition
# ==============================================================================

run_chisq_gof_chromosomes <- function(dt_subset) {
  
  contrast_id <- unique(dt_subset$contrast)
  partition_id <- unique(dt_subset$partition)
  
  if (length(contrast_id) != 1 || length(partition_id) != 1) {
    stop("run_chisq_gof_chromosomes() expects one contrast and one partition.")
  }
  
  total_covered <- sum(dt_subset$covered_CpGs, na.rm = TRUE)
  total_DMLs <- sum(dt_subset$DMLs_FDR_le_0_01, na.rm = TRUE)
  n_chromosomes <- length(unique(dt_subset$chr))
  
  if (n_chromosomes < 2) {
    
    global_result <- data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = total_covered,
      total_DMLs = total_DMLs,
      chisq_statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      min_expected_DMLs = NA_real_,
      n_chromosomes_expected_lt_5 = NA_integer_,
      test_status = "not_tested_single_chromosome"
    )
    
    residual_result <- copy(dt_subset)
    residual_result[, expected_DMLs := NA_real_]
    residual_result[, pearson_residual := NA_real_]
    residual_result[, standardized_residual := NA_real_]
    residual_result[, residual_direction := "not_tested_single_chromosome"]
    
    return(list(global = global_result, residuals = residual_result))
  }
  
  if (total_covered == 0 || total_DMLs == 0) {
    
    global_result <- data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = total_covered,
      total_DMLs = total_DMLs,
      chisq_statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      min_expected_DMLs = NA_real_,
      n_chromosomes_expected_lt_5 = NA_integer_,
      test_status = "not_tested_insufficient_DMLs_or_coverage"
    )
    
    residual_result <- copy(dt_subset)
    residual_result[, expected_DMLs := NA_real_]
    residual_result[, pearson_residual := NA_real_]
    residual_result[, standardized_residual := NA_real_]
    residual_result[, residual_direction := "not_tested_insufficient_DMLs_or_coverage"]
    
    return(list(global = global_result, residuals = residual_result))
  }
  
  expected_prob <- dt_subset$covered_CpGs / total_covered
  
  chisq_obj <- suppressWarnings(
    chisq.test(
      x = dt_subset$DMLs_FDR_le_0_01,
      p = expected_prob,
      rescale.p = TRUE
    )
  )
  
  expected_counts <- as.numeric(chisq_obj$expected)
  
  global_result <- data.table(
    contrast = contrast_id,
    partition = partition_id,
    n_chromosomes = n_chromosomes,
    total_covered_CpGs = total_covered,
    total_DMLs = total_DMLs,
    chisq_statistic = unname(chisq_obj$statistic),
    df = unname(chisq_obj$parameter),
    p_value = chisq_obj$p.value,
    p_adjust_BH = NA_real_,
    min_expected_DMLs = min(expected_counts, na.rm = TRUE),
    n_chromosomes_expected_lt_5 = sum(expected_counts < 5, na.rm = TRUE),
    test_status = fifelse(
      any(expected_counts < 5),
      "tested_but_some_expected_counts_lt_5",
      "tested"
    )
  )
  
  residual_result <- copy(dt_subset)
  
  residual_result[, expected_DMLs := expected_counts]
  
  residual_result[, pearson_residual :=
                    (DMLs_FDR_le_0_01 - expected_DMLs) / sqrt(expected_DMLs)]
  
  residual_result[, standardized_residual := as.numeric(chisq_obj$stdres)]
  
  residual_result[, residual_direction := fifelse(
    standardized_residual > 0,
    "DML_excess",
    fifelse(
      standardized_residual < 0,
      "DML_depletion",
      "as_expected"
    )
  )]
  
  return(list(global = global_result, residuals = residual_result))
}

chisq_chromosome_results_list <- lapply(
  split(
    chromosome_summary_by_partition,
    by = c("contrast", "partition"),
    keep.by = TRUE
  ),
  run_chisq_gof_chromosomes
)

chisq_chromosome_global <- rbindlist(
  lapply(chisq_chromosome_results_list, function(x) x$global),
  use.names = TRUE,
  fill = TRUE
)

chisq_chromosome_residuals <- rbindlist(
  lapply(chisq_chromosome_results_list, function(x) x$residuals),
  use.names = TRUE,
  fill = TRUE
)

chisq_chromosome_global[, p_adjust_BH := safe_p_adjust(p_value, method = "BH")]

chisq_chromosome_global <- add_partition_order(chisq_chromosome_global)
setorder(chisq_chromosome_global, contrast, partition_order)
chisq_chromosome_global[, partition_order := NULL]

chisq_chromosome_residuals <- add_partition_order(chisq_chromosome_residuals)
chisq_chromosome_residuals <- add_chromosome_number(chisq_chromosome_residuals)
setorder(chisq_chromosome_residuals, contrast, partition_order, chr_number)
chisq_chromosome_residuals[, c("partition_order", "chr_number") := NULL]

chisq_chromosome_global_file <- file.path(
  output_dir,
  "08_DNMTools_chisq_chromosome_distribution_by_partition_contrast.tsv"
)

chisq_chromosome_residuals_file <- file.path(
  output_dir,
  "09_DNMTools_chisq_chromosome_residuals_by_partition_contrast.tsv"
)

fwrite(
  chisq_chromosome_global,
  chisq_chromosome_global_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  chisq_chromosome_residuals,
  chisq_chromosome_residuals_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Saved chromosome chi-square global results to: ", chisq_chromosome_global_file)
message("Saved chromosome chi-square residuals to: ", chisq_chromosome_residuals_file)

# ==============================================================================
# 12. AUTO-vs-CH12 partition-level chi-square goodness-of-fit
# ==============================================================================

run_chisq_gof_partitions <- function(dt_subset) {
  
  contrast_id <- unique(dt_subset$contrast)
  
  if (length(contrast_id) != 1) {
    stop("run_chisq_gof_partitions() expects one contrast at a time.")
  }
  
  dt_subset <- dt_subset[partition %in% c("AUTO", "CH12")]
  dt_subset <- dt_subset[covered_CpGs > 0]
  
  total_covered <- sum(dt_subset$covered_CpGs, na.rm = TRUE)
  total_DMLs <- sum(dt_subset$DMLs_FDR_le_0_01, na.rm = TRUE)
  n_partitions <- length(unique(dt_subset$partition))
  
  if (n_partitions < 2 || total_covered == 0 || total_DMLs == 0) {
    
    global_result <- data.table(
      contrast = contrast_id,
      n_partitions = n_partitions,
      total_covered_CpGs = total_covered,
      total_DMLs = total_DMLs,
      chisq_statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      min_expected_DMLs = NA_real_,
      n_partitions_expected_lt_5 = NA_integer_,
      test_status = "not_tested_insufficient_partitions_DMLs_or_coverage"
    )
    
    residual_result <- copy(dt_subset)
    residual_result[, expected_DMLs := NA_real_]
    residual_result[, pearson_residual := NA_real_]
    residual_result[, standardized_residual := NA_real_]
    residual_result[, residual_direction := "not_tested"]
    
    return(list(global = global_result, residuals = residual_result))
  }
  
  expected_prob <- dt_subset$covered_CpGs / total_covered
  
  chisq_obj <- suppressWarnings(
    chisq.test(
      x = dt_subset$DMLs_FDR_le_0_01,
      p = expected_prob,
      rescale.p = TRUE
    )
  )
  
  expected_counts <- as.numeric(chisq_obj$expected)
  
  global_result <- data.table(
    contrast = contrast_id,
    n_partitions = n_partitions,
    total_covered_CpGs = total_covered,
    total_DMLs = total_DMLs,
    chisq_statistic = unname(chisq_obj$statistic),
    df = unname(chisq_obj$parameter),
    p_value = chisq_obj$p.value,
    p_adjust_BH = NA_real_,
    min_expected_DMLs = min(expected_counts, na.rm = TRUE),
    n_partitions_expected_lt_5 = sum(expected_counts < 5, na.rm = TRUE),
    test_status = fifelse(
      any(expected_counts < 5),
      "tested_but_some_expected_counts_lt_5",
      "tested"
    )
  )
  
  residual_result <- copy(dt_subset)
  
  residual_result[, expected_DMLs := expected_counts]
  
  residual_result[, pearson_residual :=
                    (DMLs_FDR_le_0_01 - expected_DMLs) / sqrt(expected_DMLs)]
  
  residual_result[, standardized_residual := as.numeric(chisq_obj$stdres)]
  
  residual_result[, residual_direction := fifelse(
    standardized_residual > 0,
    "DML_excess",
    fifelse(
      standardized_residual < 0,
      "DML_depletion",
      "as_expected"
    )
  )]
  
  return(list(global = global_result, residuals = residual_result))
}

partition_chisq_results_list <- lapply(
  split(summary_dt[partition %in% c("AUTO", "CH12")], by = "contrast", keep.by = TRUE),
  run_chisq_gof_partitions
)

partition_chisq_global <- rbindlist(
  lapply(partition_chisq_results_list, function(x) x$global),
  use.names = TRUE,
  fill = TRUE
)

partition_chisq_residuals <- rbindlist(
  lapply(partition_chisq_results_list, function(x) x$residuals),
  use.names = TRUE,
  fill = TRUE
)

partition_chisq_global[, p_adjust_BH := safe_p_adjust(p_value, method = "BH")]

setorder(partition_chisq_global, contrast)

partition_chisq_residuals <- add_partition_order(partition_chisq_residuals)
setorder(partition_chisq_residuals, contrast, partition_order)
partition_chisq_residuals[, partition_order := NULL]

partition_chisq_global_file <- file.path(
  output_dir,
  "10_DNMTools_AUTO_vs_CH12_partition_chisq_by_contrast.tsv"
)

partition_chisq_residuals_file <- file.path(
  output_dir,
  "11_DNMTools_AUTO_vs_CH12_partition_residuals_by_contrast.tsv"
)

fwrite(
  partition_chisq_global,
  partition_chisq_global_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  partition_chisq_residuals,
  partition_chisq_residuals_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Saved AUTO-vs-CH12 partition chi-square results to: ", partition_chisq_global_file)
message("Saved AUTO-vs-CH12 partition residuals to: ", partition_chisq_residuals_file)

# ==============================================================================
# 13. Generate window-level DML summary by contrast and partition
# ==============================================================================

window_summary <- dt_partitioned[
  covered_both_groups == TRUE,
  .(
    covered_CpGs = .N,
    DMLs_FDR_le_0_01 = sum(is_DML, na.rm = TRUE)
  ),
  by = .(
    contrast,
    partition,
    chr,
    window_start = floor(pos / WINDOW_SIZE_BP) * WINDOW_SIZE_BP
  )
]

window_summary[, window_end := window_start + WINDOW_SIZE_BP]

window_summary[, non_DMLs := covered_CpGs - DMLs_FDR_le_0_01]

window_summary[, percentage_DMLs := fifelse(
  covered_CpGs > 0,
  100 * DMLs_FDR_le_0_01 / covered_CpGs,
  NA_real_
)]

window_summary <- window_summary[
  covered_CpGs >= MIN_COVERED_CPGS_PER_WINDOW
]

window_summary <- add_partition_order(window_summary)
window_summary <- add_chromosome_number(window_summary)

setorder(
  window_summary,
  contrast,
  partition_order,
  chr_number,
  window_start
)

window_summary[, c("partition_order", "chr_number") := NULL]

window_summary_file <- file.path(
  output_dir,
  "12_DNMTools_window_level_DML_summary_by_partition.tsv"
)

fwrite(
  window_summary,
  window_summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Saved window-level DML summary to: ", window_summary_file)

# ==============================================================================
# 14. Window-level binomial GLM: chromosome effect within each partition
# ==============================================================================

run_window_binomial_glm <- function(dt_subset) {
  
  contrast_id <- unique(dt_subset$contrast)
  partition_id <- unique(dt_subset$partition)
  
  if (length(contrast_id) != 1 || length(partition_id) != 1) {
    stop("run_window_binomial_glm() expects one contrast and one partition.")
  }
  
  dt_model <- copy(dt_subset)
  
  dt_model <- dt_model[
    covered_CpGs > 0 &
      !is.na(DMLs_FDR_le_0_01) &
      !is.na(non_DMLs)
  ]
  
  total_DMLs <- sum(dt_model$DMLs_FDR_le_0_01, na.rm = TRUE)
  total_non_DMLs <- sum(dt_model$non_DMLs, na.rm = TRUE)
  n_chromosomes <- length(unique(dt_model$chr))
  
  if (n_chromosomes < 2) {
    
    return(data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_windows = nrow(dt_model),
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
      total_DMLs = total_DMLs,
      glm_chisq_statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      overdispersion_ratio = NA_real_,
      test_status = "not_tested_single_chromosome"
    ))
  }
  
  if (
    nrow(dt_model) < 2 ||
    total_DMLs == 0 ||
    total_non_DMLs == 0
  ) {
    
    return(data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_windows = nrow(dt_model),
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
      total_DMLs = total_DMLs,
      glm_chisq_statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      overdispersion_ratio = NA_real_,
      test_status = "not_tested_insufficient_variation"
    ))
  }
  
  dt_model[, chr := factor(chr)]
  
  result <- tryCatch({
    
    null_model <- glm(
      cbind(DMLs_FDR_le_0_01, non_DMLs) ~ 1,
      family = binomial,
      data = dt_model
    )
    
    full_model <- glm(
      cbind(DMLs_FDR_le_0_01, non_DMLs) ~ chr,
      family = binomial,
      data = dt_model
    )
    
    if (is.na(full_model$df.residual) || full_model$df.residual <= 0) {
      return(data.table(
        contrast = contrast_id,
        partition = partition_id,
        n_windows = nrow(dt_model),
        n_chromosomes = n_chromosomes,
        total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
        total_DMLs = total_DMLs,
        glm_chisq_statistic = NA_real_,
        df = NA_real_,
        p_value = NA_real_,
        p_adjust_BH = NA_real_,
        overdispersion_ratio = NA_real_,
        test_status = "not_tested_saturated_model"
      ))
    }
    
    glm_test <- anova(
      null_model,
      full_model,
      test = "Chisq"
    )
    
    overdispersion_ratio <- sum(
      residuals(full_model, type = "pearson")^2,
      na.rm = TRUE
    ) / full_model$df.residual
    
    data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_windows = nrow(dt_model),
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
      total_DMLs = total_DMLs,
      glm_chisq_statistic = glm_test$Deviance[2],
      df = glm_test$Df[2],
      p_value = glm_test$`Pr(>Chi)`[2],
      p_adjust_BH = NA_real_,
      overdispersion_ratio = overdispersion_ratio,
      test_status = fifelse(
        overdispersion_ratio > 2,
        "tested_overdispersion_detected_use_quasibinomial",
        "tested"
      )
    )
    
  }, error = function(e) {
    
    data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_windows = nrow(dt_model),
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
      total_DMLs = total_DMLs,
      glm_chisq_statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      overdispersion_ratio = NA_real_,
      test_status = paste0("error: ", conditionMessage(e))
    )
  })
  
  return(result)
}

window_binomial_glm_results <- rbindlist(
  lapply(
    split(window_summary, by = c("contrast", "partition"), keep.by = TRUE),
    run_window_binomial_glm
  ),
  use.names = TRUE,
  fill = TRUE
)

window_binomial_glm_results[, p_adjust_BH := safe_p_adjust(p_value, method = "BH")]

window_binomial_glm_results <- add_partition_order(window_binomial_glm_results)
setorder(window_binomial_glm_results, contrast, partition_order)
window_binomial_glm_results[, partition_order := NULL]

window_binomial_glm_file <- file.path(
  output_dir,
  "13_DNMTools_window_binomial_GLM_chromosome_effect_by_partition_contrast.tsv"
)

fwrite(
  window_binomial_glm_results,
  window_binomial_glm_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Saved window-level binomial GLM results to: ", window_binomial_glm_file)

# ==============================================================================
# 15. Window-level quasibinomial GLM: chromosome effect within each partition
# ==============================================================================

run_window_quasibinomial_glm <- function(dt_subset) {
  
  contrast_id <- unique(dt_subset$contrast)
  partition_id <- unique(dt_subset$partition)
  
  if (length(contrast_id) != 1 || length(partition_id) != 1) {
    stop("run_window_quasibinomial_glm() expects one contrast and one partition.")
  }
  
  dt_model <- copy(dt_subset)
  
  dt_model <- dt_model[
    covered_CpGs > 0 &
      !is.na(DMLs_FDR_le_0_01) &
      !is.na(non_DMLs)
  ]
  
  total_DMLs <- sum(dt_model$DMLs_FDR_le_0_01, na.rm = TRUE)
  total_non_DMLs <- sum(dt_model$non_DMLs, na.rm = TRUE)
  n_chromosomes <- length(unique(dt_model$chr))
  
  if (n_chromosomes < 2) {
    
    return(data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_windows = nrow(dt_model),
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
      total_DMLs = total_DMLs,
      quasi_F_statistic = NA_real_,
      df_numerator = NA_real_,
      df_denominator = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      test_status = "not_tested_single_chromosome"
    ))
  }
  
  if (
    nrow(dt_model) < 2 ||
    total_DMLs == 0 ||
    total_non_DMLs == 0
  ) {
    
    return(data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_windows = nrow(dt_model),
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
      total_DMLs = total_DMLs,
      quasi_F_statistic = NA_real_,
      df_numerator = NA_real_,
      df_denominator = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      test_status = "not_tested_insufficient_variation"
    ))
  }
  
  dt_model[, chr := factor(chr)]
  
  result <- tryCatch({
    
    null_model <- glm(
      cbind(DMLs_FDR_le_0_01, non_DMLs) ~ 1,
      family = quasibinomial,
      data = dt_model
    )
    
    full_model <- glm(
      cbind(DMLs_FDR_le_0_01, non_DMLs) ~ chr,
      family = quasibinomial,
      data = dt_model
    )
    
    if (is.na(full_model$df.residual) || full_model$df.residual <= 0) {
      return(data.table(
        contrast = contrast_id,
        partition = partition_id,
        n_windows = nrow(dt_model),
        n_chromosomes = n_chromosomes,
        total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
        total_DMLs = total_DMLs,
        quasi_F_statistic = NA_real_,
        df_numerator = NA_real_,
        df_denominator = full_model$df.residual,
        p_value = NA_real_,
        p_adjust_BH = NA_real_,
        test_status = "not_tested_saturated_model"
      ))
    }
    
    quasi_test <- anova(
      null_model,
      full_model,
      test = "F"
    )
    
    data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_windows = nrow(dt_model),
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
      total_DMLs = total_DMLs,
      quasi_F_statistic = quasi_test$F[2],
      df_numerator = quasi_test$Df[2],
      df_denominator = full_model$df.residual,
      p_value = quasi_test$`Pr(>F)`[2],
      p_adjust_BH = NA_real_,
      test_status = "tested"
    )
    
  }, error = function(e) {
    
    data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_windows = nrow(dt_model),
      n_chromosomes = n_chromosomes,
      total_covered_CpGs = sum(dt_model$covered_CpGs, na.rm = TRUE),
      total_DMLs = total_DMLs,
      quasi_F_statistic = NA_real_,
      df_numerator = NA_real_,
      df_denominator = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      test_status = paste0("error: ", conditionMessage(e))
    )
  })
  
  return(result)
}

window_quasibinomial_glm_results <- rbindlist(
  lapply(
    split(window_summary, by = c("contrast", "partition"), keep.by = TRUE),
    run_window_quasibinomial_glm
  ),
  use.names = TRUE,
  fill = TRUE
)

window_quasibinomial_glm_results[, p_adjust_BH := safe_p_adjust(p_value, method = "BH")]

window_quasibinomial_glm_results <- add_partition_order(window_quasibinomial_glm_results)
setorder(window_quasibinomial_glm_results, contrast, partition_order)
window_quasibinomial_glm_results[, partition_order := NULL]

window_quasibinomial_glm_file <- file.path(
  output_dir,
  "14_DNMTools_window_quasibinomial_GLM_chromosome_effect_by_partition_contrast.tsv"
)

fwrite(
  window_quasibinomial_glm_results,
  window_quasibinomial_glm_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Saved window-level quasibinomial GLM results to: ", window_quasibinomial_glm_file)

# ==============================================================================
# 16. Within-CH12 window-level chi-square goodness-of-fit
# ==============================================================================

run_chisq_gof_windows <- function(dt_subset) {
  
  contrast_id <- unique(dt_subset$contrast)
  partition_id <- unique(dt_subset$partition)
  
  if (length(contrast_id) != 1 || length(partition_id) != 1) {
    stop("run_chisq_gof_windows() expects one contrast and one partition.")
  }
  
  total_covered <- sum(dt_subset$covered_CpGs, na.rm = TRUE)
  total_DMLs <- sum(dt_subset$DMLs_FDR_le_0_01, na.rm = TRUE)
  n_windows <- nrow(dt_subset)
  
  if (n_windows < 2 || total_covered == 0 || total_DMLs == 0) {
    
    global_result <- data.table(
      contrast = contrast_id,
      partition = partition_id,
      n_windows = n_windows,
      total_covered_CpGs = total_covered,
      total_DMLs = total_DMLs,
      chisq_statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_,
      p_adjust_BH = NA_real_,
      min_expected_DMLs = NA_real_,
      n_windows_expected_lt_5 = NA_integer_,
      test_status = "not_tested_insufficient_windows_DMLs_or_coverage"
    )
    
    residual_result <- copy(dt_subset)
    residual_result[, expected_DMLs := NA_real_]
    residual_result[, pearson_residual := NA_real_]
    residual_result[, standardized_residual := NA_real_]
    residual_result[, residual_direction := "not_tested"]
    
    return(list(global = global_result, residuals = residual_result))
  }
  
  expected_prob <- dt_subset$covered_CpGs / total_covered
  
  chisq_obj <- suppressWarnings(
    chisq.test(
      x = dt_subset$DMLs_FDR_le_0_01,
      p = expected_prob,
      rescale.p = TRUE
    )
  )
  
  expected_counts <- as.numeric(chisq_obj$expected)
  
  global_result <- data.table(
    contrast = contrast_id,
    partition = partition_id,
    n_windows = n_windows,
    total_covered_CpGs = total_covered,
    total_DMLs = total_DMLs,
    chisq_statistic = unname(chisq_obj$statistic),
    df = unname(chisq_obj$parameter),
    p_value = chisq_obj$p.value,
    p_adjust_BH = NA_real_,
    min_expected_DMLs = min(expected_counts, na.rm = TRUE),
    n_windows_expected_lt_5 = sum(expected_counts < 5, na.rm = TRUE),
    test_status = fifelse(
      any(expected_counts < 5),
      "tested_but_some_expected_counts_lt_5",
      "tested"
    )
  )
  
  residual_result <- copy(dt_subset)
  
  residual_result[, expected_DMLs := expected_counts]
  
  residual_result[, pearson_residual :=
                    (DMLs_FDR_le_0_01 - expected_DMLs) / sqrt(expected_DMLs)]
  
  residual_result[, standardized_residual := as.numeric(chisq_obj$stdres)]
  
  residual_result[, residual_direction := fifelse(
    standardized_residual > 0,
    "DML_excess",
    fifelse(
      standardized_residual < 0,
      "DML_depletion",
      "as_expected"
    )
  )]
  
  return(list(global = global_result, residuals = residual_result))
}

ch12_window_summary <- window_summary[partition == "CH12"]

if (nrow(ch12_window_summary) > 0) {
  
  ch12_window_chisq_results_list <- lapply(
    split(ch12_window_summary, by = "contrast", keep.by = TRUE),
    run_chisq_gof_windows
  )
  
  ch12_window_chisq_global <- rbindlist(
    lapply(ch12_window_chisq_results_list, function(x) x$global),
    use.names = TRUE,
    fill = TRUE
  )
  
  ch12_window_chisq_residuals <- rbindlist(
    lapply(ch12_window_chisq_results_list, function(x) x$residuals),
    use.names = TRUE,
    fill = TRUE
  )
  
  ch12_window_chisq_global[, p_adjust_BH := safe_p_adjust(p_value, method = "BH")]
  
  setorder(ch12_window_chisq_global, contrast)
  setorder(ch12_window_chisq_residuals, contrast, window_start)
  
} else {
  
  ch12_window_chisq_global <- data.table(
    contrast = character(),
    partition = character(),
    n_windows = integer(),
    total_covered_CpGs = numeric(),
    total_DMLs = numeric(),
    chisq_statistic = numeric(),
    df = numeric(),
    p_value = numeric(),
    p_adjust_BH = numeric(),
    min_expected_DMLs = numeric(),
    n_windows_expected_lt_5 = integer(),
    test_status = character()
  )
  
  ch12_window_chisq_residuals <- data.table()
}

ch12_window_chisq_global_file <- file.path(
  output_dir,
  "15_DNMTools_CH12_window_chisq_by_contrast.tsv"
)

ch12_window_chisq_residuals_file <- file.path(
  output_dir,
  "16_DNMTools_CH12_window_chisq_residuals_by_contrast.tsv"
)

fwrite(
  ch12_window_chisq_global,
  ch12_window_chisq_global_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  ch12_window_chisq_residuals,
  ch12_window_chisq_residuals_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("Saved CH12 window chi-square results to: ", ch12_window_chisq_global_file)
message("Saved CH12 window chi-square residuals to: ", ch12_window_chisq_residuals_file)

# ==============================================================================
# 17. Print key summaries to terminal
# ==============================================================================

message("\nMain DML summary:")
print(summary_dt)

message("\nMean DML percentage by partition:")
print(partition_percentage_summary)

message("\nMean DML percentage by comparison class and partition:")
print(comparison_partition_summary)

message("\nChromosome-level chi-square results:")
print(chisq_chromosome_global)

message("\nAUTO-vs-CH12 partition chi-square results:")
print(partition_chisq_global)

message("\nWindow-level binomial GLM results:")
print(window_binomial_glm_results)

message("\nWindow-level quasibinomial GLM results:")
print(window_quasibinomial_glm_results)

message("\nCH12 window chi-square results:")
print(ch12_window_chisq_global)

message("\nAll tasks completed.")
