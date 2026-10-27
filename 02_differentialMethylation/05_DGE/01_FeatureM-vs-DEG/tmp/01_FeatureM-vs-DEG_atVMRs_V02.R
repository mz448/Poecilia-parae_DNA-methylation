#!/usr/bin/env Rscript

# DATE: 20260809
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Test whether the magnitude of methylation differentiation at specific
#   genetic features is linearly associated with the magnitude of differential
#   gene expression.
#
#   Analysis:
#
#       X = absolute difference in morph-level methylation
#           |meth.morph.A - meth.morph.B|
#
#       Y = absolute differential expression effect size
#           |log2FoldChange|
#
#   The analysis is performed independently for:
#       - u2000
#       - geneBody
#       - d2000
#
#   and for the six biological contrasts:
#       - i_vs_y_01
#       - p_vs_i_01
#       - y_vs_p_01
#       - f_vs_i
#       - f_vs_p
#       - f_vs_y
#
#   Each feature produces one 6-panel scatterplot.
#
#   Linear model:
#       abs(log2FC) ~ abs(delta methylation)
#
#   Plot annotation:
#       R²
#       P
#       n
#
#   Regression P-values are additionally corrected across all 18
#   feature × contrast tests using Benjamini-Hochberg FDR.
#
# DETAILS:
#   - DEG significance is NOT used as a filter.
#   - Coverage is NOT used as a filter.
#   - Missing data are removed independently within each contrast.
#   - Morph-level methylation values are used.
#   - gene + feature must be unique in the methylation table.
#   - gene_id must be unique in the DEG table.
#   - Duplicate rows are exported and cause the script to stop rather than
#     silently averaging or duplicating genes.
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
#
# NOTE ON GENE × FEATURE COLLAPSING:
#   The methylation input can contain multiple genomic intervals assigned to the
#   same gene and feature (u2000, geneBody, or d2000), reflecting multiple
#   annotated instances of that feature, for example from alternative transcript
#   isoforms.
#
#   Because the differential-expression table contains one log2FoldChange value
#   per gene, these multiple feature intervals are collapsed to a single
#   gene-level methylation value before comparison with gene expression.
#
#   For each gene × feature × morph combination, the arithmetic mean of the
#   available morph-level methylation values across all associated intervals is
#   calculated. Missing interval values are ignored when at least one valid value
#   is available; if all intervals are missing for a morph, the collapsed value
#   remains NA.
#
#   Coverage is not used for weighting. Therefore, each annotated interval
#   contributes equally to the gene-level feature mean.
#
#   Contrast-specific methylation differences are calculated only AFTER this
#   collapsing step:
#
#       abs_delta_meth = | mean_meth_morph_A - mean_meth_morph_B |
#
#   Thus, the analysis evaluates the difference between the mean methylation
#   states of a gene-associated feature, rather than the mean of interval-specific
#   methylation differences.
#
#   The number of intervals contributing to each gene × feature combination is
#   retained as n_intervals for reporting and quality control.
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# CHROMOSOME-SPECIFIC ANALYSIS:
#   The gene-level analysis is performed separately for autosomal genes and
#   genes located on chromosome 12 (Parae_12).
#
#   Multiple intervals are first collapsed to one mean methylation value per
#   gene × feature, as described above. Genes are then assigned to one of two
#   genomic partitions:
#
#       Autosome       = features annotated as autosomal
#       Chromosome_12  = features located on Parae_12
#
#   Regressions are performed independently within each genomic partition.
#   Therefore, the complete analysis contains:
#
#       3 features × 6 contrasts × 2 genomic partitions = 36 regressions.
#
#   Benjamini-Hochberg correction is applied across all 36 regression tests.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

# 1) Libraries ------------------------------------------------------------------

library(tidyverse)


# 2) Inputs and outputs ----------------------------------------------------------

METH_FILE <- "/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/04_GE/data/00_GE-and-DGE/gene_methylation_means/final_tables/meth.combined.tsv"

DEG_FILE <- "/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/04_GE/data/00_GE-and-DGE/Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv"

VMR_FILE <- paste0(
  "/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/",
  "02_differentialMethylation/05_DGE/data/VMR_diff0.overlap_to_mRNA.tsv"
)

OUTDIR <- "../plots/feature_methylation-vs-DGE_atVMR_FDR"

dir.create(
  OUTDIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# Features to analyze
FEATURES <- c(
  "u2000",
  "geneBody",
  "d2000"
)


# Contrast definitions
#
# morph_A and morph_B are only used to select the corresponding
# morph-level methylation columns.
#
# Because X uses the absolute methylation difference, contrast direction
# does not affect X.

CONTRAST_INFO <- tribble(
  ~contrast,    ~morph_A,       ~morph_B,       ~l2fc_column,
  "i_vs_y_01",  "immaculata",   "yellow",       "i_vs_y_01_l2fc",
  "p_vs_i_01",  "parae",        "immaculata",   "p_vs_i_01_l2fc",
  "y_vs_p_01",  "yellow",       "parae",        "y_vs_p_01_l2fc",
  "f_vs_i",     "female",       "immaculata",   "f_vs_i_l2fc",
  "f_vs_p",     "female",       "parae",        "f_vs_p_l2fc",
  "f_vs_y",     "female",       "yellow",       "f_vs_y_l2fc"
)


# 3) Read input tables -----------------------------------------------------------

meth <- read_tsv(
  METH_FILE,
  show_col_types = FALSE,
  na = c("NA", "NaN", "")
)

deg <- read_tsv(
  DEG_FILE,
  show_col_types = FALSE,
  na = c("NA", "NaN", "")
)

vmr_overlap <- read_tsv(
  VMR_FILE,
  show_col_types = FALSE,
  na = c("NA", "NaN", "")
)


# 4) Check required columns ------------------------------------------------------

required_meth_columns <- c(
  "chr",
  "start",
  "end",
  "gene",
  "feature",
  "chromosome_type",
  "meth.morph.female",
  "meth.morph.immaculata",
  "meth.morph.parae",
  "meth.morph.yellow"
)

required_deg_columns <- c(
  "gene_id",
  "gene_name",
  CONTRAST_INFO$l2fc_column
)


missing_meth_columns <- setdiff(
  required_meth_columns,
  colnames(meth)
)

missing_deg_columns <- setdiff(
  required_deg_columns,
  colnames(deg)
)


if (length(missing_meth_columns) > 0) {
  stop(
    "Missing required columns from methylation table:\n",
    paste(missing_meth_columns, collapse = "\n")
  )
}


if (length(missing_deg_columns) > 0) {
  stop(
    "Missing required columns from DEG table:\n",
    paste(missing_deg_columns, collapse = "\n")
  )
}


# 5) Select requested features ---------------------------------------------------

meth_selected <- meth %>%
  filter(feature %in% FEATURES)


message(
  "Rows retained for selected features: ",
  nrow(meth_selected)
)

# 5.1) Generate list of genes overlapping VMRs ------------------------------------

# The VMR overlap table can contain multiple rows for the same gene because:
#   - multiple VMRs may overlap the same gene
#   - one VMR may potentially contribute multiple overlap records
#
# For the methylation–DGE analysis, only gene membership is required.
# Therefore, retain one unique gene ID per VMR-associated gene.

vmr_genes <- vmr_overlap %>%
  
  filter(
    !is.na(gene),
    gene != ""
  ) %>%
  
  distinct(
    gene
  ) %>%
  
  arrange(
    gene
  )


message(
  "Unique genes overlapping at least one VMR: ",
  nrow(vmr_genes)
)


write_tsv(
  vmr_genes,
  file.path(
    OUTDIR,
    "VMR-overlapping_genes.tsv"
  )
)

# # 6) Check gene + feature uniqueness --------------------------------------------
# 
# # One gene should contribute only one observation for each feature.
# #
# # If duplicates exist, they are reported and the script stops. This avoids
# # artificially weighting genes multiple times in the regression.
# 
# meth_duplicates <- meth_selected %>%
#   count(
#     gene,
#     feature,
#     name = "n_rows"
#   ) %>%
#   filter(n_rows > 1) %>%
#   arrange(
#     feature,
#     desc(n_rows),
#     gene
#   )
# 
# 
# write_tsv(
#   meth_duplicates,
#   file.path(
#     OUTDIR,
#     "methylation_gene-feature_duplicates.tsv"
#   )
# )
# 
# 
# if (nrow(meth_duplicates) > 0) {
#   
#   message("")
#   message(
#     "Found ",
#     nrow(meth_duplicates),
#     " duplicated gene + feature combinations."
#   )
#   
#   message(
#     "Duplicates were written to: ",
#     file.path(
#       OUTDIR,
#       "methylation_gene-feature_duplicates.tsv"
#     )
#   )
#   
#   stop(
#     "Methylation gene + feature combinations are not unique. ",
#     "Resolve these duplicates before continuing."
#   )
#   
# } else {
#   
#   message(
#     "Methylation uniqueness check passed: ",
#     "gene + feature combinations are unique."
#   )
# }
# 
# # 6) Diagnose gene + feature duplicates -----------------------------------------
# 
# # Number of rows associated with each gene × feature combination
# meth_duplicate_counts <- meth_selected %>%
#   count(
#     gene,
#     feature,
#     name = "n_rows"
#   ) %>%
#   filter(n_rows > 1) %>%
#   arrange(
#     feature,
#     desc(n_rows),
#     gene
#   )
# 
# 
# write_tsv(
#   meth_duplicate_counts,
#   file.path(
#     OUTDIR,
#     "methylation_gene-feature_duplicate_counts.tsv"
#   )
# )
# 
# 
# # Summary by feature
# duplicate_summary <- meth_selected %>%
#   group_by(
#     feature
#   ) %>%
#   summarise(
#     total_rows = n(),
#     unique_genes = n_distinct(gene),
#     
#     gene_feature_combinations = n_distinct(
#       paste(gene, feature, sep = "|")
#     ),
#     
#     duplicated_gene_feature_combinations = sum(
#       duplicated(gene) |
#         duplicated(gene, fromLast = TRUE)
#     ),
#     
#     max_rows_per_gene = max(
#       table(gene)
#     ),
#     
#     .groups = "drop"
#   )
# 
# 
# write_tsv(
#   duplicate_summary,
#   file.path(
#     OUTDIR,
#     "methylation_gene-feature_duplicate_summary.tsv"
#   )
# )
# 
# 
# # Export the ACTUAL rows involved in duplicated gene × feature combinations
# meth_duplicate_rows <- meth_selected %>%
#   semi_join(
#     meth_duplicate_counts,
#     by = c(
#       "gene",
#       "feature"
#     )
#   ) %>%
#   arrange(
#     feature,
#     gene,
#     chr,
#     start,
#     end
#   )
# 
# 
# write_tsv(
#   meth_duplicate_rows,
#   file.path(
#     OUTDIR,
#     "methylation_gene-feature_duplicate_rows.tsv"
#   )
# )
# 
# 
# # Determine whether duplicated gene-feature rows also have duplicated coordinates
# duplicate_coordinate_summary <- meth_duplicate_rows %>%
#   group_by(
#     gene,
#     feature
#   ) %>%
#   summarise(
#     n_rows = n(),
#     
#     n_unique_intervals = n_distinct(
#       paste(
#         chr,
#         start,
#         end,
#         sep = ":"
#       )
#     ),
#     
#     .groups = "drop"
#   ) %>%
#   mutate(
#     duplicate_type = case_when(
#       n_unique_intervals == 1 ~ "same_interval_repeated",
#       n_unique_intervals > 1 ~ "multiple_intervals",
#       TRUE ~ "unknown"
#     )
#   )
# 
# 
# write_tsv(
#   duplicate_coordinate_summary,
#   file.path(
#     OUTDIR,
#     "methylation_gene-feature_duplicate_coordinate_summary.tsv"
#   )
# )
# 
# 
# # Global summary
# print(
#   duplicate_summary
# )
# 
# print(
#   duplicate_coordinate_summary %>%
#     count(
#       feature,
#       duplicate_type
#     )
# )
# 
# 
# if (nrow(meth_duplicate_counts) > 0) {
#   
#   stop(
#     "Gene + feature combinations are not unique. ",
#     "Duplicate diagnostics were exported for inspection."
#   )
# }


# 6) Collapse multiple intervals to one gene × feature --------------------------

# A gene may contain multiple annotated intervals for the same feature,
# for example due to multiple transcript isoforms.
#
# For this analysis, all intervals belonging to the same gene × feature
# are collapsed by taking the arithmetic mean of their morph-level
# methylation values.
#
# This produces exactly one methylation value per:
#
#     gene × feature × morph
#
# Coverage is NOT used for weighting.
#
# If some intervals contain NA methylation values, the mean is calculated
# from the available intervals. If ALL intervals are NA for a morph, the
# collapsed value remains NA.


# First report how many intervals exist per gene × feature
meth_duplicate_counts <- meth_selected %>%
  
  count(
    gene,
    feature,
    name = "n_intervals"
  ) %>%
  
  arrange(
    feature,
    desc(n_intervals),
    gene
  )


write_tsv(
  meth_duplicate_counts,
  file.path(
    OUTDIR,
    "methylation_gene-feature_interval_counts.tsv"
  )
)


# Report only gene × feature combinations with >1 interval
meth_duplicates <- meth_duplicate_counts %>%
  
  filter(
    n_intervals > 1
  )


write_tsv(
  meth_duplicates,
  file.path(
    OUTDIR,
    "methylation_gene-feature_duplicates.tsv"
  )
)


message(
  "Gene × feature combinations with multiple intervals: ",
  nrow(meth_duplicates)
)


# Check that all intervals assigned to the same gene × feature
# occur on the same chromosome.
#
# This should normally be true. If not, stop because that would indicate
# a potentially problematic gene annotation.

chromosome_conflicts <- meth_selected %>%
  
  group_by(
    gene,
    feature
  ) %>%
  
  summarise(
    n_chr = n_distinct(chr),
    chromosomes = paste(
      sort(unique(chr)),
      collapse = ","
    ),
    .groups = "drop"
  ) %>%
  
  filter(
    n_chr > 1
  )


write_tsv(
  chromosome_conflicts,
  file.path(
    OUTDIR,
    "methylation_gene-feature_chromosome_conflicts.tsv"
  )
)


if (nrow(chromosome_conflicts) > 0) {
  
  stop(
    "Some gene × feature combinations occur on multiple chromosomes. ",
    "See methylation_gene-feature_chromosome_conflicts.tsv"
  )
  
} else {
  
  message(
    "Chromosome consistency check passed."
  )
}


# Function:
# mean of available methylation values;
# return NA if all values are missing.

mean_nonNA <- function(x) {
  
  if (all(is.na(x))) {
    return(NA_real_)
  }
  
  mean(
    x,
    na.rm = TRUE
  )
}


# Collapse to one row per gene × feature
meth_gene_feature <- meth_selected %>%
  
  group_by(
    gene,
    feature
  ) %>%
  
  summarise(
    
    chr = first(chr),
    
    chromosome_type = first(chromosome_type),
    
    # Number of annotated intervals being collapsed
    n_intervals = n(),
    
    # Coordinate range retained only for reference.
    # This is NOT interpreted as a newly defined genomic feature.
    feature_span_start = min(
      start,
      na.rm = TRUE
    ),
    
    feature_span_end = max(
      end,
      na.rm = TRUE
    ),
    
    # Mean morph-level methylation
    meth.morph.female = mean_nonNA(
      meth.morph.female
    ),
    
    meth.morph.immaculata = mean_nonNA(
      meth.morph.immaculata
    ),
    
    meth.morph.parae = mean_nonNA(
      meth.morph.parae
    ),
    
    meth.morph.yellow = mean_nonNA(
      meth.morph.yellow
    ),
    
    # Number of intervals contributing methylation information
    # for each morph
    n_nonNA.female = sum(
      !is.na(meth.morph.female)
    ),
    
    n_nonNA.immaculata = sum(
      !is.na(meth.morph.immaculata)
    ),
    
    n_nonNA.parae = sum(
      !is.na(meth.morph.parae)
    ),
    
    n_nonNA.yellow = sum(
      !is.na(meth.morph.yellow)
    ),
    
    .groups = "drop"
  )

# Assign genomic partition -------------------------------------------------------

meth_gene_feature <- meth_gene_feature %>%
  mutate(
    genomic_partition = case_when(
      chr == "Parae_12" ~ "Chromosome_12",
      chromosome_type == "Autosome" ~ "Autosome",
      TRUE ~ NA_character_
    )
  )


# Report excluded chromosome classes, if any
excluded_partition <- meth_gene_feature %>%
  filter(
    is.na(genomic_partition)
  )


message(
  "Gene × feature combinations outside Autosome/Parae_12: ",
  nrow(excluded_partition)
)


write_tsv(
  excluded_partition,
  file.path(
    OUTDIR,
    "methylation_gene-feature.excluded_chromosomes.tsv"
  )
)


# Keep only autosomes and chromosome 12
meth_gene_feature <- meth_gene_feature %>%
  filter(
    !is.na(genomic_partition)
  )


meth_gene_feature <- meth_gene_feature %>%
  mutate(
    genomic_partition = factor(
      genomic_partition,
      levels = c(
        "Autosome",
        "Chromosome_12"
      )
    )
  )



# Export collapsed methylation table
write_tsv(
  meth_gene_feature,
  file.path(
    OUTDIR,
    "methylation_gene-feature.collapsed_mean.tsv"
  )
)


# Verify that collapsing produced exactly one row per gene × feature
stopifnot(
  !anyDuplicated(
    meth_gene_feature[
      c(
        "gene",
        "feature"
      )
    ]
  )
)


message(
  "Collapsed methylation table: ",
  nrow(meth_gene_feature),
  " unique gene × feature combinations."
)

# 7) Check DEG gene_id uniqueness ------------------------------------------------

# A duplicated gene_id in the DEG table would also duplicate observations
# after joining to methylation data.

deg_duplicates <- deg %>%
  count(
    gene_id,
    name = "n_rows"
  ) %>%
  filter(n_rows > 1) %>%
  arrange(
    desc(n_rows),
    gene_id
  )


write_tsv(
  deg_duplicates,
  file.path(
    OUTDIR,
    "DEG_gene_id_duplicates.tsv"
  )
)


if (nrow(deg_duplicates) > 0) {
  
  message("")
  message(
    "Found ",
    nrow(deg_duplicates),
    " duplicated gene_id values in the DEG table."
  )
  
  message(
    "Duplicates were written to: ",
    file.path(
      OUTDIR,
      "DEG_gene_id_duplicates.tsv"
    )
  )
  
  stop(
    "DEG gene_id values are not unique. ",
    "Resolve these duplicates before continuing."
  )
  
} else {
  
  message(
    "DEG uniqueness check passed: gene_id values are unique."
  )
}


# 8) Build feature × contrast analysis table ------------------------------------

# Each row in the resulting table represents:
#
#     one gene
#     × one feature
#     × one contrast
#
# Missing methylation or log2FC information is removed only from the
# corresponding contrast.

analysis_list <- vector(
  "list",
  nrow(CONTRAST_INFO)
)


for (i in seq_len(nrow(CONTRAST_INFO))) {
  
  contrast <- CONTRAST_INFO$contrast[i]
  morph_A  <- CONTRAST_INFO$morph_A[i]
  morph_B  <- CONTRAST_INFO$morph_B[i]
  l2fc_col <- CONTRAST_INFO$l2fc_column[i]
  
  meth_A_col <- paste0(
    "meth.morph.",
    morph_A
  )
  
  meth_B_col <- paste0(
    "meth.morph.",
    morph_B
  )
  
  
  # analysis_list[[i]] <- meth_selected %>%
    # transmute(
    #   chr,
    #   start,
    #   end,
    #   gene,
    #   feature,
    #   chromosome_type,
  
  analysis_list[[i]] <- meth_gene_feature %>%
    
    transmute(
      chr,
      feature_span_start,
      feature_span_end,
      gene,
      feature,
      chromosome_type,
      genomic_partition,
      n_intervals,
      
      contrast = contrast,
      morph_A = morph_A,
      morph_B = morph_B,
      
      meth_A = .data[[meth_A_col]],
      meth_B = .data[[meth_B_col]],
      
      abs_delta_meth = abs(
        .data[[meth_A_col]] -
          .data[[meth_B_col]]
      )
    ) %>%
    
    left_join(
      deg %>%
        transmute(
          gene_id,
          gene_name,
          log2FC = .data[[l2fc_col]]
        ),
      by = c(
        "gene" = "gene_id"
      )
    ) %>%
    
    mutate(
      abs_log2FC = abs(log2FC)
    )
}


analysis_all <- bind_rows(
  analysis_list
)


# Preserve requested ordering
analysis_all <- analysis_all %>%
  mutate(
    feature = factor(
      feature,
      levels = FEATURES
    ),
    contrast = factor(
      contrast,
      levels = CONTRAST_INFO$contrast
    )
  )


# 9) Summarize missing observations ---------------------------------------------

# Report how many genes are available before and after removing missing
# contrast-specific information.

sample_counts <- analysis_all %>%
  
  group_by(
    genomic_partition,
    feature,
    contrast
  ) %>%
  
  summarise(
    total_feature_genes = n(),
    
    missing_methylation = sum(
      is.na(meth_A) |
        is.na(meth_B)
    ),
    
    missing_log2FC = sum(
      is.na(log2FC)
    ),
    
    complete_cases = sum(
      !is.na(meth_A) &
        !is.na(meth_B) &
        !is.na(log2FC)
    ),
    
    .groups = "drop"
  )


write_tsv(
  sample_counts,
  file.path(
    OUTDIR,
    "feature_methylation-vs-DGE.sample_counts.tsv"
  )
)

# Report effect of VMR gene filtering -------------------------------------------

vmr_filter_summary <- analysis_all %>%
  
  mutate(
    overlaps_VMR = gene %in% vmr_genes$gene
  ) %>%
  
  group_by(
    genomic_partition,
    feature,
    contrast
  ) %>%
  
  summarise(
    genes_before_VMR_filter = n_distinct(gene),
    
    genes_overlapping_VMR = n_distinct(
      gene[overlaps_VMR]
    ),
    
    percent_overlapping_VMR = 100 *
      genes_overlapping_VMR /
      genes_before_VMR_filter,
    
    .groups = "drop"
  )


write_tsv(
  vmr_filter_summary,
  file.path(
    OUTDIR,
    "VMR_gene-filter_summary.tsv"
  )
)



# 10) Remove contrast-specific missing observations -----------------------------

# analysis_complete <- analysis_all %>%
#   
#   filter(
#     !is.na(meth_A),
#     !is.na(meth_B),
#     !is.na(log2FC),
#     !is.na(abs_delta_meth),
#     !is.na(abs_log2FC)
#   )

analysis_complete <- analysis_all %>%
  
  # Keep only genes overlapping at least one VMR
  semi_join(
    vmr_genes,
    by = "gene"
  ) %>%
  
  # Remove missing information independently for each contrast
  filter(
    !is.na(meth_A),
    !is.na(meth_B),
    !is.na(log2FC),
    !is.na(abs_delta_meth),
    !is.na(abs_log2FC)
  )

write_tsv(
  analysis_complete,
  file.path(
    OUTDIR,
    "feature_methylation-vs-DGE.long.tsv"
  )
)


message(
  "Complete gene × feature × contrast observations: ",
  nrow(analysis_complete)
)


# 11) Run linear regressions -----------------------------------------------------

# Model:
#
#     |log2FC| = beta0 + beta1 * |delta methylation|
#
# Statistics exported:
#     n
#     slope
#     Pearson r
#     R²
#     P
#     BH-adjusted P
#
# R², raw P, and n are shown in the figure.

regression_stats <- analysis_complete %>%
  
  group_by(
    genomic_partition,
    feature,
    contrast
  ) %>%
  
  group_modify(
    ~ {
      
      dat <- .x
      
      n_obs <- nrow(dat)
      
      
      # Regression requires at least 3 observations and variation in X and Y.
      if (
        n_obs < 3 ||
        n_distinct(dat$abs_delta_meth) < 2 ||
        n_distinct(dat$abs_log2FC) < 2
      ) {
        
        return(
          tibble(
            n = n_obs,
            slope = NA_real_,
            pearson_r = NA_real_,
            r_squared = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      
      fit <- lm(
        abs_log2FC ~ abs_delta_meth,
        data = dat
      )
      
      fit_summary <- summary(fit)
      
      
      tibble(
        n = n_obs,
        
        slope = unname(
          coef(fit)[["abs_delta_meth"]]
        ),
        
        pearson_r = cor(
          dat$abs_delta_meth,
          dat$abs_log2FC,
          method = "pearson"
        ),
        
        r_squared = fit_summary$r.squared,
        
        p_value = coef(
          fit_summary
        )[
          "abs_delta_meth",
          "Pr(>|t|)"
        ]
      )
    }
  ) %>%
  
  ungroup() %>%
  
  # Correct all feature × contrast regression tests together.
  mutate(
    p_adj_BH = p.adjust(
      p_value,
      method = "BH"
    )
  )


write_tsv(
  regression_stats,
  file.path(
    OUTDIR,
    "feature_methylation-vs-DGE.regression_stats.tsv"
  )
)


# 12) Function for clean plot P-value labels ------------------------------------

format_plot_p <- function(x) {
  
  case_when(
    is.na(x) ~ "NA",
    x < 0.0001 ~ "<0.0001",
    TRUE ~ formatC(
      x,
      format = "f",
      digits = 4
    )
  )
}


# 13) Create annotation table ----------------------------------------------------

plot_annotations <- regression_stats %>%
  
  mutate(
    label = paste0(
      "R² = ",
      if_else(
        is.na(r_squared),
        "NA",
        formatC(
          r_squared,
          format = "f",
          digits = 3
        )
      ),
      # "\nP = ",
      # format_plot_p(p_value),
      "\nFDR = ",
      format_plot_p(p_adj_BH),
      "\nn = ",
      n
    )
  )


# 14) Plot each feature ----------------------------------------------------------

# facet_wrap() uses fixed scales by default.
#
# Therefore all six contrasts within a feature use the same X and Y scale,
# allowing direct visual comparison across contrasts.

# for (feature_i in FEATURES) {
#   
#   plot_data <- analysis_complete %>%
#     filter(
#       feature == feature_i
#     )
#   
#   annotation_data <- plot_annotations %>%
#     filter(
#       feature == feature_i
#     )
#   
#   
#   p <- ggplot(
#     plot_data,
#     aes(
#       x = abs_delta_meth,
#       y = abs_log2FC
#     )
#   ) +
#     
#     geom_point(
#       alpha = 0.40,
#       size = 1.4
#     ) +
#     
#     geom_smooth(
#       method = "lm",
#       formula = y ~ x,
#       se = TRUE,
#       linewidth = 0.8
#     ) +
#     
#     geom_text(
#       data = annotation_data,
#       aes(
#         x = Inf,
#         y = Inf,
#         label = label
#       ),
#       inherit.aes = FALSE,
#       hjust = 1.08,
#       vjust = 1.15,
#       size = 3.4
#     ) +
#     
#     facet_wrap(
#       ~ contrast,
#       ncol = 3,
#       scales = "fixed"
#     ) +
#     
#     labs(
#       title = paste0(
#         "Feature: ",
#         feature_i
#       ),
#       x = expression(
#         "Absolute methylation difference (" *
#           "|" * Delta * " methylation" * "|" * ")"
#       ),
#       y = expression(
#         "Absolute differential expression (" *
#           "|" * log[2] * "FC" * "|" * ")"
#       )
#     ) +
#     
#     theme_bw(
#       base_size = 12
#     ) +
#     
#     theme(
#       panel.grid.minor = element_blank(),
#       
#       strip.text = element_text(
#         face = "bold"
#       ),
#       
#       plot.title = element_text(
#         face = "bold",
#         hjust = 0.5
#       )
#     )
#   
#   
#   # PDF
#   ggsave(
#     filename = file.path(
#       OUTDIR,
#       paste0(
#         feature_i,
#         ".absDeltaMeth-vs-absLFC.pdf"
#       )
#     ),
#     plot = p,
#     width = 11,
#     height = 7.5
#   )
#   
#   
#   # PNG
#   ggsave(
#     filename = file.path(
#       OUTDIR,
#       paste0(
#         feature_i,
#         ".absDeltaMeth-vs-absLFC.png"
#       )
#     ),
#     plot = p,
#     width = 11,
#     height = 7.5,
#     dpi = 300
#   )
# }

GENOMIC_PARTITIONS <- c(
  "Autosome",
  "Chromosome_12"
)


for (partition_i in GENOMIC_PARTITIONS) {
  
  for (feature_i in FEATURES) {
    
    plot_data <- analysis_complete %>%
      filter(
        genomic_partition == partition_i,
        feature == feature_i
      )
    
    
    annotation_data <- plot_annotations %>%
      filter(
        genomic_partition == partition_i,
        feature == feature_i
      )
    
    
    p <- ggplot(
      plot_data,
      aes(
        x = abs_delta_meth,
        y = abs_log2FC
      )
    ) +
      
      geom_point(
        alpha = 0.40,
        size = 1.4
      ) +
      
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = TRUE,
        linewidth = 0.8
      ) +
      
      geom_text(
        data = annotation_data,
        aes(
          x = Inf,
          y = Inf,
          label = label
        ),
        inherit.aes = FALSE,
        hjust = 1.08,
        vjust = 1.15,
        size = 3.4
      ) +
      
      facet_wrap(
        ~ contrast,
        ncol = 3,
        scales = "fixed"
      ) +
      
      labs(
        title = paste0(
          "Feature: ",
          feature_i,
          " — ",
          partition_i
        ),
        x = expression(
          "Absolute methylation difference (" *
            "|" * Delta * " methylation" * "|" * ")"
        ),
        y = expression(
          "Absolute differential expression (" *
            "|" * log[2] * "FC" * "|" * ")"
        )
      ) +
      
      theme_bw(
        base_size = 12
      ) +
      
      theme(
        panel.grid.minor = element_blank(),
        
        strip.text = element_text(
          face = "bold"
        ),
        
        plot.title = element_text(
          face = "bold",
          hjust = 0.5
        )
      )
    
    
    # PDF
    ggsave(
      filename = file.path(
        OUTDIR,
        paste0(
          partition_i,
          ".",
          feature_i,
          ".absDeltaMeth-vs-absLFC.pdf"
        )
      ),
      plot = p,
      width = 11,
      height = 7.5
    )
    
    
    # PNG
    ggsave(
      filename = file.path(
        OUTDIR,
        paste0(
          partition_i,
          ".",
          feature_i,
          ".absDeltaMeth-vs-absLFC.png"
        )
      ),
      plot = p,
      width = 11,
      height = 7.5,
      dpi = 300
    )
  }
}

# 15) Print regression results ---------------------------------------------------

print(
  regression_stats %>%
    arrange(
      feature,
      contrast
    )
)


message("")
message("Analysis completed successfully.")
message("Outputs written to: ", OUTDIR)