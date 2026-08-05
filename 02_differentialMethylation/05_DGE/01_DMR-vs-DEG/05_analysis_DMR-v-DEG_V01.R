#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
# DATE:       2026-03-18
# AUTHOR:     MZF & ChatGPT
# SCRIPT:     03_DMR_DE_geneLevel_firstPass.R
# VERSION:    01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Test whether differential DNA methylation is associated with differential
#   gene expression across real pairwise morph comparisons in P. parae muscle.
#
# QUESTION:
#   For each real pairwise comparison, do genes overlapping one or more DMRs:
#     (1) show greater expression divergence than genes without DMR overlap?
#     (2) show larger expression divergence when the magnitude of methylation
#         difference across their overlapping DMRs is larger?
#
# APPROACH:
#   Build one analysis row per gene x comparison by:
#     - using DESeq2 long-format results as the gene-expression universe
#     - restricting to real morph comparisons only
#     - collapsing all overlapping DMRs for each gene x comparison into summary
#       methylation metrics (presence/absence, count, mean/max/weighted |diff|)
#     - avoiding pseudoreplication by not treating multiple DMRs from the same
#       gene as independent expression observations
#   Then perform:
#     A) Wilcoxon tests on |shrunken log2FC| for DMR-overlap vs no-overlap genes
#     B) Spearman correlations and linear models relating |shrunken log2FC| to
#        gene-level methylation-difference summaries among DMR-overlap genes
#     C) Optional signed exploratory analyses after aligning DMR direction to the
#        DESeq2 log2FC sign convention
#
# INPUTS:
#   1) DESeq2_allComparisons.tsv
#      - long-format DESeq2 results with comparison identity, group_ref/group_alt,
#        baseMean, raw log2FC, and shrunken log2FC
#   2) compiled_intersect_DEGs_AUTO_CPMfilter.tsv
#      - compiled gene–DMR overlap table with DMR coordinates, methylation values,
#        differential methylation statistics, comparison labels, and gene IDs
#
# OUTPUTS:
#   1) DMR_DE_gene_level_firstpass.tsv
#      - gene x comparison analysis table with DE and summarized DMR metrics
#   2) DMR_DE_anyOverlap_wilcox.tsv
#      - per-comparison Wilcoxon tests for |log2FC| in DMR vs non-DMR genes
#   3) DMR_DE_doseResponse_abs.tsv
#      - per-comparison dose-response tests using absolute methylation difference
#   4) DMR_DE_doseResponse_signed.tsv
#      - optional exploratory signed analyses using sign-aligned methylation
#
# NOTES:
#   - First-pass analysis uses any gene overlap, regardless of feature class.
#   - Main response variable is |shrunken log2FC| from DESeq2.
#   - Signed analyses are exploratory because feature context is not yet split
#     into promoter, gene body, exon, intron, etc.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

# ============================================================================
# INPUTS
# ============================================================================
de_path <- "DESeq2_allComparisons.tsv"
ix_path <- "compiled_intersect_DEGs_CPMfilter_labeled.tsv"

real_comparisons <- c(
  "i_vs_y_01",
  "p_vs_i_01",
  "y_vs_p_01",
  "f_vs_i",
  "f_vs_p",
  "f_vs_y"
)


# feature_keep <- "CDS"
# feature_keep <- "d2000"
# feature_keep <- "exon"
# feature_keep <- "five_prime_UTR"
# feature_keep <- "five_prime_UTR_exon"
# feature_keep <- "gene"
# feature_keep <- "geneBody"
# feature_keep <- "geneBody_exon"
# feature_keep <- "geneBody_intron"
# feature_keep <- "intron"
# feature_keep <- "mRNA"
# feature_keep <- "three_prime_UTR"
# feature_keep <- "three_prime_UTR_exon"
# feature_keep <- "u2000"
feature_keep <- "any" # if commenting the feature filter in line 156

# ============================================================================
# 1) READ DE TABLE
# ============================================================================
de <- read_tsv(de_path, show_col_types = FALSE)

# Basic sanity checks
stopifnot(ncol(de) == 14)

required_de_cols <- c(
  "gene_id", "condition", "comparison", "group_ref", "group_alt",
  "tissue", "treatment", "baseMean", "log2FoldChange",
  "pvalue", "padj", "log2FoldChange_shrunk"
)
missing_de <- setdiff(required_de_cols, colnames(de))
if (length(missing_de) > 0) {
  stop("Missing DE columns: ", paste(missing_de, collapse = ", "))
}

de_use <- de %>%
  filter(
    condition %in% c("real_morph", "real_sex"),
    comparison %in% real_comparisons,
    tissue == "Muscle",
    treatment == "veh"
  ) %>%
  transmute(
    gene_id,
    comparison,
    group_ref,
    group_alt,
    baseMean,
    l2fc_raw    = log2FoldChange,
    l2fc_shrunk = log2FoldChange_shrunk,
    l2fc_use    = coalesce(log2FoldChange_shrunk, log2FoldChange),
    abs_l2fc    = abs(l2fc_use),
    pvalue,
    padj
  )

# ============================================================================
# 2) READ INTERSECT TABLE
# ============================================================================
ix <- read_tsv(ix_path, show_col_types = FALSE)

required_ix_cols <- c(
  "gene", "dmr_chr", "dmr_start", "dmr_end",
  "dmr_meth1", "dmr_meth2", "dmr_diff",
  "dmr_condition", "dmr_comparison"
)
missing_ix <- setdiff(required_ix_cols, colnames(ix))
if (length(missing_ix) > 0) {
  stop("Missing intersect columns: ", paste(missing_ix, collapse = ", "))
}

ix_use <- ix %>%
  filter(
    dmr_condition %in% c("real_morph", "real_sex"),
    dmr_comparison %in% real_comparisons,
    # feature == feature_keep # comment to exclude filtering
  ) %>%
  mutate(
    gene_id    = gene,
    comparison = dmr_comparison,
    dmr_width  = dmr_end - dmr_start + 1L,
    dmr_id     = paste(dmr_chr, dmr_start, dmr_end, sep = ":")
  ) %>%
  # Remove repeated feature-level hits for the same gene x DMR x comparison
  distinct(gene_id, comparison, dmr_id, .keep_all = TRUE) %>%
  # Align DMR sign to DE sign:
  # DE is group_alt - group_ref
  # DMR diff appears to be meth1 - meth2 = first group - second group
  # Therefore aligned methylation difference is second - first = -dmr_diff
  mutate(
    dmr_diff_aligned = -dmr_diff
  )

# Optional consistency check on the first few rows:
# dmr_diff should equal dmr_meth1 - dmr_meth2
check_diff <- ix_use %>%
  mutate(diff_check = dmr_meth1 - dmr_meth2) %>%
  summarise(max_abs_error = max(abs(dmr_diff - diff_check), na.rm = TRUE)) %>%
  pull(max_abs_error)

message("Max |dmr_diff - (dmr_meth1 - dmr_meth2)| = ", signif(check_diff, 6))

# ============================================================================
# 3) SUMMARIZE DMRs TO ONE ROW PER gene x comparison
# ============================================================================
dmr_gene <- ix_use %>%
  group_by(gene_id, comparison) %>%
  summarise(
    dmr_any = 1L,
    n_dmr   = n(),
    
    mean_abs_dmr_diff     = mean(abs(dmr_diff), na.rm = TRUE),
    max_abs_dmr_diff      = max(abs(dmr_diff), na.rm = TRUE),
    sum_abs_dmr_diff      = sum(abs(dmr_diff), na.rm = TRUE),
    weighted_abs_dmr_diff = weighted.mean(abs(dmr_diff), w = dmr_width, na.rm = TRUE),
    
    mean_signed_dmr_diff_raw     = mean(dmr_diff, na.rm = TRUE),
    weighted_signed_dmr_diff_raw = weighted.mean(dmr_diff, w = dmr_width, na.rm = TRUE),
    
    mean_signed_dmr_diff_aligned     = mean(dmr_diff_aligned, na.rm = TRUE),
    weighted_signed_dmr_diff_aligned = weighted.mean(dmr_diff_aligned, w = dmr_width, na.rm = TRUE),
    
    mean_dmr_meth1 = weighted.mean(dmr_meth1, w = dmr_width, na.rm = TRUE),
    mean_dmr_meth2 = weighted.mean(dmr_meth2, w = dmr_width, na.rm = TRUE),
    
    .groups = "drop"
  )

# ============================================================================
# 4) BUILD FINAL ANALYSIS TABLE
#    Keep all genes from DE universe; genes with no DMR get zeros
# ============================================================================
analysis_df <- de_use %>%
  left_join(dmr_gene, by = c("gene_id", "comparison")) %>%
  mutate(
    dmr_any = if_else(is.na(dmr_any), 0L, dmr_any),
    n_dmr   = if_else(is.na(n_dmr), 0L, n_dmr),
    
    across(
      c(
        mean_abs_dmr_diff, max_abs_dmr_diff, sum_abs_dmr_diff,
        weighted_abs_dmr_diff,
        mean_signed_dmr_diff_raw, weighted_signed_dmr_diff_raw,
        mean_signed_dmr_diff_aligned, weighted_signed_dmr_diff_aligned,
        mean_dmr_meth1, mean_dmr_meth2
      ),
      ~ replace_na(.x, 0)
    )
  )

# ============================================================================
# 5) TEST 1: ANY OVERLAP vs NO OVERLAP
#    Response = |shrunken LFC|
# ============================================================================
res_any <- analysis_df %>%
  group_by(comparison) %>%
  group_modify(~{
    d <- .x
    
    wt <- wilcox.test(abs_l2fc ~ dmr_any, data = d, exact = FALSE)
    
    tibble(
      n_genes = nrow(d),
      n_dmr_genes = sum(d$dmr_any == 1, na.rm = TRUE),
      median_abs_l2fc_noDMR = median(d$abs_l2fc[d$dmr_any == 0], na.rm = TRUE),
      median_abs_l2fc_DMR   = median(d$abs_l2fc[d$dmr_any == 1], na.rm = TRUE),
      p_wilcox = wt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(p_wilcox_fdr = p.adjust(p_wilcox, method = "BH"))

# ============================================================================
# 6) TEST 2: DOSE-RESPONSE AMONG DMR GENES
#    Response = |shrunken LFC|
#    Predictor = weighted absolute methylation difference
# ============================================================================
res_dose_abs <- analysis_df %>%
  filter(dmr_any == 1) %>%
  group_by(comparison) %>%
  group_modify(~{
    d <- .x
    
    ct <- suppressWarnings(
      cor.test(d$abs_l2fc, d$weighted_abs_dmr_diff,
               method = "spearman", exact = FALSE)
    )
    
    fit <- lm(
      abs_l2fc ~ weighted_abs_dmr_diff + log10(baseMean + 1) + n_dmr,
      data = d
    )
    
    coefs <- summary(fit)$coefficients
    
    tibble(
      n_dmr_genes = nrow(d),
      spearman_rho = unname(ct$estimate),
      p_spearman   = ct$p.value,
      beta_weighted_abs_dmr_diff = coefs["weighted_abs_dmr_diff", "Estimate"],
      p_lm_weighted_abs_dmr_diff = coefs["weighted_abs_dmr_diff", "Pr(>|t|)"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_spearman_fdr = p.adjust(p_spearman, method = "BH"),
    p_lm_fdr       = p.adjust(p_lm_weighted_abs_dmr_diff, method = "BH")
  )

# ============================================================================
# 7) TEST 3: EXPLORATORY SIGNED ANALYSIS
#    Response = signed shrunken LFC
#    Predictor = aligned signed methylation difference
#
#    Interpretation:
#    positive predictor  -> more methylation in group_alt
#    positive response   -> more expression in group_alt
#    negative beta       -> methylation opposes expression overall
# ============================================================================
res_dose_signed <- analysis_df %>%
  filter(dmr_any == 1) %>%
  group_by(comparison) %>%
  group_modify(~{
    d <- .x
    
    ct <- suppressWarnings(
      cor.test(d$l2fc_use, d$weighted_signed_dmr_diff_aligned,
               method = "spearman", exact = FALSE)
    )
    
    fit <- lm(
      l2fc_use ~ weighted_signed_dmr_diff_aligned + log10(baseMean + 1) + n_dmr,
      data = d
    )
    
    coefs <- summary(fit)$coefficients
    
    tibble(
      n_dmr_genes = nrow(d),
      spearman_rho = unname(ct$estimate),
      p_spearman   = ct$p.value,
      beta_weighted_signed_dmr_diff_aligned =
        coefs["weighted_signed_dmr_diff_aligned", "Estimate"],
      p_lm_weighted_signed_dmr_diff_aligned =
        coefs["weighted_signed_dmr_diff_aligned", "Pr(>|t|)"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_spearman_fdr = p.adjust(p_spearman, method = "BH"),
    p_lm_fdr       = p.adjust(p_lm_weighted_signed_dmr_diff_aligned, method = "BH")
  )

# ============================================================================
# 8) WRITE OUTPUTS
# ============================================================================
write_tsv(analysis_df,      paste0(feature_keep,".DMR_DE_gene_level_firstpass.tsv"))
write_tsv(res_any,          paste0(feature_keep,".DMR_DE_anyOverlap_wilcox.tsv"))
write_tsv(res_dose_abs,     paste0(feature_keep,".DMR_DE_doseResponse_abs.tsv"))
write_tsv(res_dose_signed,  paste0(feature_keep,".DMR_DE_doseResponse_signed.tsv"))

message("Done.")
