#!/usr/bin/env Rscript

# DATE: 20260603
# AUTHOR: MZF
# GOAL:
#   Test whether cDMR / VMR occupancy in mRNA features differs from
#   expected occupancy based on genome gene coverage.
#
# TESTS:
#   1) Pooled genome-wide chi-square test
#   2) Chromosome-aware expected chi-square test
#
# USAGE:
#   Rscript 04_test_cDMR_gene_occupancy_chisq.R \
#     cDMRs_diff0.overlap_summary.byChromosome.tsv \
#     mRNA_geneModel_coverage.summary.byChromosome.tsv \
#     WG
#
#   Rscript 04_test_cDMR_gene_occupancy_chisq.R \
#     cDMRs_diff0.overlap_summary.byChromosome.tsv \
#     mRNA_geneModel_coverage.summary.byChromosome.tsv \
#     AUTO

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

# cdmr_tsv <- args[1]
# gene_tsv <- args[2]
# PARTITION <- args[3]

cdmr_tsv <- "../plots_0/cDMRs_diff0.overlap_summary.byChromosome.tsv"
gene_tsv <- "../plots_0/mRNA_geneModel_coverage_output/mRNA_geneModel_coverage.summary.byChromosome.tsv"
PARTITION <- "WG"
# PARTITION <- "AUTO"


OUTDIR <- "cDMR_gene_occupancy_chisq_output"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Load tables
# ------------------------------------------------------------

cdmr <- fread(cdmr_tsv, sep = "\t", header = TRUE)
gene <- fread(gene_tsv, sep = "\t", header = TRUE)

# ------------------------------------------------------------
# Merge tables
# ------------------------------------------------------------

dt <- merge(
  cdmr,
  gene,
  by = c("chromosome", "chromosome_type")
)

# ------------------------------------------------------------
# Filter genome partition
# ------------------------------------------------------------

if (PARTITION == "WG") {
  dt_use <- copy(dt)
}

if (PARTITION == "AUTO") {
  dt_use <- dt[chromosome_type == "Autosome"]
}

# ------------------------------------------------------------
# 1) Pooled genome-wide chi-square test
# ------------------------------------------------------------

observed_intragenic <- sum(dt_use$count_overlap_to_mRNA)
observed_intergenic <- sum(dt_use$count_intergenic)
observed_total <- observed_intragenic + observed_intergenic

genome_mRNA_bp <- sum(dt_use$mRNA_covered_bp)
genome_non_mRNA_bp <- sum(dt_use$non_mRNA_bp)
genome_total_bp <- genome_mRNA_bp + genome_non_mRNA_bp

expected_p_intragenic <- genome_mRNA_bp / genome_total_bp
expected_p_intergenic <- genome_non_mRNA_bp / genome_total_bp

expected_intragenic <- observed_total * expected_p_intragenic
expected_intergenic <- observed_total * expected_p_intergenic

pooled_chisq <- sum(
  (observed_intragenic - expected_intragenic)^2 / expected_intragenic,
  (observed_intergenic - expected_intergenic)^2 / expected_intergenic
)

pooled_df <- 1
pooled_pvalue <- pchisq(pooled_chisq, df = pooled_df, lower.tail = FALSE)

pooled_report <- data.table(
  test = "pooled_genome_wide_chisq",
  partition = PARTITION,
  n_chromosomes = nrow(dt_use),
  observed_intragenic = observed_intragenic,
  observed_intergenic = observed_intergenic,
  observed_total = observed_total,
  expected_p_intragenic = expected_p_intragenic,
  expected_p_intergenic = expected_p_intergenic,
  expected_intragenic = expected_intragenic,
  expected_intergenic = expected_intergenic,
  chisq = pooled_chisq,
  df = pooled_df,
  pvalue = pooled_pvalue
)

pooled_outfile <- file.path(
  OUTDIR,
  paste0("cDMR_gene_occupancy.pooled_chisq.", PARTITION, ".tsv")
)

fwrite(
  pooled_report,
  pooled_outfile,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote: ", pooled_outfile, "\n", sep = "")

# ------------------------------------------------------------
# 2) Chromosome-aware expected chi-square test
# ------------------------------------------------------------

dt_use[, expected_p_intragenic_chr := mRNA_covered_bp / chromosome_length_bp]
dt_use[, expected_p_intergenic_chr := non_mRNA_bp / chromosome_length_bp]

dt_use[, expected_intragenic_chr := total_cDMRs * expected_p_intragenic_chr]
dt_use[, expected_intergenic_chr := total_cDMRs * expected_p_intergenic_chr]

dt_use[, chisq_component_intragenic := 
         (count_overlap_to_mRNA - expected_intragenic_chr)^2 / expected_intragenic_chr]

dt_use[, chisq_component_intergenic := 
         (count_intergenic - expected_intergenic_chr)^2 / expected_intergenic_chr]

dt_use[, chisq_component_total := 
         chisq_component_intragenic + chisq_component_intergenic]

chromosome_aware_chisq <- sum(dt_use$chisq_component_total)

chromosome_aware_df <- nrow(dt_use)

chromosome_aware_pvalue <- pchisq(
  chromosome_aware_chisq,
  df = chromosome_aware_df,
  lower.tail = FALSE
)

chromosome_aware_report <- data.table(
  test = "chromosome_aware_expected_chisq",
  partition = PARTITION,
  n_chromosomes = nrow(dt_use),
  observed_intragenic = observed_intragenic,
  observed_intergenic = observed_intergenic,
  observed_total = observed_total,
  chisq = chromosome_aware_chisq,
  df = chromosome_aware_df,
  pvalue = chromosome_aware_pvalue
)

chromosome_aware_outfile <- file.path(
  OUTDIR,
  paste0("cDMR_gene_occupancy.chromosome_aware_chisq.", PARTITION, ".tsv")
)

fwrite(
  chromosome_aware_report,
  chromosome_aware_outfile,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote: ", chromosome_aware_outfile, "\n", sep = "")

# ------------------------------------------------------------
# Save chromosome-level expected/observed table
# ------------------------------------------------------------

chromosome_detail <- dt_use[, .(
  chromosome,
  chromosome_type,
  total_cDMRs,
  observed_intragenic = count_overlap_to_mRNA,
  observed_intergenic = count_intergenic,
  expected_p_intragenic_chr,
  expected_p_intergenic_chr,
  expected_intragenic_chr,
  expected_intergenic_chr,
  chisq_component_intragenic,
  chisq_component_intergenic,
  chisq_component_total
)]

chromosome_detail_outfile <- file.path(
  OUTDIR,
  paste0("cDMR_gene_occupancy.chromosome_aware_details.", PARTITION, ".tsv")
)

fwrite(
  chromosome_detail,
  chromosome_detail_outfile,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote: ", chromosome_detail_outfile, "\n", sep = "")