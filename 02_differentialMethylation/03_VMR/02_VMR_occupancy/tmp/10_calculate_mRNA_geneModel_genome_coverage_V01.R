#!/usr/bin/env Rscript

# DATE: 20260521
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Calculate the percentage of the genome covered by gene-model features,
#   defined here as annotation rows where the feature column is equal to "mRNA".
#
#   The script calculates genome coverage at two scales:
#
#     1) Genome partitions:
#          a) Whole genome
#          b) Autosomes only
#          c) X chromosome only / Chromosome 12 only
#
#     2) Per chromosome:
#          percentage of each chromosome covered by mRNA features
#
#   The script merges overlapping mRNA intervals before calculating coverage.
#   This is critical because multiple transcript isoforms or overlapping gene
#   models can otherwise cause the same genomic base to be counted more than once.
#
# INPUT:
#   1) BED-like annotation file with genomic features.
#      Expected columns:
#        col 1 = chromosome
#        col 2 = start
#        col 3 = end
#        col 4 = gene / feature ID
#        col 5 = feature type
#        col 6 = strand
#        col 7 = feature length
#
#      The script filters rows where:
#        feature type == "mRNA"
#
#   2) Genome index file in .fai format.
#      Expected columns:
#        col 1 = chromosome
#        col 2 = chromosome length
#
# OUTPUT FILES:
#   Tables:
#     - mRNA_geneModel_coverage.summary.byPartition.tsv
#     - mRNA_geneModel_coverage.summary.byChromosome.tsv
#     - mRNA_geneModel_coverage.merged_mRNA_intervals.bed
#
#   Figures:
#     - mRNA_geneModel_coverage.percentBars.stacked.byPartition.pdf
#     - mRNA_geneModel_coverage.percentBars.stacked.byChromosome.pdf
#
# DEFINITIONS:
#   mRNA covered:
#     genomic bases overlapping at least one merged mRNA feature.
#
#   non-mRNA covered:
#     genomic bases not overlapping any mRNA feature.
#
#   Whole genome:
#     all Parae_XX chromosomes in the .fai file.
#
#   Autosomes:
#     Parae_01–Parae_23 excluding Parae_12.
#
#   X chromosome / Chromosome 12:
#     Parae_12.
#
# DEPENDENCIES:
#   data.table
#   ggplot2
#
# USAGE:
#   Rscript 12_calculate_mRNA_geneModel_genome_coverage.R \
#     allFeatures.feature.promoters_u400-TSS-d200.bed \
#     WG.PparFemVer2024.fasta.fai \
#     mRNA_geneModel_coverage_output
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "\nUsage:\n",
    "  Rscript 10_calculate_mRNA_geneModel_genome_coverage_V01.R ",
    "<annotation_features.bed> <genome.fai> [output_dir]\n\n",
    "Example:\n",
    "  Rscript 10_calculate_mRNA_geneModel_genome_coverage_V01.R ",
    "allFeatures.feature.promoters_u400-TSS-d200.bed ",
    "WG.PparFemVer2024.fasta.fai ",
    "mRNA_geneModel_coverage_output\n"
  )
}

BED_FILE <- args[1]
FAI_FILE <- args[2]
OUTDIR   <- ifelse(length(args) >= 3, args[3], "mRNA_geneModel_coverage_output")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# User-defined parameters
# ------------------------------------------------------------

FEATURE_NAME <- "mRNA"
X_CHROM <- "Parae_12"

# Restrict calculations to Parae_XX chromosomes.
# This avoids including mitochondria, controls, or extra scaffolds if present.
PARAE_CHR_REGEX <- "^Parae_[0-9]{2}$"

# BED column definitions
CHR_COL     <- 1
START_COL   <- 2
END_COL     <- 3
ID_COL      <- 4
FEATURE_COL <- 5

# ------------------------------------------------------------
# Output file names
# ------------------------------------------------------------

merged_mRNA_bed <- file.path(
  OUTDIR,
  "mRNA_geneModel_coverage.merged_mRNA_intervals.bed"
)

summary_partition_tsv <- file.path(
  OUTDIR,
  "mRNA_geneModel_coverage.summary.byPartition.tsv"
)

summary_chromosome_tsv <- file.path(
  OUTDIR,
  "mRNA_geneModel_coverage.summary.byChromosome.tsv"
)

partition_pdf <- file.path(
  OUTDIR,
  "mRNA_geneModel_coverage.percentBars.stacked.byPartition.pdf"
)

chromosome_pdf <- file.path(
  OUTDIR,
  "mRNA_geneModel_coverage.percentBars.stacked.byChromosome.pdf"
)

# ------------------------------------------------------------
# Read genome index
# ------------------------------------------------------------

fai <- fread(
  FAI_FILE,
  sep = "\t",
  header = FALSE,
  fill = TRUE
)

if (ncol(fai) < 2) {
  stop("The .fai file must contain at least two columns: chromosome and length.")
}

fai <- fai[, .(
  chr = as.character(V1),
  chr_length_bp = as.numeric(V2)
)]

fai[, chr_order := .I]

# Keep only Parae_XX chromosomes
fai <- fai[grepl(PARAE_CHR_REGEX, chr)]

if (nrow(fai) == 0) {
  stop("No chromosomes matching ", PARAE_CHR_REGEX, " were found in the .fai file.")
}

if (!X_CHROM %in% fai$chr) {
  stop("The X chromosome specified in X_CHROM = ", X_CHROM, " was not found in the .fai file.")
}

fai[, chromosome_type := fifelse(
  chr == X_CHROM,
  "X_chromosome_CH12",
  "Autosome"
)]

# ------------------------------------------------------------
# Read annotation BED file
# ------------------------------------------------------------

bed <- fread(
  BED_FILE,
  sep = "\t",
  header = FALSE,
  fill = TRUE,
  quote = "",
  na.strings = c("NA", "", ".")
)

if (ncol(bed) < FEATURE_COL) {
  stop("BED file must contain at least ", FEATURE_COL, " columns.")
}

# Rename only the columns needed by the script
setnames(
  bed,
  old = c(CHR_COL, START_COL, END_COL, ID_COL, FEATURE_COL),
  new = c("chr", "start", "end", "feature_id", "feature")
)

bed[, chr := as.character(chr)]
bed[, start := as.numeric(start)]
bed[, end := as.numeric(end)]
bed[, feature := as.character(feature)]

if (any(is.na(bed$start)) || any(is.na(bed$end))) {
  stop("Some BED start/end coordinates could not be converted to numeric values.")
}

# ------------------------------------------------------------
# Filter mRNA features
# ------------------------------------------------------------

mrna <- bed[
  feature == FEATURE_NAME,
  .(chr, start, end, feature_id)
]

if (nrow(mrna) == 0) {
  stop("No rows with feature == '", FEATURE_NAME, "' were found.")
}

# Keep only chromosomes present in the filtered .fai
mrna <- mrna[chr %chin% fai$chr]

if (nrow(mrna) == 0) {
  stop("No mRNA features overlap chromosomes retained from the .fai file.")
}

# ------------------------------------------------------------
# Clean and bound coordinates to chromosome sizes
# ------------------------------------------------------------

mrna <- merge(
  mrna,
  fai[, .(chr, chr_length_bp)],
  by = "chr",
  all.x = TRUE
)

# BED coordinates should be 0-based half-open.
# These safeguards avoid intervals outside chromosome bounds.
mrna[start < 0, start := 0]
mrna[end > chr_length_bp, end := chr_length_bp]

# Remove invalid or zero-length intervals
mrna <- mrna[end > start]

if (nrow(mrna) == 0) {
  stop("All mRNA intervals were invalid or zero-length after coordinate filtering.")
}

# Remove exact duplicate intervals before merging
mrna <- unique(mrna[, .(chr, start, end)])

# ------------------------------------------------------------
# Function: merge overlapping intervals by chromosome
# ------------------------------------------------------------

merge_intervals_by_chr <- function(dt) {
  if (nrow(dt) == 0) {
    return(data.table(chr = character(), start = numeric(), end = numeric()))
  }
  
  setorder(dt, chr, start, end)
  
  merged <- dt[, {
    s <- start
    e <- end
    
    out_start <- numeric(0)
    out_end   <- numeric(0)
    
    current_start <- s[1]
    current_end   <- e[1]
    
    if (.N > 1) {
      for (i in 2:.N) {
        if (s[i] <= current_end) {
          current_end <- max(current_end, e[i])
        } else {
          out_start <- c(out_start, current_start)
          out_end   <- c(out_end, current_end)
          
          current_start <- s[i]
          current_end   <- e[i]
        }
      }
    }
    
    out_start <- c(out_start, current_start)
    out_end   <- c(out_end, current_end)
    
    .(start = out_start, end = out_end)
  }, by = chr]
  
  return(merged)
}

# ------------------------------------------------------------
# Merge mRNA intervals and calculate covered bp per chromosome
# ------------------------------------------------------------

merged_mrna <- merge_intervals_by_chr(mrna)

setorder(merged_mrna, chr, start, end)

# Save merged mRNA intervals as BED
fwrite(
  merged_mrna,
  merged_mRNA_bed,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  na = "NA"
)

cat("Wrote: ", merged_mRNA_bed, "\n", sep = "")

covered_by_chr <- merged_mrna[, .(
  mRNA_covered_bp = sum(end - start)
), by = chr]

# ------------------------------------------------------------
# Per-chromosome summary
# ------------------------------------------------------------

chr_summary <- merge(
  fai,
  covered_by_chr,
  by = "chr",
  all.x = TRUE
)

chr_summary[is.na(mRNA_covered_bp), mRNA_covered_bp := 0]

chr_summary[, non_mRNA_bp := chr_length_bp - mRNA_covered_bp]

chr_summary[, pct_mRNA_covered := 100 * mRNA_covered_bp / chr_length_bp]
chr_summary[, pct_non_mRNA := 100 * non_mRNA_bp / chr_length_bp]

setorder(chr_summary, chr_order)

chr_summary_out <- chr_summary[, .(
  chromosome = chr,
  chromosome_type,
  chromosome_length_bp = chr_length_bp,
  mRNA_covered_bp,
  non_mRNA_bp,
  pct_mRNA_covered,
  pct_non_mRNA
)]

fwrite(
  chr_summary_out,
  summary_chromosome_tsv,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote: ", summary_chromosome_tsv, "\n", sep = "")

# ------------------------------------------------------------
# Genome partition summary
# ------------------------------------------------------------

partition_list <- list(
  Whole_genome = chr_summary$chr,
  Autosomes_only = chr_summary[chromosome_type == "Autosome", chr],
  X_chromosome_only_CH12 = chr_summary[chr == X_CHROM, chr]
)

partition_summary <- rbindlist(lapply(names(partition_list), function(partition_name) {
  chrs <- partition_list[[partition_name]]
  x <- chr_summary[chr %chin% chrs]
  
  data.table(
    partition = partition_name,
    n_chromosomes = nrow(x),
    genome_length_bp = sum(x$chr_length_bp),
    mRNA_covered_bp = sum(x$mRNA_covered_bp),
    non_mRNA_bp = sum(x$non_mRNA_bp)
  )
}))

partition_summary[, pct_mRNA_covered := 100 * mRNA_covered_bp / genome_length_bp]
partition_summary[, pct_non_mRNA := 100 * non_mRNA_bp / genome_length_bp]

partition_summary[, partition := factor(
  partition,
  levels = c("Whole_genome", "Autosomes_only", "X_chromosome_only_CH12")
)]

setorder(partition_summary, partition)

fwrite(
  partition_summary,
  summary_partition_tsv,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote: ", summary_partition_tsv, "\n", sep = "")

# ------------------------------------------------------------
# Plot 1: stacked percentage bar plot by genome partition
# ------------------------------------------------------------

partition_plot_dt <- melt(
  partition_summary,
  id.vars = c("partition"),
  measure.vars = c("pct_mRNA_covered", "pct_non_mRNA"),
  variable.name = "coverage_class",
  value.name = "pct"
)

partition_plot_dt[, coverage_class := fifelse(
  coverage_class == "pct_mRNA_covered",
  "Covered by mRNA gene models",
  "Not covered by mRNA gene models"
)]

partition_plot_dt[, coverage_class := factor(
  coverage_class,
  levels = c(
    "Covered by mRNA gene models",
    "Not covered by mRNA gene models"
  )
)]

partition_plot_dt[, label := fifelse(
  pct >= 2,
  sprintf("%.1f%%", pct),
  ""
)]

p_partition <- ggplot(
  partition_plot_dt,
  aes(x = partition, y = pct, fill = coverage_class)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white"
  ) +
  scale_fill_manual(
    values = c(
      "Covered by mRNA gene models" = "darkblue",
      "Not covered by mRNA gene models" = "lightblue"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = "Genome partition",
    y = "Percentage of genome",
    fill = "Coverage class",
    title = "Genome percentage covered by mRNA gene models"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

print(p_partition)

ggsave(
  partition_pdf,
  plot = p_partition,
  width = 7,
  height = 5
)

cat("Wrote: ", partition_pdf, "\n", sep = "")

# ------------------------------------------------------------
# Plot 2: stacked percentage bar plot by chromosome
# ------------------------------------------------------------

chromosome_plot_dt <- melt(
  chr_summary_out,
  id.vars = c("chromosome", "chromosome_type"),
  measure.vars = c("pct_mRNA_covered", "pct_non_mRNA"),
  variable.name = "coverage_class",
  value.name = "pct"
)

chromosome_plot_dt[, coverage_class := fifelse(
  coverage_class == "pct_mRNA_covered",
  "Covered by mRNA gene models",
  "Not covered by mRNA gene models"
)]

chromosome_plot_dt[, coverage_class := factor(
  coverage_class,
  levels = c(
    "Covered by mRNA gene models",
    "Not covered by mRNA gene models"
  )
)]

chromosome_plot_dt[, chromosome := factor(
  chromosome,
  levels = chr_summary_out$chromosome
)]

chromosome_plot_dt[, label := fifelse(
  pct >= 5,
  sprintf("%.1f%%", pct),
  ""
)]

p_chromosome <- ggplot(
  chromosome_plot_dt,
  aes(x = chromosome, y = pct, fill = coverage_class)
) +
  geom_col(width = 0.85) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 2.4,
    color = "white"
  ) +
  scale_fill_manual(
    values = c(
      "Covered by mRNA gene models" = "darkblue",
      "Not covered by mRNA gene models" = "lightblue"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = "Chromosome",
    y = "Percentage of chromosome",
    fill = "Coverage class",
    title = "Chromosome percentage covered by mRNA gene models"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_chromosome)

ggsave(
  chromosome_pdf,
  plot = p_chromosome,
  width = 12,
  height = 5
)

cat("Wrote: ", chromosome_pdf, "\n", sep = "")

# ------------------------------------------------------------
# Final message
# ------------------------------------------------------------

cat("\nDone.\n")
cat("Summary tables and plots were written to: ", OUTDIR, "\n", sep = "")