#!/usr/bin/env Rscript

# DATE: 20260315
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Visualize and summarize the genomic position of VMR
#   relative to the closest mRNA feature.
#
#   The script generates:
#     1) Histogram of signed distance between VMR and the closest mRNA
#     2) Scatter-style histogram of counts by signed distance
#     3) Percentage bar plot for three positional classes:
#          dist_signed < 0  → upstream of gene
#          dist_signed == 0 → overlap with gene body
#          dist_signed > 0  → downstream of gene
#     4) Stacked percentage bar plot by chromosome type
#          - Autosome
#          - Sex_Ch
#
#   The script also exports lists of VMR by genomic relationship to genes:
#     • VMR overlapping a gene (dist_signed == 0)
#     • intergenic VMR (dist_signed != 0)
#
#   Finally, it generates a summary table with counts and percentages of
#   VMR by chromosome type and positional category.
#
# INPUT:
#   TSV file produced by the closest-mRNA annotation script:
#     VMR_diff<DIFF_THRESH>.closest_mRNA.tsv
#
#   Required columns:
#     dist_signed
#     chromosome_type
#
# OUTPUT FILES:
#   Figures
#     - VMR_diff<DIFF_THRESH>.distTo-mRNA.histogram.pdf
#     - VMR_diff<DIFF_THRESH>.distToTSS.scatterHist.pdf
#     - VMR_diff<DIFF_THRESH>.distTo-mRNA.direction_percentBars.pdf
#     - VMR_diff<DIFF_THRESH>.distTo-mRNA.direction_percentBars.stacked.byChromType.pdf
#
#   Tables
#     - VMR_diff<DIFF_THRESH>.overlap_to_mRNA.tsv
#     - VMR_diff<DIFF_THRESH>.intergenic.tsv
#     - VMR_diff<DIFF_THRESH>.overlap_summary.byChromType.tsv
#
# DEFINITIONS:
#   dist_signed interpretation:
#     negative = VMR upstream of gene
#     zero     = VMR overlaps gene body
#     positive = VMR downstream of gene
#
#   genomic classes:
#     upstream   → dist_signed < 0
#     overlap    → dist_signed == 0
#     downstream → dist_signed > 0
#
#   chromosome_type:
#     Autosome → Parae_01–Parae_23 except Parae_12
#     Sex_Ch   → Parae_12
#
# DEPENDENCIES:
#   data.table
#   ggplot2
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggbreak) 
})

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
DIFF_THRESH <- "0"
plots_folder <- "../plots/"
outdir<- "dist_VMR-to-gene"
outFolder <-  paste0( "../plots/","dist_VMR-to-gene")

# infile <- paste0("../plots/", DIFF_THRESH,"/VMR_diff", DIFF_THRESH, ".closest_mRNA.tsv")
infile <- "/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/plots/0/VMRs_diff0.closest_mRNA.tsv"
# load dataframe
dt <- fread(infile, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))

dir.create(outFolder, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# Simple histogram of signed distance to TSS
# ------------------------------------------------------------
hist_pdf <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH, "/",outdir),
  paste0("VMR_diff", DIFF_THRESH, ".distTo-mRNA.histogram.pdf")
)

p_hist <- ggplot(dt[dist_signed >= -30000 & dist_signed <= 30000], aes(x = dist_signed)) +
  geom_histogram(bins = 100, color="black", fill="grey") +
  # geom_density(color="purple")+
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Distance",
    y = "Count of VMR",
    title = "Histogram of consensus DMR distance to closest mRNA"
      ) +
  # scale_x_continuous(limits = c(-30000, 30000)) +
  # coord_cartesian(xlim = c(-500, 500))+ # Set x-axis limits from 0 to 10000
  # scale_y_break(c(300, 4900))+
  theme_classic(base_size = 10)


print (p_hist)

ggsave(hist_pdf, plot = p_hist, width = 8, height = 5)

cat("Wrote: ", hist_pdf, "\n", sep = "")


# ------------------------------------------------------------
# Percentage bar plot for upstream / overlap / downstream
# ------------------------------------------------------------
bar_pdf <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH, "/",outdir),
  paste0("VMR_diff", DIFF_THRESH, ".distTo-mRNA.direction_percentBars.pdf")
)

dt_bar <- copy(dt)
dt_bar <- dt_bar[!is.na(dist_signed)]

dt_bar[, distance_class := fifelse(
  dist_signed < 0, "Upstream (<0)",
  fifelse(dist_signed == 0, "Overlap (=0)", "Downstream (>0)")
)]

dt_bar[, distance_class := factor(
  distance_class,
  levels = c("Upstream (<0)", "Overlap (=0)", "Downstream (>0)")
)]

bar_df <- dt_bar[, .N, by = distance_class]
bar_df[, pct := 100 * N / sum(N)]
bar_df[, label := sprintf("%.2f%%", pct)]

p_bar <- ggplot(bar_df, aes(x = distance_class, y = pct)) +
  geom_col() +
  geom_text(aes(label = label), vjust = -0.4, size = 4) +
  labs(
    x = "Distance class relative to closest mRNA",
    y = "Percentage of VMR",
    title = "Consensus DMR positions relative to closest mRNA"
  ) +
  ylim(0, max(bar_df$pct) * 1.12) +
  theme_classic(base_size = 10)

print(p_bar)

ggsave(bar_pdf, plot = p_bar, width = 6, height = 5)

cat("Wrote: ", bar_pdf, "\n", sep = "")


# ------------------------------------------------------------
# Percentage bar plot for upstream / overlap / downstream by chromosome type
# ------------------------------------------------------------
bar_pdf <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH, "/",outdir),
  paste0("VMR_diff", DIFF_THRESH, ".distTo-mRNA.direction_percentBars.byChromType.pdf")
)

dt_bar <- copy(dt)
dt_bar <- dt_bar[!is.na(dist_signed)]

dt_bar[, distance_class := fifelse(
  dist_signed < 0, "Upstream (<0)",
  fifelse(dist_signed == 0, "Overlap (=0)", "Downstream (>0)")
)]

dt_bar[, distance_class := factor(
  distance_class,
  levels = c("Upstream (<0)", "Overlap (=0)", "Downstream (>0)")
)]

bar_df <- dt_bar[, .N, by = .(chromosome_type, distance_class)]
bar_df[, pct := 100 * N / sum(N), by = chromosome_type]
bar_df[, label := sprintf("%.2f%%", pct)]

p_bar <- ggplot(bar_df, aes(x = distance_class, y = pct)) +
  geom_col() +
  geom_text(aes(label = label), vjust = -0.4, size = 3) +
  facet_wrap(~ chromosome_type) +
  labs(
    x = "Distance class relative to closest mRNA",
    y = "Percentage of VMR",
    title = "Consensus DMR positions relative to closest mRNA"
  ) +
  theme_classic(base_size = 10)

print(p_bar)

ggsave(bar_pdf, plot = p_bar, width = 8, height = 5)

cat("Wrote: ", bar_pdf, "\n", sep = "")


# ------------------------------------------------------------
# Percentage stacked bar plot for upstream / overlap / downstream
# by chromosome type
#   - one stacked bar for Autosome
#   - one stacked bar for Sex_Ch
#   - colors:
#       upstream   = dark blue
#       overlap    = green
#       downstream = light blue
# ------------------------------------------------------------
bar_pdf <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH, "/",outdir),
  paste0("VMR_diff", DIFF_THRESH, ".distTo-mRNA.direction_percentBars.stacked.byChromType.pdf")
)

dt_bar <- copy(dt)
dt_bar <- dt_bar[!is.na(dist_signed) & !is.na(chromosome_type)]

dt_bar[, distance_class := fifelse(
  dist_signed < 0, "Upstream (<0)",
  fifelse(dist_signed == 0, "Overlap (=0)", "Downstream (>0)")
)]

dt_bar[, distance_class := factor(
  distance_class,
  levels = c("Upstream (<0)", "Overlap (=0)", "Downstream (>0)")
)]

dt_bar[, chromosome_type := factor(
  chromosome_type,
  levels = c("Autosome", "Sex_Ch")
)]

bar_df <- dt_bar[, .N, by = .(chromosome_type, distance_class)]
bar_df[, pct := 100 * N / sum(N), by = chromosome_type]

p_bar <- ggplot(bar_df, aes(x = chromosome_type, y = pct, fill = distance_class)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct)),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "white"
  ) +
  scale_fill_manual(
    values = c(
      "Upstream (<0)"   = "darkblue",
      "Overlap (=0)"    = "green",
      "Downstream (>0)" = "lightblue"
    )
  ) +
  labs(
    x = "Chromosome type",
    y = "Percentage of VMR",
    fill = "Distance class",
    title = "Consensus DMR positions relative to closest mRNA"
  ) +
  theme_classic(base_size = 10)

print(p_bar)

ggsave(bar_pdf, plot = p_bar, width = 6, height = 5)

cat("Wrote: ", bar_pdf, "\n", sep = "")


# scatter histogram -----------------------------------------------------------      
# Optional bin size in bp.
# Use 1 for exact-distance counts.
BIN_SIZE <- 1000
#
# Optional x-axis limits (set to NULL to keep full range)
XMIN <- NULL
XMAX <- NULL
#
# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------
dt <- fread(infile, sep = "\t", header = TRUE, na.strings = c("NA", "", "."))

if (!"dist_signed" %in% names(dt)) {
  stop("Missing required column: dist_signed")
}

dt <- dt[!is.na(dist_signed)]
dt[, dist_signed := as.numeric(dist_signed)]

# ------------------------------------------------------------
# Build count table
# ------------------------------------------------------------
if (BIN_SIZE == 1) {
  counts <- dt[, .N, by = .(x = dist_signed)]
  xlab_txt <- "Signed distance to TSS (bp)"
} else {
  dt[, x := floor(dist_signed / BIN_SIZE) * BIN_SIZE]
  counts <- dt[, .N, by = x]
  xlab_txt <- paste0("Signed distance to TSS (bp; binned by ", BIN_SIZE, " bp)")
}

setorder(counts, x)

# ------------------------------------------------------------
# Optional x-axis filtering
# ------------------------------------------------------------
if (!is.null(XMIN)) counts <- counts[x >= XMIN]
if (!is.null(XMAX)) counts <- counts[x <= XMAX]
#
# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

out_pdf <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH),
  paste0("VMR_diff", DIFF_THRESH, ".distTomRNA.scatterHist.pdf")
  # paste0("VMR_diff", DIFF_THRESH, ".distToTSS.scatterHist_zoom.pdf")
)


p <- ggplot(counts, aes(x = x, y = N)) +
  geom_point(size = .5) +
  geom_vline(xintercept = 0, linetype = "dashed", color="red") +
  labs(
    x = xlab_txt,
    y = "Count of VMR",
    title = "VMR occupancy relative to closest mRNA"
  ) +
  coord_cartesian(xlim = c(-10000, 10000))+ # Set x-axis limits from -10000 to 10000
  theme_classic(base_size = 12)
print (p)

ggsave(out_pdf, plot = p, width = 8, height = 5)
cat("Wrote: ", out_pdf, "\n", sep = "")

# ------------------------------------------------------------
# Export overlap / non-overlap sVMR lists + summary by chromosome type
# Definitions:
#   overlap_to_mRNA = dist_signed == 0
#   intergenic      = dist_signed != 0
#   upstream        = dist_signed < 0
#   downstream      = dist_signed > 0
# ------------------------------------------------------------

overlap_tsv <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH, "/",outdir),
  paste0("VMR_diff", DIFF_THRESH, ".overlap_to_mRNA.tsv")
)

intergenic_tsv <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH, "/",outdir),
  paste0("VMR_diff", DIFF_THRESH, ".intergenic.tsv")
)

summary_tsv <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH, "/",outdir),
  paste0("VMR_diff", DIFF_THRESH, ".overlap_summary.byChromType.tsv")
)

dt2 <- copy(dt)
dt2 <- dt2[!is.na(dist_signed) & !is.na(chromosome_type)]

# overlap list
dt_overlap <- dt2[dist_signed == 0]
fwrite(dt_overlap, overlap_tsv, sep = "\t", quote = FALSE, na = "NA")

# non-overlap / intergenic list
dt_intergenic <- dt2[dist_signed != 0]
fwrite(dt_intergenic, intergenic_tsv, sep = "\t", quote = FALSE, na = "NA")

# summary
summary_dt <- dt2[, .(
  total_VMR          = .N,
  count_overlap_to_mRNA = sum(dist_signed == 0),
  count_intergenic      = sum(dist_signed != 0),
  count_overlap_upstream   = sum(dist_signed < 0),
  count_overlap_downstream = sum(dist_signed > 0)
), by = chromosome_type]

summary_dt[, pct_overlap_to_mRNA := 100 * count_overlap_to_mRNA / total_VMR]
summary_dt[, pct_intergenic      := 100 * count_intergenic / total_VMR]
summary_dt[, pct_upstream        := 100 * count_overlap_upstream / total_VMR]
summary_dt[, pct_downstream      := 100 * count_overlap_downstream / total_VMR]

setcolorder(
  summary_dt,
  c(
    "chromosome_type",
    "total_VMR",
    "count_overlap_to_mRNA",
    "count_intergenic",
    "count_overlap_upstream",
    "count_overlap_downstream",
    "pct_overlap_to_mRNA",
    "pct_intergenic",
    "pct_upstream",
    "pct_downstream"
  )
)

fwrite(summary_dt, summary_tsv, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote: ", overlap_tsv, "\n", sep = "")
cat("Wrote: ", intergenic_tsv, "\n", sep = "")
cat("Wrote: ", summary_tsv, "\n", sep = "")

# ------------------------------------------------------------
# Per-chromosome summary of VMR / VMR occupancy relative to mRNA
#
# This block is independent of Script 03.
# It summarizes, per chromosome:
#   - total number of VMR
#   - number and percentage overlapping mRNA
#   - number and percentage intergenic
#   - number and percentage upstream
#   - number and percentage downstream
#
# Definitions:
#   overlap_to_mRNA = dist_signed == 0
#   intergenic      = dist_signed != 0
#   upstream        = dist_signed < 0
#   downstream      = dist_signed > 0
#
# Output:
#   VMR_diff<DIFF_THRESH>.overlap_summary.byChromosome.tsv
# ------------------------------------------------------------

summary_by_chr_tsv <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH, "/",outdir),
  paste0("VMR_diff", DIFF_THRESH, ".overlap_summary.byChromosome.tsv")
)

dt_chr <- copy(dt2)
dt_chr <- dt_chr[!is.na(dist_signed)]

# ------------------------------------------------------------
# Identify chromosome column
# ------------------------------------------------------------

chr_col_candidates <- c(
  "chromosome",
  "chr",
  "seqnames",
  "seqname",
  "VMR_chr",
  "DMR_chr",
  "V1"
)


chr_col <- "gene_chr"
dt_chr[, chromosome := as.character(get(chr_col))]



# Use Script 03-compatible chromosome type labels
dt_chr[, chromosome_type := fifelse(
  chromosome == "Parae_12",
  "X_chromosome_CH12",
  "Autosome"
)]

# ------------------------------------------------------------
# Summarize VMR / VMRs per chromosome
# ------------------------------------------------------------

summary_by_chr <- dt_chr[, .(
  total_VMR = .N,
  count_overlap_to_mRNA = sum(dist_signed == 0),
  count_intergenic = sum(dist_signed != 0),
  count_overlap_upstream = sum(dist_signed < 0),
  count_overlap_downstream = sum(dist_signed > 0)
), by = .(chromosome, chromosome_type)]

summary_by_chr[, pct_VMR_overlap_to_mRNA := 100 * count_overlap_to_mRNA / total_VMR]
summary_by_chr[, pct_VMR_intergenic := 100 * count_intergenic / total_VMR]
summary_by_chr[, pct_VMR_upstream := 100 * count_overlap_upstream / total_VMR]
summary_by_chr[, pct_VMR_downstream := 100 * count_overlap_downstream / total_VMR]

# Natural chromosome order: Parae_01, Parae_02, ..., Parae_23
summary_by_chr[, chr_num := as.integer(gsub("Parae_", "", chromosome))]
setorder(summary_by_chr, chr_num)

summary_by_chr_out <- summary_by_chr[, .(
  chromosome,
  chromosome_type,
  total_VMR,
  count_overlap_to_mRNA,
  count_intergenic,
  count_overlap_upstream,
  count_overlap_downstream,
  pct_VMR_overlap_to_mRNA,
  pct_VMR_intergenic,
  pct_VMR_upstream,
  pct_VMR_downstream
)]

fwrite(
  summary_by_chr_out,
  summary_by_chr_tsv,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("Wrote: ", summary_by_chr_tsv, "\n", sep = "")


# ------------------------------------------------------------
# Stacked percentage bar plot by chromosome
# VMR / VMRs overlapping mRNA vs intergenic
#
# Output:
#   VMR_diff<DIFF_THRESH>.distTo-mRNA.direction_percentBars.stacked.byChromosome.pdf
# ------------------------------------------------------------

bar_by_chr_pdf <- file.path(
  "..",
  paste0("plots/", DIFF_THRESH, "/",outdir),
  paste0(
    "VMR_diff",
    DIFF_THRESH,
    ".distTo-mRNA.direction_percentBars.stacked.byChromosome.pdf"
  )
)

VMR_chr_plot_dt <- melt(
  summary_by_chr_out,
  id.vars = c("chromosome", "chromosome_type"),
  measure.vars = c(
    "pct_VMR_overlap_to_mRNA",
    "pct_VMR_intergenic"
  ),
  variable.name = "coverage_class",
  value.name = "pct"
)

VMR_chr_plot_dt[, coverage_class := fifelse(
  coverage_class == "pct_VMR_overlap_to_mRNA",
  "VMR overlapping mRNA",
  "Intergenic VMR"
)]

VMR_chr_plot_dt[, coverage_class := factor(
  coverage_class,
  levels = c(
    "VMR overlapping mRNA",
    "Intergenic VMR"
  )
)]

VMR_chr_plot_dt[, chromosome := factor(
  chromosome,
  levels = summary_by_chr_out$chromosome
)]

VMR_chr_plot_dt[, label := fifelse(
  !is.na(pct) & pct >= 5,
  sprintf("%.1f%%", pct),
  ""
)]

p_VMR_chr <- ggplot(
  VMR_chr_plot_dt,
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
      "VMR overlapping mRNA" = "darkblue",
      "Intergenic VMR" = "lightblue"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    x = "Chromosome",
    y = "Percentage of VMR",
    fill = "VMR class",
    title = "Chromosome-level VMR occupancy relative to genes"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_VMR_chr)

ggsave(
  bar_by_chr_pdf,
  plot = p_VMR_chr,
  width = 12,
  height = 5
)

cat("Wrote: ", bar_by_chr_pdf, "\n", sep = "")
