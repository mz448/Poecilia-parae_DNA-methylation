#!/usr/bin/env Rscript
# DATE:       2026-03-24
# AUTHOR:     MZF & ChatGPT
# SCRIPT:     09_diffExpression_V03.R
# VERSION:    05
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   For intragenic VMRs (overlap_to_mRNA), test the relationship between
#   contrast-specific methylation difference and contrast-specific gene
#   expression divergence.
#
#   Specifically:
#     1) Load:
#          - joined intragenic DEG table
#          - VMR methylation summary table with pairwise methylation differences
#     2) Join both tables by cons_id
#     3) Match each DEG contrast to the biologically corresponding methylation
#        difference column
#     4) Build a long-format plotting table
#     5) Plot methylation difference vs |log2FC| with one linear regression per
#        contrast
#     6) Add regression equation, R², and p-value inside each facet
#
# INPUTS:
#   1) joinedIntragenicDEG:
#      Table containing intragenic VMR-gene associations and DEG columns such as:
#        i_vs_y_01_l2fc, p_vs_i_01_l2fc, y_vs_p_01_l2fc, f_vs_i_l2fc, ...
#      and a column:
#        cons_id
#
#   2) VMR methylation table:
#      Table containing VMR methylation summaries and pairwise methylation
#      difference columns such as:
#        d_F_I, d_F_P, d_F_Y, d_I_P, d_I_Y, d_P_Y
#      and a column:
#        cons_id
#
# OUTPUTS:
#   - plot database:
#       intragenic_VMR_methDiff_vs_DEG.long.tsv
#   - faceted scatter plot:
#       intragenic_VMR_methDiff_vs_absL2FC.scatter_lm_stats.pdf
#   - regression stats:
#       intragenic_VMR_methDiff_vs_absL2FC.lm_stats.tsv
#
# NOTES:
#   - The DEG metric plotted is |log2FC|
#   - The methylation metric plotted is the contrast-matched methylation
#     difference, not the maximum pairwise difference across all groups
#   - If the joined DEG file still contains non-intragenic rows, they are removed
#     using dist_signed == 0 when that column is present
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ==============================================================================
# 1) USER INPUTS
# ==============================================================================

joinedIntragenicDEG_file <- "VMR_diff0.overlap_to_mRNA.joined-DEG.tsv"
methDiff_file          <- "vmr.morph_and_sample_observed.filtered.diffGE_0.tsv"

outdir <- "../plots/intragenic_VMR_vs_DEG_in-WG"
# outdir <- "../plots/intragenic_VMR_vs_DEG_in-AUTO"
# outdir <- "../plots/intragenic_VMR_vs_DEG_in-CH12"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 2) LOAD FILES
# ==============================================================================

message("Loading joined intragenic DEG table...")
deg <- fread(joinedIntragenicDEG_file, sep = "\t", header = TRUE, quote = "")

message("Loading VMR methylation-difference table...")
meth <- fread(methDiff_file, sep = "\t", header = TRUE, quote = "")

# ==============================================================================
# 3) BASIC FILTERING / COLUMN CHECKS
# ==============================================================================

if ("dist_signed" %in% names(deg)) {
  deg <- deg[dist_signed == 0]
}

if (!"cons_id" %in% names(deg)) {
  stop("The DEG table does not contain a 'cons_id' column.")
}
if (!"cons_id" %in% names(meth)) {
  stop("The methylation table does not contain a 'cons_id' column.")
}

pair_map <- data.table(
  contrast  = c("i_vs_y_01", "p_vs_i_01", "y_vs_p_01", "f_vs_i", "f_vs_p", "f_vs_y"),
  meth_diff = c("d_I_Y",     "d_I_P",     "d_P_Y",     "d_F_I",  "d_F_P",  "d_F_Y")
)

pair_map[, l2fc_col := paste0(contrast, "_l2fc")]

missing_l2fc <- pair_map[!l2fc_col %in% names(deg), l2fc_col]
if (length(missing_l2fc) > 0) {
  stop("Missing expected DEG columns in DEG table: ",
       paste(missing_l2fc, collapse = ", "))
}

missing_meth <- pair_map[!meth_diff %in% names(meth), meth_diff]
if (length(missing_meth) > 0) {
  stop("Missing expected methylation-difference columns in methylation table: ",
       paste(missing_meth, collapse = ", "))
}

# ==============================================================================
# 3B) FILTER TO AUTOSOMES ONLY
# ==============================================================================
# 
# # In the P. parae assembly:
# #   Parae_01–Parae_23 are nuclear chromosomes
# #   Parae_12 is the sex chromosome
# # Therefore, autosomes are Parae_01–Parae_11 and Parae_13–Parae_23

# autosome_list <- paste0("Parae_", sprintf("%02d", c(1:11, 13:23)))
# 
# message("Filtering dataset to autosomes only...")
# 
# # Filter methylation-difference table using the VMR chromosome column
# if ("chr" %in% names(meth)) {
#   meth <- meth[chr %in% autosome_list]
# } else {
#   stop("Cannot filter methylation table to autosomes because column 'chr' is missing.")
# }
# 
# # Optional but useful: also filter DEG/VMR joined table if it contains DMR chromosome info
# if ("dmr_chr" %in% names(deg)) {
#   deg <- deg[dmr_chr %in% autosome_list]
# }
# 
# message("Rows in methylation table after autosome filter: ", nrow(meth))
# message("Rows in DEG table after autosome filter: ", nrow(deg))
# message("Unique autosomal cons_id in methylation table: ", uniqueN(meth$cons_id))
# message("Unique autosomal cons_id in DEG table: ", uniqueN(deg$cons_id))


# ==============================================================================
# 3B) FILTER TO CHROMOSOME 12 ONLY
# ==============================================================================
# 
# # In the P. parae assembly, Parae_12 is the sex chromosome / CH12
# ch12_id <- "Parae_12"
# 
# message("Filtering dataset to Chromosome 12 only...")
# 
# # First filter the methylation-difference table by VMR chromosome
# if ("chr" %in% names(meth)) {
#   meth <- meth[chr == ch12_id]
# } else {
#   stop("Cannot filter methylation table to CH12 because column 'chr' is missing.")
# }
# 
# # Keep only DEG/VMR-gene associations whose cons_id is present in CH12 methylation table
# ch12_cons_id <- unique(meth$cons_id)
# 
# deg <- deg[cons_id %in% ch12_cons_id]
# 
# message("Rows in methylation table after CH12 filter: ", nrow(meth))
# message("Rows in DEG table after CH12 filter: ", nrow(deg))
# message("Unique CH12 cons_id in methylation table: ", uniqueN(meth$cons_id))
# message("Unique CH12 cons_id in DEG table: ", uniqueN(deg$cons_id))

# ==============================================================================
# 4) REDUCE METHYLATION TABLE TO NEEDED COLUMNS
# ==============================================================================

meth_keep <- unique(c(
  "cons_id",
  "chr", "start", "end", "chromosome_type",
  pair_map$meth_diff
))
meth_keep <- meth_keep[meth_keep %in% names(meth)]

meth_sub <- unique(meth[, ..meth_keep], by = "cons_id")

# ==============================================================================
# 5) JOIN TABLES BY cons_id
# ==============================================================================

message("Joining DEG and methylation tables by cons_id...")
dt <- merge(deg, meth_sub, by = "cons_id", all.x = TRUE)

message("Rows after join: ", nrow(dt))
message("Unique cons_id after join: ", uniqueN(dt$cons_id))

# ==============================================================================
# 6) BUILD LONG-FORMAT PLOTTING TABLE
# ==============================================================================

message("Building long-format plotting table...")

long_list <- lapply(seq_len(nrow(pair_map)), function(i) {
  this_contrast <- pair_map$contrast[i]
  this_l2fc_col <- pair_map$l2fc_col[i]
  this_meth_col <- pair_map$meth_diff[i]
  
  tmp <- copy(dt)
  
  tmp[, contrast := this_contrast]
  tmp[, l2fc := get(this_l2fc_col)]
  tmp[, abs_l2fc := abs(l2fc)]
  tmp[, meth_diff := abs(get(this_meth_col))]
  
  keep_cols <- c(
    "cons_id",
    "gene_id", "gene_name",
    "dmr_chr", "dmr_start", "dmr_end",
    "chr", "start", "end",
    "chromosome_type",
    "contrast", "l2fc", "abs_l2fc", "meth_diff"
  )
  keep_cols <- keep_cols[keep_cols %in% names(tmp)]
  
  tmp[, ..keep_cols]
})

plot_dt <- rbindlist(long_list, use.names = TRUE, fill = TRUE)
plot_dt <- plot_dt[!is.na(meth_diff) & !is.na(abs_l2fc)]


# ==============================================================================
# 6B) BIN METHYLATION DIFFERENCE FOR BINNED BOX-AND-WHISKER PLOT
# ==============================================================================

# Bin width for methylation difference
bin_width <- 0.05

# Fixed x-axis limits for all panels
x_min <-  0
x_max <-  1

# Breaks every 0.05 from -1 to 1
meth_breaks <- seq(x_min, x_max, by = bin_width)

# Assign each methylation-difference value to a 0.05-wide bin
plot_dt[, meth_bin_id := findInterval(
  meth_diff,
  vec = meth_breaks,
  rightmost.closed = TRUE,
  all.inside = TRUE
)]

# Prevent values exactly at x_max from creating an extra bin
plot_dt[, meth_bin_id := pmin(meth_bin_id, length(meth_breaks) - 1)]

# Calculate the midpoint of each bin for plotting
plot_dt[, meth_bin_left := meth_breaks[meth_bin_id]]
plot_dt[, meth_bin_mid  := meth_bin_left + bin_width / 2]



plot_dt[, contrast := factor(
  contrast,
  levels = c("i_vs_y_01", "p_vs_i_01", "y_vs_p_01", "f_vs_i", "f_vs_p", "f_vs_y")
)]

message("Rows in plotting table: ", nrow(plot_dt))
message("Unique cons_id in plotting table: ", uniqueN(plot_dt$cons_id))

fwrite(
  plot_dt,
  file = file.path(outdir, "intragenic_VMR_methDiff_vs_DEG.long.tsv"),
  sep = "\t"
)

# ==============================================================================
# 7) REGRESSION STATS PER CONTRAST
# ==============================================================================

lm_stats <- plot_dt[, {
  fit <- lm(abs_l2fc ~ meth_diff)
  s   <- summary(fit)
  
  intercept <- unname(coef(fit)[1])
  slope     <- unname(coef(fit)[2])
  r2        <- s$r.squared
  pval      <- coef(s)[2, 4]
  
  data.table(
    n = .N,
    intercept = intercept,
    slope = slope,
    r_squared = r2,
    adj_r_squared = s$adj.r.squared,
    p_value_slope = pval
  )
}, by = contrast]

fwrite(
  lm_stats,
  file = file.path(outdir, "intragenic_VMR_methDiff_vs_absL2FC.lm_stats.tsv"),
  sep = "\t"
)

# ==============================================================================
# 8) BUILD LABELS FOR FACETS AND PANEL ANNOTATIONS
# ==============================================================================

facet_labels <- setNames(
  paste0(as.character(lm_stats$contrast), "\n(n = ", lm_stats$n, ")"),
  lm_stats$contrast
)

# format_p <- function(x) {
#   ifelse(x < 2.2e-16, "< 2.2e-16", formatC(x, format = "e", digits = 2))
# }
format_p <- function(x) {
  formatC(x, format = "f", digits = 4)
}

lm_stats[, eq_label := sprintf("y = %.3f %s %.3fx",
                               intercept,
                               ifelse(slope < 0, "-", "+"),
                               abs(slope))]

lm_stats[, r2_label := sprintf("R² = %.3f", r_squared)]
lm_stats[, p_label  := sprintf("p = %s", format_p(p_value_slope))]

lm_stats[, label := paste(eq_label, r2_label, p_label, sep = "\n")]

# ==============================================================================
# 8B) PANEL-SPECIFIC ANNOTATION POSITIONS
# ==============================================================================


#  # panel-specific top-left annotation positions
# ann_pos <- plot_dt[, .(
#   x = quantile(meth_bin_mid, probs = 0.80, na.rm = TRUE),
#   y = quantile(abs_l2fc, probs = 0.97, na.rm = TRUE)
# ), by = contrast]
# 
# ann_dt <- merge(lm_stats, ann_pos, by = "contrast", all.x = TRUE)
# 
# # Place stats labels over the actual boxplot region rather than at x = -1
# ann_pos <- plot_dt[, .(
#   x = median(meth_bin_mid, na.rm = TRUE),
#   y = quantile(abs_l2fc, probs = 0.75, na.rm = TRUE)
# ), by = contrast]
# 
# ann_dt <- merge(lm_stats, ann_pos, by = "contrast", all.x = TRUE)

# Place stats labels in the upper-right corner of each panel
ann_pos <- plot_dt[, .(
  x = Inf,
  y = Inf
), by = contrast]

ann_dt <- merge(lm_stats, ann_pos, by = "contrast", all.x = TRUE)

# # ==============================================================================
# # 9) PLOT: Scatter Plots
# # ==============================================================================

p <- ggplot(plot_dt, aes(x = meth_diff, y = abs_l2fc)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.7) +
  geom_text(
    data = ann_dt,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  facet_wrap(
    ~ contrast,
    scales = "free",
    ncol = 3,
    labeller = as_labeller(facet_labels)
  ) +
  geom_text(
    data = ann_dt,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.15,
    vjust = 1.35,
    size = 3
  ) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.25),
    minor_breaks = seq(0, 1, by = 0.05),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    breaks = seq(0, 8, by = 1),
    minor_breaks = seq(0, 1, by = 0.05),
    limits = c(0, 8)
  ) +
  labs(
    title = "Intragenic VMRs in the Whole Genome: contrast-matched methylation difference vs |log2FC|",
    # title = "Intragenic VMRs in the AUTOSOME: contrast-matched methylation difference vs |log2FC|",
    # title = "Intragenic VMRs in the X-Chromosome: contrast-matched methylation difference vs |log2FC|",
    x = "Methylation difference for the matched contrast",
    y = "|log2 fold change|"
  ) +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "black"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(outdir, "intragenic_VMR_methDiff_vs_absL2FC.scatter_lm_stats.pdf"),
  plot = p,
  width = 10,
  height = 7
)

ggsave(
  filename = file.path(outdir, "intragenic_VMR_methDiff_vs_absL2FC.scatter_lm_stats.png"),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)

# # ==============================================================================
# # 9) PLOT: BINNED BOX-AND-WHISKER PLOT
# # ==============================================================================
# 
# p <- ggplot(plot_dt, aes(y = abs_l2fc)) +
#   geom_boxplot(
#     aes(
#       x = meth_bin_mid,
#       group = meth_bin_id
#     ),
#     width = bin_width * 0.85,
#     outlier.shape = NA
#   ) +
#   geom_smooth(
#     aes(x = meth_diff, y = abs_l2fc),
#     method = "lm",
#     formula = y ~ x,
#     se = TRUE,
#     linewidth = 0.7
#   ) +
#   # geom_text(
#   #   data = ann_dt,
#   #   aes(x = x, y = y, label = label),
#   #   inherit.aes = FALSE,
#   #   hjust = 0,
#   #   vjust = 1,
#   #   size = 3
#   # ) +
#   geom_text(
#     data = ann_dt,
#     aes(x = x, y = y, label = label),
#     inherit.aes = FALSE,
#     hjust = 1.15,
#     vjust = 1.35,
#     size = 3
#   ) +
#   # facet_wrap(
#   #   ~ contrast,
#   #   scales = "free_y",
#   #   ncol = 3,
#   #   labeller = as_labeller(facet_labels)
#   # ) +
#   facet_wrap(
#     ~ contrast,
#     scales = "free_y",
#     ncol = 3,
#     labeller = as_labeller(facet_labels),
#     axes = "all_x",
#     axis.labels = "all_x"
#   )+
#   coord_cartesian(xlim = c(0, 1)) +
#   scale_x_continuous(
#     breaks = seq(0, 1, by = 0.25),
#     minor_breaks = seq(0, 1, by = 0.05),
#     limits = c(0, 1)
#   ) +
#   labs(
#     title = "Intragenic VMRs in the Whole Genome: contrast-matched methylation difference vs |log2FC|",
#     x = "Methylation difference for the matched contrast",
#     y = "|log2 fold change|"
#   ) +
#   theme_classic(base_size = 10) +
#   theme(
#     strip.background = element_rect(fill = "grey95", colour = "black"),
#     strip.text = element_text(face = "bold"),
#     plot.title = element_text(face = "bold")
#   )
# 
# ggsave(
#   filename = file.path(outdir, "intragenic_VMR_methDiff_vs_absL2FC.binnedBox_lm_stats.pdf"),
#   plot = p,
#   width = 10,
#   height = 7
# )
# 
# ggsave(
#   filename = file.path(outdir, "intragenic_VMR_methDiff_vs_absL2FC.binnedBox_lm_stats.png"),
#   plot = p,
#   width = 10,
#   height = 7,
#   dpi = 300
# )

message("Done.")
message("Outputs written to: ", outdir)
