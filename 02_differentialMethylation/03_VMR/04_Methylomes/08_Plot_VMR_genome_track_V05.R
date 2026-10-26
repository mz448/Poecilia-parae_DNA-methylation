#!/usr/bin/env Rscript
# DATE:       2026-09-18
# AUTHOR:     MZF
# SCRIPT:     09_Plot_VMR_biologicalCategory_tracks_V05.R
# VERSION:    05
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#             Plots VMRs across chromosomes using exact chromosome
#             lengths from the genome .fai file.
#
#             Version 5:
#             - Reads the biological-category per-VMR table directly
#             - Plots one genomic track per biological VMR category
#             - Preserve biological_category_order from the input table
#             - Optionally add an all-VMR reference track
#             - Generate one plot for autosomes and one plot for chromosome 12
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

# ==============================================================================
# User settings
# ==============================================================================

# Folder containing the biological-category VMR table
setwd("/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/data")

vmr_file <- "../plots/complexupset/vmr_real.contrastSharing.biologicalCategories.perVMR.tsv"
fai_file <- "/local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/WG.PparFemVer2024.fasta.fai"

outdir <- "../plots/VMR_category_tracks"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# FALSE = exactly one track per biological VMR category
# TRUE  = add the V03-style All_VMRs reference track above the category tracks
ADD_ALL_VMRS_TRACK <- TRUE

# Plot mode:
# "segments" = plot each VMR using start-end coordinates
# "ticks"    = plot each VMR as a vertical tick at its midpoint
# PLOT_MODE <- "segments"
PLOT_MODE <- "ticks"

all_track_suffix <- if (ADD_ALL_VMRS_TRACK) ".withAllTrack" else ""
output_pdf_autosome <- file.path(outdir, paste0("VMR_biologicalCategory_tracks_AUTOSOME", all_track_suffix, ".pdf"))
output_pdf_ch12 <- file.path(outdir, paste0("VMR_biologicalCategory_tracks_CH12", all_track_suffix, ".pdf"))

# ==============================================================================
# Load genome index
# ==============================================================================

fai <- fread(fai_file, header = FALSE, select = c(1, 2), col.names = c("chr", "chr_length"))
fai <- fai[grepl("^Parae_[0-9]+$", chr)]
fai[, chr_number := as.integer(sub("^Parae_", "", chr))]
setorder(fai, chr_number)

if (!"Parae_12" %in% fai$chr) stop("Parae_12 was not found in the .fai file.")

# ==============================================================================
# Load biological-category VMR table
# ==============================================================================

vmr_dt <- fread(vmr_file)

required_cols <- c("chr", "start", "end", "biological_category", "biological_category_order")
missing_cols <- setdiff(required_cols, names(vmr_dt))
if (length(missing_cols) > 0) stop("Missing required columns from VMR file: ", paste(missing_cols, collapse = ", "))

vmr_dt[, start := as.numeric(start)]
vmr_dt[, end := as.numeric(end)]
vmr_dt[, biological_category_order := as.numeric(biological_category_order)]

n_input <- nrow(vmr_dt)
vmr_dt <- vmr_dt[!is.na(chr) & !is.na(start) & !is.na(end) & !is.na(biological_category) & !is.na(biological_category_order)]
vmr_dt[start > end, c("start", "end") := .(end, start)]

if (nrow(vmr_dt) == 0) stop("No valid VMRs remain after input validation.")
if (nrow(vmr_dt) < n_input) cat("Removed", n_input - nrow(vmr_dt), "rows with missing required values.\n")

# Confirm that each biological category has one unique order value
category_order_check <- unique(vmr_dt[, .(biological_category, biological_category_order)])[, .N, by = biological_category][N > 1]
if (nrow(category_order_check) > 0) stop("At least one biological_category has multiple biological_category_order values.")

category_info <- unique(vmr_dt[, .(biological_category, biological_category_order)])
setorder(category_info, biological_category_order, biological_category)
category_levels <- category_info$biological_category

cat("Biological-category track order, top to bottom:\n")
print(category_info)

# ==============================================================================
# Helper function: add category tracks and optionally the All_VMRs track
# ==============================================================================

make_track_table <- function(dt, category_levels, add_all_track = FALSE) {
  category_tracks <- copy(dt)
  category_tracks[, track := biological_category]
  
  if (add_all_track) {
    all_track <- copy(dt)
    all_track[, track := "All_VMRs"]
    plot_dt <- rbindlist(list(category_tracks, all_track), use.names = TRUE, fill = TRUE)
    track_levels <- c(rev(category_levels), "All_VMRs")
  } else {
    plot_dt <- category_tracks
    track_levels <- rev(category_levels)
  }
  
  plot_dt[, track := factor(track, levels = track_levels)]
  list(plot_dt = plot_dt, track_levels = track_levels)
}

# ==============================================================================
# Block 1: AUTOSOME biological-category tracks
# ==============================================================================

autosome_chrs <- fai[chr != "Parae_12", chr]
auto_chr_info <- copy(fai[chr %in% autosome_chrs])
auto_chr_info[, chr_offset := cumsum(shift(chr_length, fill = 0))]
auto_chr_info[, chr_midpoint := chr_offset + chr_length / 2]
auto_boundaries <- auto_chr_info[-1, .(boundary = chr_offset)]

auto_dt <- copy(vmr_dt[chr %in% autosome_chrs])
if (nrow(auto_dt) == 0) stop("No autosomal VMRs found in the input file.")

auto_dt <- merge(auto_dt, auto_chr_info[, .(chr, chr_offset)], by = "chr", all.x = TRUE, sort = FALSE)
auto_dt[, cumulative_start := chr_offset + start]
auto_dt[, cumulative_end := chr_offset + end]
auto_dt[, cumulative_mid := cumulative_start + ((cumulative_end - cumulative_start) / 2)]

auto_tracks <- make_track_table(auto_dt, category_levels, ADD_ALL_VMRS_TRACK)
auto_plot_dt <- auto_tracks$plot_dt
auto_track_levels <- auto_tracks$track_levels

cat("\nAutosomal VMR counts by biological category:\n")
print(auto_dt[, .N, by = .(biological_category, biological_category_order)][order(biological_category_order)])

p_auto_base <- ggplot(auto_plot_dt) +
  geom_vline(
    data = auto_boundaries,
    aes(xintercept = boundary),
    linetype = "dashed",
    linewidth = 0.25,
    color = "red",
    inherit.aes = FALSE
  ) +
  scale_x_continuous(
    breaks = auto_chr_info$chr_midpoint,
    labels = auto_chr_info$chr,
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_discrete(limits = auto_track_levels, drop = FALSE) +
  labs(
    x = "Autosomes",
    y = "VMR biological category",
    title = "Autosomal VMRs by biological category"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

if (PLOT_MODE == "segments") {
  p_auto <- p_auto_base +
    geom_segment(
      aes(x = cumulative_start, xend = cumulative_end, y = track, yend = track, color = track),
      linewidth = 1.5
    )
} else if (PLOT_MODE == "ticks") {
  p_auto <- p_auto_base +
    geom_point(aes(x = cumulative_mid, y = track, color = track), shape = "|", size = 2.8)
} else {
  stop("PLOT_MODE must be either 'segments' or 'ticks'.")
}

print(p_auto)
ggsave(
  filename = output_pdf_autosome,
  plot = p_auto,
  width = 8,
  height = 1.8 + 0.38 * length(auto_track_levels),
  units = "in"
)

# ==============================================================================
# Block 2: CH12 biological-category tracks
# ==============================================================================

ch12_chr_info <- fai[chr == "Parae_12"]
if (nrow(ch12_chr_info) != 1) stop("Parae_12 was not found uniquely in the .fai file.")

ch12_dt <- copy(vmr_dt[chr == "Parae_12"])
if (nrow(ch12_dt) == 0) stop("No CH12 VMRs found after filtering to Parae_12.")

ch12_dt[, cumulative_start := start]
ch12_dt[, cumulative_end := end]
ch12_dt[, cumulative_mid := start + ((end - start) / 2)]

ch12_tracks <- make_track_table(ch12_dt, category_levels, ADD_ALL_VMRS_TRACK)
ch12_plot_dt <- ch12_tracks$plot_dt
ch12_track_levels <- ch12_tracks$track_levels

cat("\nCH12 VMR counts by biological category:\n")
print(ch12_dt[, .N, by = .(biological_category, biological_category_order)][order(biological_category_order)])

p_ch12_base <- ggplot(ch12_plot_dt) +
  scale_x_continuous(
    limits = c(0, ch12_chr_info$chr_length),
    breaks = seq(from = 0, to = ch12_chr_info$chr_length, by = 2500000),
    labels = function(x) paste0(x / 1000000, " Mb"),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  scale_y_discrete(limits = ch12_track_levels, drop = FALSE) +
  labs(
    x = "Position on chromosome X (12)",
    y = "VMR biological category",
    title = "X chromosome VMRs by biological category"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

if (PLOT_MODE == "segments") {
  p_ch12 <- p_ch12_base +
    geom_segment(
      aes(x = cumulative_start, xend = cumulative_end, y = track, yend = track, color = track),
      linewidth = 1.5
    )
} else if (PLOT_MODE == "ticks") {
  p_ch12 <- p_ch12_base +
    geom_point(aes(x = cumulative_mid, y = track, color = track), shape = "|", size = 2.8)
} else {
  stop("PLOT_MODE must be either 'segments' or 'ticks'.")
}

print(p_ch12)
ggsave(
  filename = output_pdf_ch12,
  plot = p_ch12,
  width = 8,
  height = 1.8 + 0.38 * length(ch12_track_levels),
  units = "in"
)

cat("\nSaved outputs to:\n", normalizePath(outdir), "\n", sep = "")
cat("Autosome plot:", output_pdf_autosome, "\n")
cat("CH12 plot:", output_pdf_ch12, "\n")