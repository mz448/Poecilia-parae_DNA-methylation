#!/usr/bin/env Rscript
# DATE:       2026-05-11
# AUTHOR:     MZF
# SCRIPT:     09_Plot_VMR_cluster_tracks_V03.R
# VERSION:    03
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#             Plot consensus DMRs across chromosomes using exact chromosome
#             lengths from the genome .fai file.
#
#             Version 3:
#             - Plot one track per VMR cluster
#             - Add one additional reference track containing all VMRs
#               irrespective of cluster
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
# Go to the plot's folder of the thresholded VMRs
setwd("/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/plots/0.4")

fai_file <- "/local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/WG.PparFemVer2024.fasta.fai"
outdir      <- "VMR_tracks"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

autosome_cluster_file <- "clusters_AUTOSOME.sampleHeatmap.diffGE_0.4.tsv"
ch12_cluster_file     <- "clusters_CH12.sampleHeatmap.diffGE_0.4.tsv"

output_pdf_autosome <- file.path(outdir,"VMR_cluster_tracks_AUTOSOME.diffGE_0.4.withAllTrack.pdf")
output_pdf_ch12     <- file.path(outdir,"VMR_cluster_tracks_CH12.diffGE_0.4.withAllTrack.pdf")

# Plot mode:
# "segments" = plot each VMR using start-end coordinates
# "ticks"    = plot each VMR as a vertical tick at its midpoint
# PLOT_MODE <- "segments"
PLOT_MODE <- "ticks"

# ==============================================================================
# Load genome index
# ==============================================================================

fai <- fread(
  fai_file,
  header = FALSE,
  select = c(1, 2),
  col.names = c("chr", "chr_length")
)

fai <- fai[grepl("^Parae_[0-9]+$", chr)]

fai[, chr_number := as.integer(sub("^Parae_", "", chr))]
setorder(fai, chr_number)

# ==============================================================================
# Block 1: AUTOSOME VMR cluster tracks + all-VMR reference track
# ==============================================================================

auto_dt <- fread(autosome_cluster_file)

required_cols <- c("chr", "start", "end", "Cluster")
missing_cols <- setdiff(required_cols, names(auto_dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns from autosome cluster file: ",
    paste(missing_cols, collapse = ", ")
  )
}

auto_dt[, start := as.numeric(start)]
auto_dt[, end   := as.numeric(end)]

auto_dt <- auto_dt[
  !is.na(chr) &
    !is.na(start) &
    !is.na(end) &
    !is.na(Cluster)
]

auto_dt[start > end, c("start", "end") := .(end, start)]

# Keep autosomes only: Parae_01 to Parae_23, excluding Parae_12
autosome_chrs <- fai[chr != "Parae_12", chr]

auto_chr_info <- fai[chr %in% autosome_chrs]

auto_chr_info[, chr := factor(chr, levels = autosome_chrs)]
setorder(auto_chr_info, chr)

auto_chr_info[, chr_offset := cumsum(shift(chr_length, fill = 0))]
auto_chr_info[, chr_midpoint := chr_offset + chr_length / 2]

auto_boundaries <- auto_chr_info[
  -1,
  .(
    boundary = chr_offset
  )
]

auto_dt <- auto_dt[chr %in% autosome_chrs]
auto_dt[, chr := factor(chr, levels = autosome_chrs)]

auto_dt <- merge(
  auto_dt,
  auto_chr_info[, .(chr, chr_offset)],
  by = "chr",
  all.x = TRUE
)

auto_dt[, cumulative_start := chr_offset + start]
auto_dt[, cumulative_end   := chr_offset + end]
auto_dt[, cumulative_mid   := cumulative_start + ((cumulative_end - cumulative_start) / 2)]

# ------------------------------------------------------------------------------
# Add cluster-specific tracks plus an all-VMR reference track
# ------------------------------------------------------------------------------

auto_dt[, Cluster := as.character(Cluster)]

auto_cluster_levels <- sort(unique(as.integer(auto_dt$Cluster)))
auto_cluster_levels <- as.character(auto_cluster_levels)

auto_cluster_tracks <- copy(auto_dt)
auto_cluster_tracks[, track := paste0("Cluster_", Cluster)]

auto_all_track <- copy(auto_dt)
auto_all_track[, track := "All_VMRs"]

auto_plot_dt <- rbindlist(
  list(
    auto_cluster_tracks,
    auto_all_track
  ),
  use.names = TRUE,
  fill = TRUE
)

# This order puts All_VMRs as the top track in the plot
auto_track_levels <- c(
  paste0("Cluster_", auto_cluster_levels),
  "All_VMRs"
)

auto_plot_dt[, track := factor(track, levels = auto_track_levels)]

# ------------------------------------------------------------------------------
# Plot autosomes
# ------------------------------------------------------------------------------

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
    labels = as.character(auto_chr_info$chr),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(
    x = "Autosomes",
    y = "VMR cluster",
    title = "Autosomal VMRs by methylation-profile cluster"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "none"
  )

if (PLOT_MODE == "segments") {
  
  p_auto <- p_auto_base +
    geom_segment(
      aes(
        x = cumulative_start,
        xend = cumulative_end,
        y = track,
        yend = track,
        color = track
      ),
      linewidth = 1.5
      # ,alpha = 0.85
    )
  
} else if (PLOT_MODE == "ticks") {
  
  p_auto <- p_auto_base +
    geom_point(
      aes(
        x = cumulative_mid,
        y = track,
        color = track
      ),
      shape = "|",
      size = 2.8
      # ,alpha = 0.85
    )
  
} else {
  stop("PLOT_MODE must be either 'segments' or 'ticks'.")
}

print(p_auto)

ggsave(
  filename = output_pdf_autosome,
  plot = p_auto,
  width = 8,
  height = 3
  # ,units = "in"
)

# ==============================================================================
# Block 2: CH12 VMR cluster tracks + all-VMR reference track
# ==============================================================================

ch12_dt <- fread(ch12_cluster_file)

required_cols <- c("chr", "start", "end", "Cluster")
missing_cols <- setdiff(required_cols, names(ch12_dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns from CH12 cluster file: ",
    paste(missing_cols, collapse = ", ")
  )
}

ch12_dt[, start := as.numeric(start)]
ch12_dt[, end   := as.numeric(end)]

ch12_dt <- ch12_dt[
  !is.na(chr) &
    !is.na(start) &
    !is.na(end) &
    !is.na(Cluster)
]

ch12_dt[start > end, c("start", "end") := .(end, start)]

# Keep chromosome 12 only
ch12_chr_info <- fai[chr == "Parae_12"]

if (nrow(ch12_chr_info) != 1) {
  stop("Parae_12 was not found uniquely in the .fai file.")
}

ch12_dt <- ch12_dt[chr == "Parae_12"]

if (nrow(ch12_dt) == 0) {
  stop("No CH12 VMRs found after filtering to Parae_12.")
}

ch12_dt[, cumulative_start := start]
ch12_dt[, cumulative_end   := end]
ch12_dt[, cumulative_mid   := start + ((end - start) / 2)]

# ------------------------------------------------------------------------------
# Add cluster-specific tracks plus an all-VMR reference track
# ------------------------------------------------------------------------------

ch12_dt[, Cluster := as.character(Cluster)]

ch12_cluster_levels <- sort(unique(as.integer(ch12_dt$Cluster)))
ch12_cluster_levels <- as.character(ch12_cluster_levels)

ch12_cluster_tracks <- copy(ch12_dt)
ch12_cluster_tracks[, track := paste0("Cluster_", Cluster)]

ch12_all_track <- copy(ch12_dt)
ch12_all_track[, track := "All_VMRs"]

ch12_plot_dt <- rbindlist(
  list(
    ch12_cluster_tracks,
    ch12_all_track
  ),
  use.names = TRUE,
  fill = TRUE
)

# This order puts All_VMRs as the top track in the plot
ch12_track_levels <- c(
  paste0("Cluster_", ch12_cluster_levels),
  "All_VMRs"
)

ch12_plot_dt[, track := factor(track, levels = ch12_track_levels)]

# ------------------------------------------------------------------------------
# Plot CH12
# ------------------------------------------------------------------------------

p_ch12_base <- ggplot(ch12_plot_dt) +
  # scale_x_continuous(
  #   limits = c(0, ch12_chr_info$chr_length),
  #   expand = expansion(mult = c(0.005, 0.005))
  # scale_x_continuous(
  #   limits = c(0, ch12_chr_info$chr_length),
  #   breaks = seq(
  #     from = 0,
  #     to = ch12_chr_info$chr_length,
  #     by = 2500000
  #   ),
  #   labels = scales::comma,
  #   expand = expansion(mult = c(0.005, 0.005))
  # ) +
  scale_x_continuous(
    limits = c(0, ch12_chr_info$chr_length),
    breaks = seq(
      from = 0,
      to = ch12_chr_info$chr_length,
      by = 2500000
    ),
    labels = function(x) paste0(x / 1000000, " Mb"),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(
    x = "Position on chromosome X (12)",
    y = "VMR cluster",
    title = "X chromosome VMRs by methylation-profile cluster"
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
      aes(
        x = cumulative_start,
        xend = cumulative_end,
        y = track,
        yend = track,
        color = track
      ),
      linewidth = 1.5,
      # ,alpha = 0.85
    )
  
} else if (PLOT_MODE == "ticks") {
  
  p_ch12 <- p_ch12_base +
    geom_point(
      aes(
        x = cumulative_mid,
        y = track,
        color = track
      ),
      shape = "|",
      size = 2.8
      # ,alpha = 0.85
    )
  
} else {
  stop("PLOT_MODE must be either 'segments' or 'ticks'.")
}

print(p_ch12)

ggsave(
  filename = output_pdf_ch12,
  plot = p_ch12,
  width = 8,
  height = 2.5
  # ,units = "in"
)
