#!/usr/bin/env Rscript
# DATE:       2026-05-06
# AUTHOR:     MZF
# SCRIPT:     02_BarPlot_DMLs.R
# VERSION:    01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#             Plot Bar plots that use the same window size as the circos plots
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ==============================================================================
# User settings
# ==============================================================================

input_file <- "DNMTools_DML_summary/12_DNMTools_window_level_DML_summary_by_partition.tsv"

output_pdf <- "DNMTools_DML_summary/17_DNMTools_window_percentage_DMLs_manhattan_by_contrast.pdf"

# Choose one partition to avoid duplicate windows.
# Options: "WG", "AUTO", "CH12"
PLOT_PARTITION <- "WG"

# Optional: set to TRUE if you want all contrasts in one PDF page.
# FALSE = one faceted plot saved normally.
free_y_axis <- FALSE

# ==============================================================================
# Load data
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
    "Missing required columns from input file: ",
    paste(missing_cols, collapse = ", ")
  )
}

dt <- dt[partition == PLOT_PARTITION]

if (nrow(dt) == 0) {
  stop("No rows found for partition: ", PLOT_PARTITION)
}

# # ==============================================================================
# # Order chromosomes and calculate cumulative genomic coordinates
# # ==============================================================================
# 
# dt[, chr_number := suppressWarnings(as.integer(sub("^Parae_", "", chr)))]
# 
# if (any(is.na(dt$chr_number))) {
#   warning(
#     "Some chromosome names could not be converted to numeric order. ",
#     "They will be ordered alphabetically after numeric chromosomes."
#   )
# }
# 
# chr_order <- unique(dt[order(chr_number, chr), chr])
# 
# dt[, chr := factor(chr, levels = chr_order)]
# 
# chr_info <- dt[
#   ,
#   .(
#     chr_length = max(window_end, na.rm = TRUE)
#   ),
#   by = chr
# ]
# 
# # Chromosome boundary positions for vertical separator lines
# chr_boundaries <- chr_info[
#   ,
#   .(boundary = chr_offset[-1])
# ]
# 
# setorder(chr_info, chr)
# 
# chr_info[, chr_offset := cumsum(shift(chr_length, fill = 0))]
# 
# chr_info[, chr_midpoint := chr_offset + chr_length / 2]
# 
# dt <- merge(
#   dt,
#   chr_info[, .(chr, chr_offset)],
#   by = "chr",
#   all.x = TRUE
# )
# 
# dt[, window_midpoint := window_start + ((window_end - window_start) / 2)]
# 
# dt[, cumulative_position := chr_offset + window_midpoint]

# ==============================================================================
# Order chromosomes and calculate cumulative genomic coordinates
# ==============================================================================

dt[, chr_number := suppressWarnings(
  as.integer(sub("^Parae_", "", chr))
)]

if (any(is.na(dt$chr_number))) {
  warning(
    "Some chromosome names could not be converted to numeric order. ",
    "They will be ordered alphabetically after numeric chromosomes."
  )
}

# Put numeric chromosomes first, followed by nonnumeric sequences
chr_order <- unique(
  dt[
    order(is.na(chr_number), chr_number, chr),
    chr
  ]
)

dt[, chr := factor(chr, levels = chr_order)]

# Determine chromosome lengths from the largest window end
chr_info <- dt[
  ,
  .(
    chr_length = max(window_end, na.rm = TRUE)
  ),
  by = chr
]

# Ensure chromosomes follow the factor-level order
setorder(chr_info, chr)

# Starting cumulative position of each chromosome
chr_info[
  ,
  chr_offset := cumsum(shift(chr_length, fill = 0))
]

# Position used for chromosome labels
chr_info[
  ,
  chr_midpoint := chr_offset + chr_length / 2
]

# Boundaries are the starting positions of chromosomes 2 through n
chr_boundaries <- chr_info[
  -1,
  .(boundary = chr_offset)
]

# Add chromosome offsets to the window-level dataset
dt <- merge(
  dt,
  chr_info[, .(chr, chr_offset)],
  by = "chr",
  all.x = TRUE,
  sort = FALSE
)

dt[
  ,
  window_midpoint :=
    window_start + ((window_end - window_start) / 2)
]

dt[
  ,
  cumulative_position :=
    chr_offset + window_midpoint
]



# ==============================================================================
# Optional contrast ordering
# ==============================================================================

contrast_order <- c(
  "f-vs-i",
  "f-vs-p",
  "f-vs-y",
  "p-vs-i",
  "p-vs-y",
  "y-vs-i"
)

dt[, contrast := factor(
  contrast,
  levels = contrast_order[contrast_order %in% unique(dt$contrast)]
)]

# ==============================================================================
# Manhattan-style bar plot
# ==============================================================================

p <- ggplot(
  dt,
  aes(
    x = cumulative_position,
    y = percentage_DMLs
  )
) +
  geom_col(
    width = median(dt$window_end - dt$window_start, na.rm = TRUE) * 0.85
  ) +
  # geom_point(size = 0.4, alpha = 0.7)+ # Manhattan plot style
  geom_vline(
    data = chr_boundaries,
    aes(xintercept = boundary),
    linetype = "dashed",
    linewidth = 0.25,
    color = "red",
    inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ contrast,
    ncol = 1,
    scales = ifelse(free_y_axis, "free_y", "fixed")
  ) +
  scale_x_continuous(
    breaks = chr_info$chr_midpoint,
    labels = as.character(chr_info$chr),
    expand = expansion(mult = c(0.005, 0.005))
  ) +
  labs(
    x = "Chromosome",
    y = "DMLs (% of covered CpGs per window)",
    title = paste0(
      "Window-level percentage of DMLs across the genome: ",
      PLOT_PARTITION
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

print(p)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 14,
  height = 10,
  units = "in"
)

