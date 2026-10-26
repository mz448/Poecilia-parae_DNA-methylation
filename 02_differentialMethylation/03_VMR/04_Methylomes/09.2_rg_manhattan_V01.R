#!/usr/bin/env Rscript
# DATE:       2026-05-11
# AUTHOR:     MZF
# SCRIPT:     10_VMR_descriptivePlots_V02.R
# VERSION:    02
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:       Generates descriptive plots from VMRs:
#             1) Manhattan of position -vs- r.g with indicator Phenotype
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# V02: Better labels and script documentation 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# LIBRARIES
library(data.table)
library(ggplot2)


# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# SETTINGS

# infile <- "indicatorVMR.filtered.vmr.morph_and_sample_observed.RG0.5_M0.1_P0.05_MAX500.tsv"

infile <- "indicatorVMR.rgLabel.vmr.morph_and_sample_observed.filtered.diffGE_0.tsv"
fai_file <- paste0("/local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/WG.PparFemVer2024.fasta.fai")

outdir <- "Heatmaps.0.Filt-RG0.5_M0.1_P0.05_MAX500/VMR_tracks"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Load Dataframe:
dt <- fread(infile)
dt[, VMR_length := end - start]


# # Plot histogram of length. 
# h <- ggplot(dt, aes(x = VMR_length)) +
#   geom_histogram(bins = 90) +
#   labs(
#     x = "VMR length (bp)",
#     y = "Count"
#   ) +
#   theme_classic()
# print(h)
# 
# ggsave(
#   filename = paste0(outdir,"/hist_len_VMRs.pdf"),
#   plot = h,
#   width = 8,
#   height = 3
#   # ,units = "in"
# )
# 
# # Plot scatter of Lenght vs max_pairwise_diff
# s <- ggplot(dt, aes(x = VMR_length, y = max_pairwise_diff)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = TRUE) +
#   labs(
#     x = "VMR length (bp)",
#     y = "Maximum pairwise methylation difference"
#   ) +
#   theme_classic()
# print (s)
# 
# ggsave(
#   filename = paste0(outdir,"/scatt_len-vs-PWdiff_VMRs.pdf"),
#   plot = s,
#   width = 8,
#   height = 3
#   # ,units = "in"
# )
# 
# # Plot scatter of Lenght vs # of DMRs
# d <- ggplot(dt, aes(x = VMR_length, y = n_dmrs)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = TRUE) +
#   labs(
#     x = "VMR length (bp)",
#     y = "# of DMRs"
#   ) +
#   theme_classic()
# print (d)
# 
# ggsave(
#   filename = paste0(outdir,"/scatt_len-vs-Ndmr_VMRs.pdf"),
#   plot = d,
#   width = 8,
#   height = 3
#   # ,units = "in"
# )
# 
# 
# 
# # Plot scatter of max_pairwise_diff vs # of DMRs
# n <- ggplot(dt, aes(x = n_dmrs,y = max_pairwise_diff )) +
#   geom_point() +
#   geom_smooth(method = "lm", se = TRUE) +
#   labs(
#     x = "# of DMRs",
#     y = "Maximum pairwise methylation difference"
#   ) +
#   theme_classic()
# print (n)
# 
# ggsave(
#   filename = paste0(outdir,"/scatt_PWdiff-vs-Ndmr_VMRs.pdf"),
#   plot = n,
#   width = 8,
#   height = 3
#   # ,units = "in"
# )

## ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Prepare cumulative chromosome positions using real chromosome lengths

# Load chromosome lengths from .fai
fai <- fread(
  fai_file,
  header = FALSE,
  select = c(1, 2),
  col.names = c("chr", "chr_length")
)

# Keep only Parae_01 to Parae_23
fai <- fai[
  grepl("^Parae_[0-9]+$", chr)
]

fai[
  ,
  chr_number := as.integer(
    sub("^Parae_", "", chr)
  )
]

fai <- fai[
  chr_number %between% c(1, 23)
]

setorder(fai, chr_number)

expected_chrs <- sprintf("Parae_%02d", 1:23)

missing_chrs <- setdiff(expected_chrs, fai$chr)

if (length(missing_chrs) > 0) {
  stop(
    "Missing chromosomes from .fai file: ",
    paste(missing_chrs, collapse = ", ")
  )
}

# Preserve chromosome order
fai[, chr := factor(chr, levels = expected_chrs)]
setorder(fai, chr)

# Keep only VMRs on the main chromosomes
dt <- dt[chr %in% expected_chrs]

# Check required columns
required_cols <- c("chr", "start", "end", "abs_r_g", "indicator_phenotype")
missing_cols <- setdiff(required_cols, names(dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in input file: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Clean and format
dt[, start := as.numeric(start)]
dt[, end   := as.numeric(end)]
dt[, abs_r_g   := as.numeric(abs_r_g)]

dt <- dt[
  !is.na(chr) &
    !is.na(start) &
    !is.na(end) &
    !is.na(abs_r_g) &
    !is.na(indicator_phenotype)
]

dt[start > end, c("start", "end") := .(end, start)]

dt[, chr := factor(chr, levels = expected_chrs)]
setorder(dt, chr, start, end)

# Optional: preserve phenotype color order
dt[, indicator_phenotype := factor(
  indicator_phenotype,
  levels = c("female", "immaculata", "parae", "yellow")
)]

# Chromosome cumulative offsets
chr_info <- copy(fai)

chr_info[
  ,
  offset := cumsum(
    shift(chr_length, fill = 0)
  )
]

chr_info[
  ,
  center := offset + chr_length / 2
]

chr_boundaries <- chr_info[-1, offset]

# Merge offsets into VMR table
dt <- merge(
  dt,
  chr_info[, .(chr, offset)],
  by = "chr",
  all.x = TRUE,
  sort = FALSE
)

# Plot each VMR at its midpoint
dt[
  ,
  vmr_mid := start + ((end - start) / 2)
]

dt[
  ,
  genome_pos := vmr_mid + offset
]

setorder(dt, chr, start, end)

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Manhattan plot: genomic position vs abs_r_g
# Dot color = indicator phenotype


phenotype_colors <- c(
  female      = "#399957",
  immaculata  = "#888888",
  parae       = "#804e9f",
  yellow      = "#fdd800"
)


m_rg <- ggplot(
  dt,
  aes(
    x = genome_pos,
    y = abs_r_g,
    color = indicator_phenotype, 
    # shape = indicator_phenotype
  )
) +
  
  geom_vline(
    xintercept = chr_boundaries,
    linetype = "dashed",
    linewidth = 0.4,
    # color = "grey50"
    color = "red"
  ) +

  geom_point(
    # size = 1.2,
    # alpha = 0.85
    size = 2,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = phenotype_colors
  ) +
  
  scale_x_continuous(
    breaks = chr_info$center,
    labels = as.character(chr_info$chr),
    expand = expansion(mult = c(0.002, 0.002))
  ) +
  
  scale_y_continuous(
    limits = c(0.84, 1)
  ) +
  
  # scale_y_continuous(
  #   limits = c(0.5, 1)
  # ) +
  
  labs(
    x = "Chromosome",
    # y = "abs_r_g",
    y = "VMRs with strongest association to a phenotype (r.g)",
    color = "Indicator phenotype",
    # title = "Genome-wide distribution of indicator VMRs"
    title = "Genome-wide distribution of indicator VMRs - Point biserial correlation coefficient (r.g)"
  ) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(face = "bold")
  )

print(m_rg)

ggsave(
  filename = paste0(outdir, "/manhattan_rG-vs-pos_byIndicatorPhenotype.pdf"),
  plot = m_rg,
  width = 10,
  height = 4
)