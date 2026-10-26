#!/usr/bin/env Rscript
# DATE:       2026-09-18
# AUTHOR:     MZF
# SCRIPT:     9_1_VMR_descriptivePlots_V05.R
# VERSION:    05
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:       Generates descriptive plots from VMRs:
#             1) Histogram of VMR lengths
#             2) scatter of Length vs max_pairwise_diff
#             3) scatter of Length vs # of DMRs
#             4) scatter of max_pairwise_diff vs # of DMRs
#             5) Manhattan of position and PWdiff
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# V02: Better labels and script documentation 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# LIBRARIES
library(data.table)
library(ggplot2)

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# SETTINGS

# infile <- "vmr.morph_and_sample_observed.filtered.diffGE_0.tsv"
# outdir <- "VMR_tracks"

infile <- "../plots/0/vmr.morph_and_sample_observed.filtered.diffGE_0.tsv"
fai_file <- paste0("/local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/WG.PparFemVer2024.fasta.fai")

# outdir <- "Heatmaps.0.Filt-RG0.5_M0.1_P0.05_MAX500/VMR_tracks"
outdir <- "../plots/VMR_descriptivePlots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Load Dataframe:
dt <- fread(infile)
dt[, VMR_length := end - start]


# Plot histogram of length.  ---------------------------------------------------
h <- ggplot(dt, aes(x = VMR_length)) +
  geom_histogram(bins = 90, fill="darkgreen"
                 # ,color="grey"
                 ) +
  labs(
    x = "VMR length (bp)",
    y = "VMR count"
  ) +
  theme_classic()
print(h)

ggsave(
  filename = paste0(outdir,"/hist_len_VMRs.pdf"),
  plot = h,
  width = 8,
  height = 3
  # ,units = "in"
)

# {By chromosome}  Plot histogram of length.  ---------------------------------------------------
h <- ggplot(dt, aes(x = VMR_length, color=chromosome_type)) +
  geom_histogram(bins = 90) +
  labs(
    x = "VMR length (bp)",
    y = "Count"
  ) +
  theme_classic()
print(h)

ggsave(
  filename = paste0(outdir,"/hist_len_VMRs_chr.pdf"),
  plot = h,
  width = 8,
  height = 3
  # ,units = "in"
)

# Plot scatter of VMR Lenght vs max_pairwise_diff  ---------------------------------
s <- ggplot(dt, aes(x = VMR_length, y = max_pairwise_diff)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "VMR length (bp)",
    y = "Maximum pairwise methylation difference"
  ) +
  theme_classic()
print (s)

ggsave(
  filename = paste0(outdir,"/scatt_len-vs-PWdiff_VMRs.pdf"),
  plot = s,
  width = 8,
  height = 3
  # ,units = "in"
)


# {By chromosome} Plot scatter of VMR Lenght vs max_pairwise_diff  ---------------------------------
s <- ggplot(dt, aes(x = VMR_length, y = max_pairwise_diff, color=chromosome_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "VMR length (bp)",
    y = "Maximum pairwise methylation difference"
  ) +
  theme_classic()
print (s)

ggsave(
  filename = paste0(outdir,"/scatt_len-vs-PWdiff_VMRs_chr.pdf"),
  plot = s,
  width = 8,
  height = 3
  # ,units = "in"
)


# Plot scatter of Lenght vs # of DMRs  -----------------------------------------
d <- ggplot(dt, aes(x = VMR_length, y = n_dmrs)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "VMR length (bp)",
    y = "# of DMRs"
  ) +
  theme_classic()
print (d)

ggsave(
  filename = paste0(outdir,"/scatt_len-vs-Ndmr_VMRs.pdf"),
  plot = d,
  width = 8,
  height = 3
  # ,units = "in"
)

# {By chromosome} Plot scatter of Lenght vs # of DMRs  -------------------------
d <- ggplot(dt, aes(x = VMR_length, y = n_dmrs, color=chromosome_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "VMR length (bp)",
    y = "# of DMRs"
  ) +
  theme_classic()
print (d)

ggsave(
  filename = paste0(outdir,"/scatt_len-vs-Ndmr_VMRs_chr.pdf"),
  plot = d,
  width = 8,
  height = 3
  # ,units = "in"
)



# Plot scatter of max_pairwise_diff vs # of DMRs -------------------------------
n <- ggplot(dt, aes(x = n_dmrs,y = max_pairwise_diff )) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "# of DMRs",
    y = "Maximum pairwise methylation difference"
  ) +
  theme_classic()
print (n)

ggsave(
  filename = paste0(outdir,"/scatt_PWdiff-vs-Ndmr_VMRs.pdf"),
  plot = n,
  width = 8,
  height = 3
  # ,units = "in"
)

# {By chromosome} plot scatter of max_pairwise_diff vs # of DMRs ---------------
n <- ggplot(dt, aes(x = n_dmrs,y = max_pairwise_diff, color = chromosome_type)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "# of DMRs",
    y = "Maximum pairwise methylation difference"
  ) +
  theme_classic()
print (n)

ggsave(
  filename = paste0(outdir,"/scatt_PWdiff-vs-Ndmr_VMRs_chr.pdf"),
  plot = n,
  width = 8,
  height = 3
  # ,units = "in"
)


# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Prepare cumulative chromosome positions using true chromosome lengths
# from the genome FASTA index (.fai).
#
# Chromosomes are explicitly ordered Parae_01 -> Parae_23.

fai <- fread(
  fai_file,
  header = FALSE,
  select = c(1, 2),
  col.names = c(
    "chr",
    "chr_length"
  )
)

# Keep only the 23 main chromosomes
fai <- fai[
  grepl(
    "^Parae_[0-9]+$",
    chr
  )
]

# Extract chromosome number and retain Parae_01 -> Parae_23
fai[
  ,
  chr_number := as.integer(
    sub(
      "^Parae_",
      "",
      chr
    )
  )
]

fai <- fai[
  chr_number %between% c(1, 23)
]

setorder(
  fai,
  chr_number
)

# Check that all 23 chromosomes are present
expected_chrs <- sprintf(
  "Parae_%02d",
  1:23
)

missing_chrs <- setdiff(
  expected_chrs,
  fai$chr
)

if (
  length(missing_chrs) > 0
) {
  
  stop(
    "Missing chromosomes from .fai file: ",
    paste(
      missing_chrs,
      collapse = ", "
    )
  )
}

# Explicitly preserve chromosome order
fai[
  ,
  chr := factor(
    chr,
    levels = expected_chrs
  )
]

setorder(
  fai,
  chr
)

# -------------------------------------------------------------------------------
# Keep VMRs located on the 23 main chromosomes
# -------------------------------------------------------------------------------

dt <- dt[
  chr %in% expected_chrs
]

dt[
  ,
  chr := factor(
    chr,
    levels = expected_chrs
  )
]

setorder(
  dt,
  chr,
  start,
  end
)

# -------------------------------------------------------------------------------
# Calculate genome-wide cumulative chromosome coordinates
# -------------------------------------------------------------------------------

chr_info <- copy(fai)

chr_info[
  ,
  offset := cumsum(
    shift(
      chr_length,
      fill = 0
    )
  )
]

chr_info[
  ,
  center := offset + chr_length / 2
]

# Chromosome boundaries for vertical lines
chr_boundaries <- chr_info[
  -1,
  offset
]

# -------------------------------------------------------------------------------
# Add chromosome offsets to each VMR
# -------------------------------------------------------------------------------

dt <- merge(
  dt,
  chr_info[
    ,
    .(
      chr,
      offset
    )
  ],
  by = "chr",
  all.x = TRUE,
  sort = FALSE
)

dt[
  ,
  genome_pos := start + offset
]

# Restore genomic order after merge
setorder(
  dt,
  chr,
  start,
  end
)

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Manhattan plot: genomic position vs maximum pairwise methylation difference

m <- ggplot(
  dt,
  aes(
    x = genome_pos,
    y = max_pairwise_diff
  )
) +
  
  geom_vline(
    xintercept = chr_boundaries,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_point(
    aes(
      color = chr
    ),
    size = 1
  ) +
  
  scale_x_continuous(
    breaks = chr_info$center,
    labels = as.character(
      chr_info$chr
    ),
    expand = expansion(
      mult = c(
        0.002,
        0.002
      )
    )
  ) +
  
  labs(
    x = "Chromosome",
    y = "Maximum pairwise methylation difference"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(m)

ggsave(
  filename = paste0(
    outdir,
    "/manhattan_PWdiff-vs-pos_VMRs.pdf"
  ),
  plot = m,
  width = 8,
  height = 3
)
