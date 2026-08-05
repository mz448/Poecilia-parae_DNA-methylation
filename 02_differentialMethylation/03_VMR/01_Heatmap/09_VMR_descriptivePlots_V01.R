#!/usr/bin/env Rscript
# DATE:       2026-05-11
# AUTHOR:     MZF
# SCRIPT:     10_VMR_descriptivePlots_V02.R
# VERSION:    02
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
infile <- "vmr.morph_and_sample_observed.filtered.diffGE_0.4.tsv"   
outdir      <- "VMR_tracks"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Load Dataframe:
dt <- fread(infile)
dt[, VMR_length := end - start]


# Plot histogram of length. 
h <- ggplot(dt, aes(x = VMR_length)) +
  geom_histogram(bins = 90) +
  labs(
    x = "VMR length (bp)",
    y = "Count"
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

# Plot scatter of Lenght vs max_pairwise_diff
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

# Plot scatter of Lenght vs # of DMRs
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



# Plot scatter of max_pairwise_diff vs # of DMRs
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




# Prepare cumulative chromosome positions
chr_info <- dt[, .(chr_len = max(end)), by = chr]
chr_info[, offset := c(0, cumsum(chr_len)[- .N])]
chr_info[, center := offset + chr_len / 2]

dt <- merge(dt, chr_info[, .(chr, offset)], by = "chr")
dt[, genome_pos := start + offset]

# Manhattan of position and PWdiff
m <- ggplot(dt, aes(x = genome_pos, y = max_pairwise_diff)) +
  geom_vline(
    xintercept = chr_info$offset[-1],
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_point(aes(color = chr), size = 1) +
  scale_x_continuous(
    breaks = chr_info$center,
    labels = chr_info$chr
  ) +
  labs(
    x = "Chromosome",
    y = "Maximum pairwise methylation difference"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print (m)

ggsave(
  filename = paste0(outdir,"/manhattan_PWdiff-vs-pos_VMRs.pdf"),
  plot = m,
  width = 8,
  height = 3
  # ,units = "in"
)
