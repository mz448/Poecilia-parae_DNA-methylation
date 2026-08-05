#!/usr/bin/env Rscript

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# AUTHOR: MZF
# DATE: 20260720
# SCRIPT: 04_plotCoverage_V01.R
# VERSION: 01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# DESCRIPTION:
#   - The script reads coverage_CpG.tsv and calculates the global weighted mean 
#     CpG coverage across all samples.
#   - Each sample’s mean CpG depth as a box-shaped point, with whiskers extending 
#     from Q1 to Q3.
#   - Dashed horizontal line at the global mean coverage 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

library(ggplot2)

# ==============================================================================
# Input
# ==============================================================================

coverage <- read.delim("coverage_CpG_V01.tsv", header = TRUE)

coverage$sample <- factor(coverage$sample, levels = coverage$sample)

# Global mean coverage across the library
global_mean <- weighted.mean(
  coverage$mean_CpG_depth,
  coverage$reported_CpG_positions
)

coverage$sample <- sub("^Sym\\.forSym\\.", "", coverage$sample)
coverage$sample <- sub("\\.CpG_report\\.txt$", "", coverage$sample)

coverage$sample <- factor(
  coverage$sample,
  levels = coverage$sample
)

# ==============================================================================
# Plot
# ==============================================================================

p <- ggplot(
  coverage,
  aes(
    # x = reorder(sample, mean_CpG_depth),
    x = sample,
    y = mean_CpG_depth
  )
) +
  # geom_errorbar(
  #   aes(
  #     ymin = Q1_CpG_depth,
  #     ymax = Q3_CpG_depth
  #   ),
  #   width = 0.3,
  #   linewidth = 0.6
  # ) +
  # geom_point(
  #   shape = 22,
  #   size = 3
  # ) +
  geom_crossbar(
    aes(
      ymin = Q1_CpG_depth,
      ymax = Q3_CpG_depth
    ),
    width = 0.6,
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = global_mean,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  labs(
    x = "Sample",
    y = "CpG coverage depth"
  ) +
  coord_flip(
    ylim = c(0, 40)
  ) +
  scale_y_continuous(breaks = seq(0, 40, by = 5))+
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )

  )

print (p)


ggsave(
  "CpG_coverage_by_sample.pdf",
  plot = p,
  width = 5.7,
  height =8.5
)
