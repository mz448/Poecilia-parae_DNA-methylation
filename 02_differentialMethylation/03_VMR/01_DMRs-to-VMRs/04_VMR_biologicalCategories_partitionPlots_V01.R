#!/usr/bin/env Rscript
# DATE: 20260814
# AUTHOR: MZF
# SCRIPT: 03_VMR_biologicalCategories_partitionPlots_V01.R
# VERSION: 01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Visualize the distribution of biologically defined VMR contrast-sharing
#   categories across genomic partitions.
#
#   The script generates:
#
#     1) Pie charts showing the percentage of VMRs assigned to each biological
#        contrast-sharing category for:
#          - Whole genome
#          - Autosomes
#          - Sex chromosome (Parae_12)
#
#     2) A bar plot comparing the percentage representation of each VMR
#        category between autosomes and the sex chromosome.
#
#        Percentage-point difference is calculated as:
#
#          Autosomes (%) - Sex chromosome (%)
#
#        Therefore:
#          - Positive values indicate enrichment on autosomes.
#          - Negative values indicate enrichment on the sex chromosome.
#          - Zero indicates equal proportional representation.
#
# BIOLOGICAL CATEGORIES:
#   - Sex-associated
#   - Sex + morph
#   - Cross-phenotype
#   - Morph-specific
#   - Pairwise-specific
#   - Other contrast combinations
#
# INPUT:
#   vmr_real.contrastSharing.biologicalCategories.summary.tsv
#
# EXPECTED INPUT COLUMNS:
#   partition
#   category_order
#   biological_category
#   N
#   percentage
#   total_vmrs
#
# OUTPUT:
#   vmr_real.biologicalCategories.partitionPie.pdf
#   vmr_real.biologicalCategories.AUTOminusChr12.pdf
#
# REQUIRED R PACKAGES:
#   data.table
#   ggplot2
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%


library(data.table)
library(ggplot2)


# To start fresh ---------------------------------------------------------------
graphics.off()
rm(list = ls())

# Input ------------------------------------------------------------------------


# input_file <- "../plots/complexupset/vmr_real.contrastSharing.biologicalCategories.summary.tsv"
input_file <- "../plots/complexupset/vmr_shuffled_all_morph.contrastSharing.biologicalCategories.summary.tsv"

dt <- fread(
  input_file
)

outdir <- "../plots/complexupset/"


# Category order


category_levels <- c(
  "Sex-associated",
  "Sex + morph",
  "Cross-phenotype",
  "Morph-specific",
  "Pairwise-specific",
  "Other contrast combinations"
)

dt[
  ,
  biological_category := factor(
    biological_category,
    levels = category_levels
  )
]



# Partition labels


dt[
  ,
  partition_label := factor(
    partition,
    levels = c(
      "Whole_genome",
      "Autosomes",
      "Sex_chromosome"
    ),
    labels = c(
      "Whole genome",
      "Autosomes",
      "Chromosome 12"
    )
  )
]



# Category colors --------------------------------------------------------------


category_colors <- c(
  "Sex-associated" = "#4C78A8",
  "Sex + morph" = "#F58518",
  "Cross-phenotype" = "#54A24B",
  "Morph-specific" = "#B279A2",
  "Pairwise-specific" = "#E45756",
  "Other contrast combinations" = "#B8B8B8"
)



# Pie chart --------------------------------------------------------------------


p_pie <- ggplot(
  dt,
  aes(
    x = "",
    y = percentage,
    fill = biological_category
  )
) +
  
  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.4
  ) +
  
  coord_polar(
    theta = "y"
  ) +
  
  facet_wrap(
    ~ partition_label,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = category_colors,
    drop = FALSE
  ) +
  
  labs(
    fill = "VMR category"
  ) +
  
  theme_void() +
  
  theme(
    strip.text = element_text(
      size = 12,
      face = "bold"
    ),
    legend.title = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 10
    )
  )


print(
  p_pie
)


# output
ggsave(
  "../plots/complexupset/vmr_real.biologicalCategories.partitionPie.pdf",
  p_pie,
  width = 11,
  height = 4.5
)

# Bar chart of Percetn differences ---------------------------------------------


# Retain autosomes and chromosome 12


dt_compare <- dt[
  partition %in% c(
    "Autosomes",
    "Sex_chromosome"
  )
]



# Convert to wide format


dt_wide <- dcast(
  dt_compare,
  biological_category ~ partition,
  value.var = "percentage"
)



# Calculate percentage-point difference
#
# Positive = greater percentage on autosomes
# Negative = greater percentage on chromosome 12


  dt_wide[
    ,
    percentage_difference :=
      Autosomes -
      Sex_chromosome
  ]



# Category order


  dt_wide[
    ,
    biological_category := factor(
      biological_category,
      levels = category_levels
    )
  ]



# Plot -------

# 
# p_difference <- ggplot(
#   dt_wide,
#   aes(
#     x = biological_category,
#     y = percentage_difference
#   )
# ) +
#   
#   geom_hline(
#     yintercept = 0,
#     linewidth = 0.6
#   ) +
#   
#   geom_col(
#     width = 0.7,
#     fill = "grey40"
#   ) +
#   
#   labs(
#     x = "VMR category",
#     y = "Percentage-point difference\nAutosomes - X chromosome"
#   ) +
#   
#   scale_y_continuous(
#     expand = expansion(
#       mult = c(
#         0.08,
#         0.08
#       )
#     )
#   ) +
#   
#   theme_classic() +
#   
#   theme(
#     axis.title = element_text(
#       size = 12
#     ),
#     axis.text.y = element_text(
#       size = 10
#     ),
#     axis.text.x = element_text(
#       size = 10,
#       angle = 45,
#       hjust = 1
#     )
#   )
# 
# 
# print(
#   p_difference
# )
# 
# 
# # Optional output
# ggsave(
#   "../plots/complexupset/vmr_real.biologicalCategories.AUTOminusChr12.pdf",
#   p_difference,
#   width = 8,
#   height = 5
# )


# Plot Bars --------------------------------------------------------------------
  
p_difference <- ggplot(
  dt_wide,
  aes(
    x = biological_category,
    y = percentage_difference,
    fill = biological_category
  )
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.6
  ) +
  
  geom_col(
    width = 0.7
  ) +
  
  scale_fill_manual(
    values = category_colors,
    drop = FALSE
  ) +
  
  scale_y_continuous(
    limits = c(
      -20,
      20
    ),
    breaks = seq(
      -20,
      20,
      by = 2.5
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  labs(
    x = "VMR category",
    y = "Percentage-point difference\nAutosomes - X chromosome"
  ) +
  
  theme_classic() +
  
  theme(
    axis.title = element_text(
      size = 12
    ),
    axis.text.y = element_text(
      size = 10
    ),
    axis.text.x = element_text(
      size = 10,
      angle = 45,
      hjust = 1
    ),
    legend.position = "none"
  )


print(
  p_difference
)


# Output
ggsave(
  "../plots/complexupset/vmr_real.biologicalCategories.AUTOminusChr12.pdf",
  p_difference,
  width = 5,
  height = 5
)
