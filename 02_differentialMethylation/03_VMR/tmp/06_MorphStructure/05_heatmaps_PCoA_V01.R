#!/usr/bin/env Rscript

# DATE:   20260803
# AUTHOR: MZF
# SCRIPT: 05_heatmaps_PCoA_V01.R
# VERSION:V01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:   Plot Condition VMR similarity using 
#         1) PCoAs
#         2) Euclidean sample-distance Heatmaps
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

library(data.table)
library(vegan)
library(ggplot2)
library(pheatmap)
library(viridis)

# To start fresh ---------------------------------------------------------------
# Clear all plots and graphics devices
graphics.off()
# Remove all objects from the workspace environment
rm(list = ls())

# Select one condition ---------------------------------------------------------
condition <- "real_morph"
# condition <- "shuffled_morph"
# condition <- "real_sex"
# condition <- "shuffled_sex"
# condition <- "shuffled_all"

# To plot the VMRs used in the rest of the analysis: 
# condition <- "real_Threshold-0" # uncomment to test the whole real VMR dataset
# condition <- "real_Threshold-0.4" # uncomment to test the filtered real VMR dataset

# Prepare data -----------------------------------------------------------------
out_dir <- paste0("../plots/dataset_Structure/",condition)
dir.create(out_dir, recursive = TRUE)

# Read VMR table
dt <- fread(paste0("./dataset_Structure/",condition,"_vmr_methylation/vmr.",condition,".morph_and_sample_observed.tsv"))

# dt <- fread("./vmr_methylation/vmr.morph_and_sample_observed.tsv") # Uncomment to test the whole real VMR dataset
# dt <- fread("../plots/0.4/vmr.morph_and_sample_observed.filtered.diffGE_0.4.tsv") # Uncomment to test the filtered real VMR dataset

# Individual-sample methylation columns
meth_cols <- grep("^meth_ppar", names(dt), value = TRUE)

# Optionally retain only VMRs measured in all 12 samples
dt <- dt[n_meth_nonNA_samples == 12]

# Convert to numeric matrix:
# current orientation = VMR rows × sample columns
meth_vmr_by_sample <- as.matrix(dt[, ..meth_cols])
storage.mode(meth_vmr_by_sample) <- "numeric"

# Transpose:
# required orientation = sample rows × VMR columns
meth_sample_by_vmr <- t(meth_vmr_by_sample)

# Clean sample names
rownames(meth_sample_by_vmr) <- sub(
  "^meth_",
  "",
  rownames(meth_sample_by_vmr)
)

# Populate Metadata ------------------------------------------------------------

# Sample metadata
samples <- rownames(meth_sample_by_vmr)

metadata <- data.frame(
  sample = samples,
  stringsAsFactors = FALSE
)

metadata$sex <- ifelse(
  grepl("^pparfmem", metadata$sample),
  "female",
  "male"
)

metadata$phenotype <- ifelse(
  grepl("^pparfmem", metadata$sample), "female",
  ifelse(
    grepl("^pparimem", metadata$sample), "immaculata",
    ifelse(
      grepl("^pparpmem", metadata$sample), "parae",
      ifelse(
        grepl("^pparymem", metadata$sample), "yellow",
        NA_character_
      )
    )
  )
)

metadata$male_morph <- ifelse(
  metadata$sex == "male",
  metadata$phenotype,
  NA_character_
)

metadata$sex <- factor(
  metadata$sex,
  levels = c("female", "male")
)

metadata$phenotype <- factor(
  metadata$phenotype,
  levels = c("female", "immaculata", "parae", "yellow")
)

metadata$male_morph <- factor(
  metadata$male_morph,
  levels = c("immaculata", "parae", "yellow")
)

stopifnot(!anyNA(metadata$sex))
stopifnot(!anyNA(metadata$phenotype))

# Save metadata
write.csv(metadata, file = file.path(out_dir,paste0("PERMANOVA.metadata.",condition,".csv")), row.names = TRUE)

# Euclidean distances between samples
d <- dist(meth_sample_by_vmr, method = "euclidean")

# Subset metadata and methylaiton matrix.
male_index <- metadata$sex == "male"
meth_males <- meth_sample_by_vmr[male_index,,drop = FALSE]
metadata_males <- droplevels(metadata[male_index, ])
# Calculate euclidean distances between male-samples
d_males <- dist(meth_males,method = "euclidean")



# Set Up custom colors
# A Distance color scale
n_colors <- 100
# color_palette <- mako(n_colors)
color_palette <- colorRampPalette(c("black", "white"))(n_colors)


# B Factor color scale
sex_colors <- c(
  female = "#399957",
  male   = "#99397B"
)

phenotype_colors <- c(
  female      = "#399957",
  immaculata  = "#888888",
  parae       = "#804e9f",
  yellow      = "#fdd800"
)

male_morph_colors <- c(
  immaculata  = "#888888",
  parae       = "#804e9f",
  yellow      = "#fdd800"
)

annotation_colors_all <- list(
  sex = sex_colors,
  phenotype = phenotype_colors
)

annotation_colors_males <- list(
  male_morph = male_morph_colors
)



# 1) PCoA and distance heatmap: all samples ------------------------------------

# Confirm that sample order matches between the distance matrix and metadata
stopifnot(identical(labels(d), metadata$sample))

# Add sample names as metadata row names for heatmap annotations
rownames(metadata) <- metadata$sample

# Classical PCoA from the same Euclidean distance matrix used in PERMANOVA
pcoa_all <- cmdscale(
  d,
  k = 2,
  eig = TRUE
)

# Percentage represented by each PCoA axis
positive_eigenvalues <- pcoa_all$eig[pcoa_all$eig > 0]

pcoa1_percent <- 100 * pcoa_all$eig[1] / sum(positive_eigenvalues)
pcoa2_percent <- 100 * pcoa_all$eig[2] / sum(positive_eigenvalues)

# Combine coordinates with sample metadata
pcoa_all_df <- data.frame(
  sample = metadata$sample,
  PCoA1 = pcoa_all$points[metadata$sample, 1],
  PCoA2 = pcoa_all$points[metadata$sample, 2],
  sex = metadata$sex,
  phenotype = metadata$phenotype
)

# Save coordinates
write.csv(
  pcoa_all_df,
  file = file.path(
    out_dir,
    paste0("PCoA.coordinates.all_samples.", condition, ".csv")
  ),
  row.names = FALSE
)

# PCoA plot
pcoa_all_plot <- ggplot(
  pcoa_all_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = phenotype,
    shape = sex
  )
) +
  geom_point(size = 4) +
  geom_polygon(
    aes(group = phenotype, fill = phenotype, color=phenotype),
    alpha = 0.5,
    color = NA,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = sample),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = phenotype_colors,
    drop = FALSE
  ) +
  scale_color_manual(
    values = phenotype_colors,
    drop = FALSE
  ) +
  labs(
    title = paste("PCoA of VMR methylation profiles:", condition),
    x = sprintf("PCoA1 (%.1f%%)", pcoa1_percent),
    y = sprintf("PCoA2 (%.1f%%)", pcoa2_percent),
    color = "Phenotype",
    shape = "Sex"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    out_dir,
    paste0("PCoA.all_samples.", condition, ".pdf")
  ),
  plot = pcoa_all_plot,
  width = 7,
  height = 5
)

# Euclidean sample-distance heatmap
distance_matrix <- as.matrix(d)

annotation_all <- metadata[, c("sex", "phenotype"), drop = FALSE]

sample_clustering <- hclust(
  d,
  method = "complete"
)

# Plot heatmap
pheatmap(
  distance_matrix,
  cluster_rows = sample_clustering,
  cluster_cols = sample_clustering,
  annotation_row = annotation_all,
  annotation_col = annotation_all,
  color = color_palette,
  annotation_colors = annotation_colors_all,
  main = paste("Euclidean sample distances:", condition),
  border_color = NA,
  filename = file.path(
    out_dir,
    paste0("Heatmap.sample_distances.", condition, ".pdf")
  ),
  width = 8,
  height = 7
)

# 2) PCoA and distance heatmap: males only -------------------------------------

stopifnot(identical(labels(d_males), metadata_males$sample))

rownames(metadata_males) <- metadata_males$sample

pcoa_males <- cmdscale(
  d_males,
  k = 2,
  eig = TRUE
)

positive_eigenvalues_males <- pcoa_males$eig[pcoa_males$eig > 0]

male_pcoa1_percent <- (
  100 * pcoa_males$eig[1] / sum(positive_eigenvalues_males)
)

male_pcoa2_percent <- (
  100 * pcoa_males$eig[2] / sum(positive_eigenvalues_males)
)

pcoa_males_df <- data.frame(
  sample = metadata_males$sample,
  PCoA1 = pcoa_males$points[metadata_males$sample, 1],
  PCoA2 = pcoa_males$points[metadata_males$sample, 2],
  male_morph = metadata_males$male_morph
)

write.csv(
  pcoa_males_df,
  file = file.path(
    out_dir,
    paste0("PCoA.coordinates.males.", condition, ".csv")
  ),
  row.names = FALSE
)

pcoa_males_plot <- ggplot(
  pcoa_males_df,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = male_morph
  )
) +
  geom_point(size = 4) +
  geom_polygon(
    aes(group = male_morph, fill = male_morph, color = male_morph),
    alpha = 0.12,
    color = NA,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = sample),
    vjust = -0.8,
    size = 3,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = male_morph_colors, drop = FALSE) +
  scale_color_manual(values = male_morph_colors, drop = FALSE) +
  labs(
    title = paste("Male-only PCoA:", condition),
    x = sprintf("PCoA1 (%.1f%%)", male_pcoa1_percent),
    y = sprintf("PCoA2 (%.1f%%)", male_pcoa2_percent),
    color = "Male morph"
  ) +
  theme_classic()

ggsave(
  filename = file.path(
    out_dir,
    paste0("PCoA.males.", condition, ".pdf")
  ),
  plot = pcoa_males_plot,
  width = 7,
  height = 5
)

male_distance_matrix <- as.matrix(d_males)

annotation_males <- metadata_males[
  ,
  "male_morph",
  drop = FALSE
]

male_clustering <- hclust(
  d_males,
  method = "complete"
)

pheatmap(
  male_distance_matrix,
  cluster_rows = male_clustering,
  cluster_cols = male_clustering,
  annotation_row = annotation_males,
  annotation_col = annotation_males,
  color = color_palette,
  annotation_colors = annotation_colors_males,
  main = paste("Male-only Euclidean distances:", condition),
  border_color = NA,
  filename = file.path(
    out_dir,
    paste0("Heatmap.male_sample_distances.", condition, ".pdf")
  ),
  width = 7,
  height = 6
)