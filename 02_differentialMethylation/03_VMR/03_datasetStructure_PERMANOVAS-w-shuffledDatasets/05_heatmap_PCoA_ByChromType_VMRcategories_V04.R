#!/usr/bin/env Rscript

# DATE:   20260814
# AUTHOR: MZF
# SCRIPT: 05_heatmap_PCoA_ByChromType_VMRcategories_V04.R
# VERSION: V04
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:   Plot condition VMR similarity independently for:
#         1) Autosomal VMRs
#         2) X-chromosome VMRs (Parae_12)
#
#         For each genomic partition:
#         A) PCoAs
#         B) Euclidean sample-distance heatmaps
#
#         Analyses are produced for:
#         - all 12 samples
#         - males only
#
# ADDITIONAL FILTER:
#         - Only VMRs that are part of the selected categories defined in the 
#           UPSET analysis
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

library(data.table)
library(vegan)
library(ggplot2)
library(pheatmap)
library(viridis)

# To start fresh ---------------------------------------------------------------
graphics.off()
rm(list = ls())


# Select VMR biological categories ---------------------------------------------
#
# These categories are defined by the UpSet contrast-sharing analysis.
# Keep one or several categories active. All selected categories are pooled
# before the sample-distance, PCoA, and heatmap analyses.

VMR_CATEGORIES <- c(
  "Sex-associated",
  "Sex + morph",
  "Cross-phenotype",
  "Morph-specific",
  "Pairwise-specific",
  "Other contrast combinations"
)

category_annotation_file <- "../plots/complexupset/vmr_real.contrastSharing.biologicalCategories.perVMR.tsv"


# Select clustering method -----------------------------------------------------
method_h_clust <- "complete"
# method_h_clust <- "average"
# method_h_clust <- "ward.D2"
# method_h_clust <- "centroid"


# Select one condition ---------------------------------------------------------
#
# Keep ONE of the input-selection blocks active
#
## A) CONDITION SPECIFIC
# condition <- "real_morph"
# condition <- "shuffled_morph"
# condition <- "real_sex"
# condition <- "shuffled_sex"
# condition <- "shuffled_all"
# th <- "0"
# dt <- fread(
#   paste0(
#     "../plots/tmp/real_contrast_VMRs/",
#     th, ".", condition,
#     "/vmr.morph_and_sample_observed.filtered.diffGE_",
#     th, ".tsv"
#   )
# )
#
## B) ALL CONDITIONS (together)
# condition0 <- "allContrasts"
# th <- "0"
# condition <- paste0(condition0, "_", th)
# dt <- fread(
#   paste0(
#     "../plots/",
#     condition,
#     "/vmr.morph_and_sample_observed.filtered.diffGE_",
#     th, ".tsv"
#   )
# )
#
# C) SUBSET OF CONDITIONS (real conditions merged)
th <- "0"
condition0 <- "real"
condition <- paste0(condition0, "_", th)
dt <- fread(
  paste0(
    "../plots/",
    th,
    "/vmr.morph_and_sample_observed.filtered.diffGE_",
    th, ".tsv"
  )
)
#
# # D) FILTERED SUBSET OF CONDITIONS (real conditions merged)
# th <- "0.4"
# condition <- "real"
# dt <- fread(
#   paste0(
#     "../plots/",
#     th,
#     "/vmr.morph_and_sample_observed.filtered.diffGE_",
#     th, ".tsv"
#   )
# )
#
##  CUSTOM EXAMPLE
# condition <- "real"
# dt <- fread(
#   "/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/plots/0/indicatorVMR.filter/indicatorVMR.filtered.vmr.morph_and_sample_observed.RG0.5_M0.1_P0.05_MAX500.tsv"
# )


# Check that an input dataset was selected -------------------------------------
if (!exists("dt") || !exists("condition")) {
  stop(
    "No input dataset selected. Activate one condition/input block above ",
    "so that both 'dt' and 'condition' are defined."
  )
}


# Join and filter VMR biological categories ------------------------------------

expected_categories <- c(
  "Sex-associated",
  "Sex + morph",
  "Cross-phenotype",
  "Morph-specific",
  "Pairwise-specific",
  "Other contrast combinations"
)

bad_categories <- setdiff(VMR_CATEGORIES, expected_categories)
if (length(bad_categories)) stop("Unknown VMR categories: ", paste(bad_categories, collapse = ", "))

if (!"cons_id" %in% names(dt)) stop("Input methylation table must contain 'cons_id'.")

# Add biological_category only when it is not already present in the input.
if (!"biological_category" %in% names(dt)) {
  anno <- fread(category_annotation_file)
  
  required_anno <- c("cons_id", "biological_category")
  missing_anno <- setdiff(required_anno, names(anno))
  if (length(missing_anno)) stop("Missing category annotation columns: ", paste(missing_anno, collapse = ", "))
  
  if (anyDuplicated(anno$cons_id)) stop("cons_id must be unique in the category annotation table.")
  if (anyDuplicated(dt$cons_id)) stop("cons_id must be unique in the methylation input table.")
  
  missing_ids <- setdiff(dt$cons_id, anno$cons_id)
  if (length(missing_ids)) {
    stop(
      length(missing_ids),
      " methylation VMRs are absent from the category annotation table."
    )
  }
  
  dt[anno, on = "cons_id", biological_category := i.biological_category]
}

if (anyNA(dt$biological_category)) stop("Some VMRs have no biological_category after annotation.")

n_vmr_before_category_filter <- nrow(dt)

category_counts_before <- dt[, .N, by = biological_category]
setorder(category_counts_before, biological_category)

dt <- dt[biological_category %in% VMR_CATEGORIES]

if (nrow(dt) == 0) stop("No VMRs remain after filtering by VMR_CATEGORIES.")

category_counts_selected <- dt[, .N, by = biological_category]
setorder(category_counts_selected, biological_category)

category_tag <- if (
  setequal(VMR_CATEGORIES, expected_categories)
) {
  "all_categories"
} else {
  paste(gsub("[^A-Za-z0-9]+", "_", VMR_CATEGORIES), collapse = "-")
}

category_title <- paste(VMR_CATEGORIES, collapse = ", ")

message("\nVMR category filter")
message("VMRs before category filtering: ", n_vmr_before_category_filter)
message("VMRs after category filtering:  ", nrow(dt))
message("Selected categories: ", category_title)


# Prepare shared data ----------------------------------------------------------

out_dir_base <- file.path(
  "../plots/dataset_Structure/",
  paste0(condition, "_", method_h_clust),
  category_tag
)

dir.create(
  out_dir_base,
  recursive = TRUE,
  showWarnings = FALSE
)

fwrite(
  category_counts_before,
  file = file.path(out_dir_base, paste0("VMR_counts.by_category.before_filter.", condition, ".tsv")),
  sep = "\t"
)

fwrite(
  category_counts_selected,
  file = file.path(out_dir_base, paste0("VMR_counts.by_category.selected.", condition, ".tsv")),
  sep = "\t"
)


# Individual-sample methylation columns
meth_cols <- grep(
  "^meth_ppar",
  names(dt),
  value = TRUE
)

if (length(meth_cols) != 12) {
  stop(
    "Expected 12 individual-sample methylation columns beginning with ",
    "'meth_ppar', but found ",
    length(meth_cols),
    ":\n",
    paste(meth_cols, collapse = "\n")
  )
}


# Populate metadata ------------------------------------------------------------

samples <- sub(
  "^meth_",
  "",
  meth_cols
)

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

rownames(metadata) <- metadata$sample


# Save metadata once
write.csv(
  metadata,
  file = file.path(
    out_dir_base,
    paste0("PERMANOVA.metadata.", condition, ".csv")
  ),
  row.names = TRUE
)


# Set up custom colors ---------------------------------------------------------

# A) Distance color scale
n_colors <- 100
color_palette <- colorRampPalette(
  c("black", "white")
)(n_colors)

# B) Factor color scales
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
  immaculata = "#888888",
  parae      = "#804e9f",
  yellow     = "#fdd800"
)

annotation_colors_all <- list(
  sex = sex_colors,
  phenotype = phenotype_colors
)

annotation_colors_males <- list(
  male_morph = male_morph_colors
)


# Define genomic partitions ----------------------------------------------------
#
# Autosomes are selected using chromosome_type.
# The X chromosome is selected directly by chromosome name so the code does not
# depend on the exact label used for chromosome_type on Parae_12.

partition_data <- list(
  Autosome = dt[
    chromosome_type == "Autosome"
  ],
  X_Parae_12 = dt[
    chr == "Parae_12"
  ]
)

partition_titles <- c(
  Autosome = "Autosomes",
  X_Parae_12 = "X chromosome (Parae_12)"
)


# Function: run PCoA + heatmaps for one genomic partition ----------------------

run_partition_analysis <- function(
    dt_partition,
    partition_label,
    partition_title
) {
  
  message(
    "\nRunning partition: ",
    partition_title
  )
  
  n_vmr_before_complete_filter <- nrow(
    dt_partition
  )
  
  if (n_vmr_before_complete_filter == 0) {
    warning(
      "No VMRs found for partition: ",
      partition_title,
      ". Skipping."
    )
    return(
      data.frame(
        partition = partition_label,
        vmr_categories = category_title,
        n_vmr_before_complete_filter = 0,
        n_vmr_complete_12_samples = 0
      )
    )
  }
  
  
  # Retain only VMRs measured in all 12 samples
  dt_partition <- dt_partition[
    n_meth_nonNA_samples == 12
  ]
  
  n_vmr_complete <- nrow(
    dt_partition
  )
  
  if (n_vmr_complete == 0) {
    warning(
      "No VMRs measured in all 12 samples for partition: ",
      partition_title,
      ". Skipping."
    )
    return(
      data.frame(
        partition = partition_label,
        vmr_categories = category_title,
        n_vmr_before_complete_filter = n_vmr_before_complete_filter,
        n_vmr_complete_12_samples = 0
      )
    )
  }
  
  
  # Partition-specific output directory
  out_dir <- file.path(
    out_dir_base,
    partition_label
  )
  
  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # Convert to numeric matrix:
  # current orientation = VMR rows × sample columns
  meth_vmr_by_sample <- as.matrix(
    dt_partition[
      ,
      ..meth_cols
    ]
  )
  
  storage.mode(
    meth_vmr_by_sample
  ) <- "numeric"
  
  if (anyNA(meth_vmr_by_sample)) {
    stop(
      "NA methylation values remain after filtering ",
      partition_title,
      " to n_meth_nonNA_samples == 12."
    )
  }
  
  
  # Transpose:
  # required orientation = sample rows × VMR columns
  meth_sample_by_vmr <- t(
    meth_vmr_by_sample
  )
  
  rownames(
    meth_sample_by_vmr
  ) <- sub(
    "^meth_",
    "",
    rownames(meth_sample_by_vmr)
  )
  
  
  # Confirm sample identities and order
  stopifnot(
    identical(
      rownames(meth_sample_by_vmr),
      metadata$sample
    )
  )
  
  
  # Euclidean distances between all samples
  d <- dist(
    meth_sample_by_vmr,
    method = "euclidean"
  )
  
  
  # Male-only subset
  male_index <- metadata$sex == "male"
  
  meth_males <- meth_sample_by_vmr[
    male_index,
    ,
    drop = FALSE
  ]
  
  metadata_males <- droplevels(
    metadata[
      male_index,
      ,
      drop = FALSE
    ]
  )
  
  d_males <- dist(
    meth_males,
    method = "euclidean"
  )
  
  
  # ---------------------------------------------------------------------------
  # 1) PCoA and distance heatmap: all samples
  # ---------------------------------------------------------------------------
  
  stopifnot(
    identical(
      labels(d),
      metadata$sample
    )
  )
  
  
  # Classical PCoA from the same partition-specific Euclidean distance matrix
  pcoa_all <- cmdscale(
    d,
    k = 2,
    eig = TRUE
  )
  
  positive_eigenvalues <- pcoa_all$eig[
    pcoa_all$eig > 0
  ]
  
  pcoa1_percent <- (
    100 *
      pcoa_all$eig[1] /
      sum(positive_eigenvalues)
  )
  
  pcoa2_percent <- (
    100 *
      pcoa_all$eig[2] /
      sum(positive_eigenvalues)
  )
  
  
  pcoa_all_df <- data.frame(
    sample = metadata$sample,
    PCoA1 = pcoa_all$points[
      metadata$sample,
      1
    ],
    PCoA2 = pcoa_all$points[
      metadata$sample,
      2
    ],
    sex = metadata$sex,
    phenotype = metadata$phenotype
  )
  
  
  write.csv(
    pcoa_all_df,
    file = file.path(
      out_dir,
      paste0(
        "PCoA.coordinates.all_samples.",
        condition,
        ".",
        partition_label,
        ".csv"
      )
    ),
    row.names = FALSE
  )
  
  
  pcoa_all_plot <- ggplot(
    pcoa_all_df,
    aes(
      x = PCoA1,
      y = PCoA2,
      color = phenotype,
      shape = sex
    )
  ) +
    geom_point(
      size = 4
    ) +
    geom_polygon(
      aes(
        group = phenotype,
        fill = phenotype,
        color = phenotype
      ),
      alpha = 0.5,
      color = NA,
      show.legend = FALSE
    ) +
    geom_text(
      aes(
        label = sample
      ),
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
      title = paste0(
        "PCoA of VMR methylation profiles: ",
        condition,
        " - ",
        partition_title,
        "\nVMR categories: ",
        category_title
      ),
      x = sprintf(
        "PCoA1 (%.1f%%)",
        pcoa1_percent
      ),
      y = sprintf(
        "PCoA2 (%.1f%%)",
        pcoa2_percent
      ),
      color = "Phenotype",
      shape = "Sex"
    ) +
    theme_classic()
  
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        "PCoA.all_samples.",
        condition,
        ".",
        partition_label,
        ".pdf"
      )
    ),
    plot = pcoa_all_plot,
    width = 7,
    height = 5
  )
  
  
  # Euclidean sample-distance heatmap
  distance_matrix <- as.matrix(
    d
  )
  
  annotation_all <- metadata[
    ,
    c("sex", "phenotype"),
    drop = FALSE
  ]
  
  sample_clustering <- hclust(
    d,
    method = method_h_clust
  )
  
  
  pheatmap(
    distance_matrix,
    cluster_rows = sample_clustering,
    cluster_cols = sample_clustering,
    annotation_row = annotation_all,
    annotation_col = annotation_all,
    color = color_palette,
    annotation_colors = annotation_colors_all,
    main = paste0(
      "Euclidean sample distances: ",
      condition,
      " - ",
      partition_title,
      "\nVMR categories: ",
      category_title
    ),
    border_color = NA,
    filename = file.path(
      out_dir,
      paste0(
        "Heatmap.sample_distances.",
        condition,
        ".",
        partition_label,
        ".pdf"
      )
    ),
    width = 8,
    height = 7
  )
  
  
  # ---------------------------------------------------------------------------
  # 2) PCoA and distance heatmap: males only
  # ---------------------------------------------------------------------------
  
  stopifnot(
    identical(
      labels(d_males),
      metadata_males$sample
    )
  )
  
  rownames(
    metadata_males
  ) <- metadata_males$sample
  
  
  pcoa_males <- cmdscale(
    d_males,
    k = 2,
    eig = TRUE
  )
  
  positive_eigenvalues_males <- pcoa_males$eig[
    pcoa_males$eig > 0
  ]
  
  male_pcoa1_percent <- (
    100 *
      pcoa_males$eig[1] /
      sum(positive_eigenvalues_males)
  )
  
  male_pcoa2_percent <- (
    100 *
      pcoa_males$eig[2] /
      sum(positive_eigenvalues_males)
  )
  
  
  pcoa_males_df <- data.frame(
    sample = metadata_males$sample,
    PCoA1 = pcoa_males$points[
      metadata_males$sample,
      1
    ],
    PCoA2 = pcoa_males$points[
      metadata_males$sample,
      2
    ],
    male_morph = metadata_males$male_morph
  )
  
  
  write.csv(
    pcoa_males_df,
    file = file.path(
      out_dir,
      paste0(
        "PCoA.coordinates.males.",
        condition,
        ".",
        partition_label,
        ".csv"
      )
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
    geom_point(
      size = 4
    ) +
    geom_polygon(
      aes(
        group = male_morph,
        fill = male_morph,
        color = male_morph
      ),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    ) +
    geom_text(
      aes(
        label = sample
      ),
      vjust = -0.8,
      size = 3,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = male_morph_colors,
      drop = FALSE
    ) +
    scale_color_manual(
      values = male_morph_colors,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "Male-only PCoA: ",
        condition,
        " - ",
        partition_title,
        "\nVMR categories: ",
        category_title
      ),
      x = sprintf(
        "PCoA1 (%.1f%%)",
        male_pcoa1_percent
      ),
      y = sprintf(
        "PCoA2 (%.1f%%)",
        male_pcoa2_percent
      ),
      color = "Male morph"
    ) +
    theme_classic()
  
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        "PCoA.males.",
        condition,
        ".",
        partition_label,
        ".pdf"
      )
    ),
    plot = pcoa_males_plot,
    width = 7,
    height = 5
  )
  
  
  male_distance_matrix <- as.matrix(
    d_males
  )
  
  annotation_males <- metadata_males[
    ,
    "male_morph",
    drop = FALSE
  ]
  
  male_clustering <- hclust(
    d_males,
    method = method_h_clust
  )
  
  
  pheatmap(
    male_distance_matrix,
    cluster_rows = male_clustering,
    cluster_cols = male_clustering,
    annotation_row = annotation_males,
    annotation_col = annotation_males,
    color = color_palette,
    annotation_colors = annotation_colors_males,
    main = paste0(
      "Male-only Euclidean distances: ",
      condition,
      " - ",
      partition_title,
      "\nVMR categories: ",
      category_title
    ),
    border_color = NA,
    filename = file.path(
      out_dir,
      paste0(
        "Heatmap.male_sample_distances.",
        condition,
        ".",
        partition_label,
        ".pdf"
      )
    ),
    width = 7,
    height = 6
  )
  
  
  # Return partition VMR counts
  data.frame(
    partition = partition_label,
    vmr_categories = category_title,
    n_vmr_before_complete_filter = n_vmr_before_complete_filter,
    n_vmr_complete_12_samples = n_vmr_complete
  )
}


# Run both partitions -----------------------------------------------------------


partition_summary <- rbindlist(
  lapply(
    names(partition_data),
    function(partition_label) {
      run_partition_analysis(
        dt_partition = partition_data[[partition_label]],
        partition_label = partition_label,
        partition_title = partition_titles[[partition_label]]
      )
    }
  ),
  fill = TRUE
)


# Save number of VMRs used in each partition ----------------------------------

fwrite(
  partition_summary,
  file = file.path(
    out_dir_base,
    paste0(
      "VMR_counts.by_genomic_partition.",
      condition,
      ".tsv"
    )
  ),
  sep = "\t"
)

print(
  partition_summary
)

cat("\nOutputs saved to:\n", normalizePath(out_dir_base), "\n")
