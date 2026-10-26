#!/usr/bin/env Rscript
# DATE: 20260917
# AUTHOR: MZF
# SCRIPT: 06_heatmap_VMRcategories_Zscore_V10.R
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Visualize selected biological VMR categories while preserving the FIXED
#   whole-dataset sample dendrogram/column order from V09.
#
# CORE DATASET SEPARATION:
#
#   dt_all
#       Complete filtered VMR dataset. Used for:
#         - phenotype-association (r.g) calculation
#         - fixed whole-dataset sample-tree construction
#
#   dt_plot
#       Category-selected subset. Used ONLY for displayed heatmap rows.
#
#   Therefore VMR_CATEGORIES affects which VMRs are plotted, but DOES NOT
#   affect the sample dendrogram or heatmap column order.
#
# V10 ROW ORDER:
#
#   biological_category
#       -> indicator_phenotype
#           -> decreasing signed r_g
#
#   Category order follows VMR_CATEGORIES.
#   Indicator phenotype order follows PHENOTYPE_ORDER.
#   Within each indicator phenotype, r_g is ordered maximum -> minimum.
#   Genomic coordinates are used only as tie-breakers.
#
# INDICATOR PHENOTYPE:
#   indicator_phenotype is the phenotype with the strongest absolute r.g.
#   abs_r_g is retained for QC/output but is NOT used for ordering.
#
# r.g CALCULATION:
#   Default: strassoc() calculates observed phenotype-specific r.g directly
#   without permutations because r.g is used descriptively for annotation
#   and ordering.
#
#   Optional: USE_PERMUTATION_TEST <- TRUE uses multipatt().
#
# FIXED COLUMN TREE:
#   build_sample_tree() is run on dt_all BEFORE plotting.
#
#   For each genomic partition:
#     1) all VMRs in dt_all for that partition are used;
#     2) VMRs are row-Z-scored across samples;
#     3) the existing optional phenotype-consistent Z filter is applied;
#     4) Euclidean sample distances are calculated;
#     5) distances are divided by sqrt(number of retained VMRs);
#     6) complete-linkage clustering generates a fixed hclust object.
#
#   hc_cols_auto and hc_cols_x are then passed directly into plot_heatmap().
#   Sample clustering is NEVER recalculated from dt_plot.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(pheatmap)
  library(viridis)
  library(indicspecies)
  library(permute)
})

setwd("/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/data")

# Parameters -------------------------------------------------------------------

MIN_SAMPLE_RANGE <- 0.10

# Phenotype-consistent Z-score filter
# Set either parameter to 0 to disable this filter
MIN_ABS_Z <- 0
MIN_N_SAMPLES_Z <- 0

# Indicator Species Analysis
ISA_FUNC <- "r.g"

# FALSE = observed phenotype-specific r.g only using strassoc()
# TRUE  = multipatt() with permutation testing
USE_PERMUTATION_TEST <- FALSE

# Used only when USE_PERMUTATION_TEST = TRUE
MAX_ORDER <- 1
ALLOW_NEGATIVE <- TRUE
N_PERM <- 9999

set.seed(20260811)


# Samples to exclude
# EXCLUDE_SAMPLES <- c(
#   "pparfmem001",
#   "pparfmem002",
#   "pparfmem003"
# )
EXCLUDE_SAMPLES <- c()


# Categories displayed in heatmap
# Does NOT affect ISA or sample dendrogram
VMR_CATEGORIES <- c(
  "Sex-associated",
  # "Sex + morph",
  # "Cross-phenotype",
  "Morph-specific"
  # "Pairwise-specific",
  # "Other contrast combinations"
)


annotation_file <- "../plots/complexupset/vmr_real.contrastSharing.biologicalCategories.perVMR.tsv"
methylation_file <- "vmr_methylation/vmr.morph_and_sample_observed.tsv"
outdir <- "../plots/VMRcategory_heatmaps_Zscore_V10_fixedWholeTree_indicatorPhenotype_rG"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# Read and join ----------------------------------------------------------------

anno <- fread(annotation_file)
meth <- fread(methylation_file)

required_anno <- c(
  "cons_id",
  "biological_category",
  "computed_contrast_combo",
  "biological_intersection_order"
)

if (length(setdiff(required_anno, names(anno))))
  stop(
    "Missing annotation columns: ",
    paste(setdiff(required_anno, names(anno)), collapse = ", ")
  )

if (anyDuplicated(anno$cons_id) || anyDuplicated(meth$cons_id))
  stop("cons_id must be unique in both input tables.")

missing_ids <- setdiff(meth$cons_id, anno$cons_id)

if (length(missing_ids))
  stop(
    length(missing_ids),
    " methylation VMRs are absent from the UpSet annotation table."
  )

valid_categories <- unique(anno$biological_category)
bad_categories <- setdiff(VMR_CATEGORIES, valid_categories)

if (length(bad_categories))
  stop(
    "Unknown VMR categories: ",
    paste(bad_categories, collapse = ", ")
  )

anno <- anno[, .(
  cons_id,
  biological_category,
  computed_contrast_combo,
  biological_intersection_order
)]

# Do NOT subset VMR_CATEGORIES here
dt <- merge(
  meth,
  anno,
  by = "cons_id",
  all.x = TRUE,
  sort = FALSE
)


# Identify and optionally exclude samples --------------------------------------

all_sample_cols <- sort(
  grep("^meth_ppar.*mem[0-9]+$", names(dt), value = TRUE)
)

if (!length(all_sample_cols))
  stop("No per-sample methylation columns found.")

exclude_cols <- paste0("meth_", EXCLUDE_SAMPLES)
unknown_exclusions <- setdiff(exclude_cols, all_sample_cols)

if (length(unknown_exclusions))
  warning(
    "Requested excluded samples were not found: ",
    paste(sub("^meth_", "", unknown_exclusions), collapse = ", ")
  )

sample_cols <- setdiff(all_sample_cols, exclude_cols)

if (!length(sample_cols))
  stop("No per-sample methylation columns remain after sample exclusion.")

sample_labels <- sub("^meth_", "", sample_cols)

cat("Samples included: ", paste(sample_labels, collapse = ", "), "\n", sep = "")
cat(
  "Samples excluded: ",
  if (length(EXCLUDE_SAMPLES)) paste(EXCLUDE_SAMPLES, collapse = ", ") else "none",
  "\n\n",
  sep = ""
)


# Assign samples to phenotypes -------------------------------------------------

sample_phenotype <- ifelse(
  grepl("^pparf", sample_labels), "Female",
  ifelse(
    grepl("^pparim", sample_labels), "Immaculata",
    ifelse(
      grepl("^pparpm", sample_labels), "Parae",
      ifelse(
        grepl("^pparym", sample_labels), "Yellow",
        NA_character_
      )
    )
  )
)

if (anyNA(sample_phenotype))
  stop(
    "Could not assign phenotype to samples: ",
    paste(sample_labels[is.na(sample_phenotype)], collapse = ", ")
  )

names(sample_phenotype) <- sample_labels

phenotype_levels <- c(
  "Female",
  "Immaculata",
  "Parae",
  "Yellow"
)

phenotype_levels <- phenotype_levels[
  phenotype_levels %in% unique(sample_phenotype)
]

# Explicit order used inside each biological-category row block
PHENOTYPE_ORDER <- phenotype_levels

sample_metadata <- data.table(
  sample_col = sample_cols,
  phenotype = factor(
    sample_phenotype[sample_labels],
    levels = phenotype_levels
  )
)

phenotype_counts <- table(sample_metadata$phenotype)

cat("Sample phenotypes:\n")
print(
  data.frame(
    Sample = sample_labels,
    Phenotype = unname(sample_phenotype)
  )
)

cat("\nSamples per phenotype:\n")
print(phenotype_counts)
cat("\n")


if (
  MIN_N_SAMPLES_Z > 0 &&
  MIN_N_SAMPLES_Z > max(phenotype_counts)
)
  stop(
    "MIN_N_SAMPLES_Z = ",
    MIN_N_SAMPLES_Z,
    " exceeds the maximum number of included samples in any phenotype (",
    max(phenotype_counts),
    ")."
  )


# Completeness and absolute methylation-range filters --------------------------

mat_all <- as.matrix(dt[, ..sample_cols])
storage.mode(mat_all) <- "double"

dt[, complete_meth :=
     rowSums(is.na(mat_all) | !is.finite(mat_all)) == 0]

dt[, sample_range :=
     apply(
       mat_all,
       1,
       function(x) if (all(is.finite(x))) max(x) - min(x) else NA_real_
     )]

n_before <- nrow(dt)
n_missing <- dt[complete_meth == FALSE, .N]

n_lowrange <- dt[
  complete_meth &
    !is.na(sample_range) &
    sample_range < MIN_SAMPLE_RANGE,
  .N
]

dt <- dt[
  complete_meth &
    !is.na(sample_range) &
    sample_range >= MIN_SAMPLE_RANGE
]

cat("VMRs before completeness/range filter: ", n_before, "\n", sep = "")
cat("Removed for missing methylation:        ", n_missing, "\n", sep = "")
cat(
  "Removed for sample range < ",
  MIN_SAMPLE_RANGE,
  ":        ",
  n_lowrange,
  "\n",
  sep = ""
)
cat("VMRs retained for ISA:                  ", nrow(dt), "\n\n", sep = "")


# Indicator Species Analysis ---------------------------------------------------

methylation_matrix <- as.matrix(dt[, ..sample_cols])
storage.mode(methylation_matrix) <- "numeric"

# samples = rows; VMRs = columns
methylation_matrix <- t(methylation_matrix)

rownames(methylation_matrix) <- sample_metadata$sample_col
colnames(methylation_matrix) <- dt$cons_id

cat(
  "Calculating phenotype association statistics on ",
  ncol(methylation_matrix),
  " VMRs...\n",
  sep = ""
)


# Fast observed r.g calculation ------------------------------------------------

if (!USE_PERMUTATION_TEST) {
  
  cat(
    "Using strassoc(): observed phenotype-specific r.g only; ",
    "no permutations.\n",
    sep = ""
  )
  
  cluster_numeric <- as.integer(sample_metadata$phenotype)
  
  r_g_mat <- indicspecies::strassoc(
    X = methylation_matrix,
    cluster = cluster_numeric,
    func = ISA_FUNC
  )
  
  if (!is.matrix(r_g_mat))
    stop("strassoc() did not return the expected association matrix.")
  
  if (nrow(r_g_mat) != nrow(dt))
    stop(
      "Unexpected number of VMR rows returned by strassoc(): ",
      nrow(r_g_mat),
      " instead of ",
      nrow(dt),
      "."
    )
  
  if (ncol(r_g_mat) != length(phenotype_levels))
    stop(
      "Unexpected number of phenotype columns returned by strassoc(): ",
      ncol(r_g_mat),
      " instead of ",
      length(phenotype_levels),
      "."
    )
  
  # Force explicit labels/order matching the phenotype factor levels.
  rownames(r_g_mat) <- dt$cons_id
  colnames(r_g_mat) <- phenotype_levels
  
  # Indicator phenotype = strongest association by magnitude.
  indicator_group_index <- max.col(
    abs(r_g_mat),
    ties.method = "first"
  )
  
  isa_stats <- data.table(
    cons_id = rownames(r_g_mat),
    indicator_phenotype = phenotype_levels[indicator_group_index],
    r_g = r_g_mat[
      cbind(
        seq_len(nrow(r_g_mat)),
        indicator_group_index
      )
    ],
    p_value = NA_real_,
    p_adj_BH = NA_real_
  )
}


# Optional multipatt permutation analysis --------------------------------------

if (USE_PERMUTATION_TEST) {
  
  cat(
    "Using multipatt() with ",
    N_PERM,
    " permutations...\n",
    sep = ""
  )
  
  isa_result <- multipatt(
    x = methylation_matrix,
    cluster = sample_metadata$phenotype,
    func = ISA_FUNC,
    max.order = MAX_ORDER,
    allow.negative = ALLOW_NEGATIVE,
    control = how(
      nperm = N_PERM
    )
  )
  
  isa_stats <- as.data.table(
    isa_result$sign,
    keep.rownames = "cons_id"
  )
  
  group_columns_named <- paste0("s.", phenotype_levels)
  group_columns_numeric <- paste0("s.", seq_along(phenotype_levels))
  
  if (all(group_columns_named %in% names(isa_stats))) {
    group_columns <- group_columns_named
  } else if (all(group_columns_numeric %in% names(isa_stats))) {
    group_columns <- group_columns_numeric
  } else {
    stop("Could not identify phenotype-membership columns in ISA output.")
  }
  
  indicator_group_index <- max.col(
    as.matrix(isa_stats[, ..group_columns]),
    ties.method = "first"
  )
  
  isa_stats[
    ,
    indicator_phenotype := phenotype_levels[indicator_group_index]
  ]
  
  setnames(
    isa_stats,
    old = c("stat", "p.value"),
    new = c("r_g", "p_value")
  )
  
  isa_stats[
    ,
    p_adj_BH := p.adjust(
      p_value,
      method = "BH"
    )
  ]
  
  isa_stats <- isa_stats[
    ,
    .(
      cons_id,
      indicator_phenotype,
      r_g,
      p_value,
      p_adj_BH
    )
  ]
}


# Derived r.g variables --------------------------------------------------------

isa_stats[, abs_r_g := abs(r_g)]

isa_stats[
  ,
  direction := fifelse(
    r_g > 0,
    "higher_in_indicator",
    fifelse(
      r_g < 0,
      "lower_in_indicator",
      "none"
    )
  )
]

cat("Phenotype association calculation completed.\n\n")


# Match ISA output to VMR order ------------------------------------------------

isa_stats <- isa_stats[
  match(
    dt$cons_id,
    isa_stats$cons_id
  )
]

if (any(dt$cons_id != isa_stats$cons_id))
  stop("ISA results could not be aligned to the VMR input order.")


# Append ISA statistics --------------------------------------------------------

dt[
  ,
  `:=`(
    indicator_phenotype = isa_stats$indicator_phenotype,
    r_g = isa_stats$r_g,
    abs_r_g = isa_stats$abs_r_g,
    direction = isa_stats$direction,
    p_value = isa_stats$p_value,
    p_adj_BH = isa_stats$p_adj_BH
  )
]

cat(
  "Indicator Species Analysis completed.\n",
  "r.g range: ",
  min(dt$r_g, na.rm = TRUE),
  " to ",
  max(dt$r_g, na.rm = TRUE),
  "\n",
  "|r.g| range: ",
  min(dt$abs_r_g, na.rm = TRUE),
  " to ",
  max(dt$abs_r_g, na.rm = TRUE),
  "\n\n",
  sep = ""
)


# Define full and plotting datasets --------------------------------------------

dt_all <- copy(dt)

dt_plot <- dt_all[
  biological_category %in% VMR_CATEGORIES
]

if (!nrow(dt_plot))
  stop("No VMRs remain after applying VMR_CATEGORIES.")

cat("VMRs in full ISA-annotated dataset: ", nrow(dt_all), "\n", sep = "")
cat("VMRs selected for heatmap:           ", nrow(dt_plot), "\n\n", sep = "")


# Output tag -------------------------------------------------------------------

tag <- paste0(
  paste(
    gsub("[^A-Za-z0-9]+", "_", VMR_CATEGORIES),
    collapse = "-"
  ),
  ".rangeGE_", MIN_SAMPLE_RANGE,
  ".zGE_", MIN_ABS_Z,
  ".samePhenoNGE_", MIN_N_SAMPLES_Z,
  ".rowOrder_category_indicatorPhenotype_rGdesc"
)


# Export ISA-annotated datasets ------------------------------------------------

fwrite(
  dt_all,
  file.path(
    outdir,
    paste0(
      "VMRs.fullFiltered.ISA.rangeGE_",
      MIN_SAMPLE_RANGE,
      ".tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  dt_plot,
  file.path(
    outdir,
    paste0(
      "VMRs.selected.preZ.",
      tag,
      ".tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# QC ---------------------------------------------------------------------------

cat(
  "Selected heatmap categories: ",
  paste(VMR_CATEGORIES, collapse = ", "),
  "\n\n",
  sep = ""
)

cat("Indicator phenotype counts:\n")

print(
  dt_all[
    ,
    .N,
    by = indicator_phenotype
  ][
    order(indicator_phenotype)
  ]
)

cat("\n")


# Helper: phenotype-consistent Z-score filter ----------------------------------

evaluate_z_filter <- function(mat_z) {
  
  if (MIN_ABS_Z <= 0 || MIN_N_SAMPLES_Z <= 0) {
    return(
      data.table(
        z_pass = rep(TRUE, nrow(mat_z)),
        z_pass_phenotypes = rep("Filter_disabled", nrow(mat_z)),
        z_pass_pattern = rep("Filter_disabled", nrow(mat_z))
      )
    )
  }
  
  phenotypes <- unique(sample_phenotype)
  
  z_pass_info <- lapply(
    seq_len(nrow(mat_z)),
    function(i) {
      
      z <- mat_z[i, ]
      passing_phenotypes <- character()
      passing_patterns <- character()
      
      for (pheno in phenotypes) {
        
        idx <- sample_phenotype[colnames(mat_z)] == pheno
        
        n_high <- sum(z[idx] >= MIN_ABS_Z)
        n_low <- sum(z[idx] <= -MIN_ABS_Z)
        
        if (n_high >= MIN_N_SAMPLES_Z) {
          passing_phenotypes <- c(passing_phenotypes, pheno)
          passing_patterns <- c(
            passing_patterns,
            paste0(pheno, ":high(", n_high, ")")
          )
        }
        
        if (n_low >= MIN_N_SAMPLES_Z) {
          passing_phenotypes <- c(passing_phenotypes, pheno)
          passing_patterns <- c(
            passing_patterns,
            paste0(pheno, ":low(", n_low, ")")
          )
        }
      }
      
      list(
        keep = length(passing_patterns) > 0,
        phenotypes = paste(unique(passing_phenotypes), collapse = ";"),
        pattern = paste(passing_patterns, collapse = ";")
      )
    }
  )
  
  data.table(
    z_pass = vapply(z_pass_info, function(x) x$keep, logical(1)),
    z_pass_phenotypes = vapply(z_pass_info, function(x) x$phenotypes, character(1)),
    z_pass_pattern = vapply(z_pass_info, function(x) x$pattern, character(1))
  )
}


# Build fixed whole-dataset sample tree ----------------------------------------

build_sample_tree <- function(chrom_type, label, prefix) {
  
  tree_sub <- dt_all[
    chromosome_type == chrom_type
  ]
  
  setorder(tree_sub, chr, start, end)
  
  if (nrow(tree_sub) < 2)
    stop(
      "Fewer than 2 VMRs available for sample-tree construction in ",
      label,
      "."
    )
  
  
  # Raw methylation matrix
  
  tree_raw <- as.matrix(tree_sub[, ..sample_cols])
  storage.mode(tree_raw) <- "double"
  
  rownames(tree_raw) <- tree_sub$cons_id
  colnames(tree_raw) <- sample_labels
  
  if (anyNA(tree_raw) || any(!is.finite(tree_raw)))
    stop(
      "Missing/non-finite methylation values in sample-tree dataset for ",
      label,
      "."
    )
  
  
  # Row Z-score
  
  tree_z <- t(scale(t(tree_raw)))
  
  if (anyNA(tree_z) || any(!is.finite(tree_z)))
    stop(
      "Non-finite values produced during whole-dataset Z-scoring in ",
      label,
      "."
    )
  
  
  # Optional phenotype-consistent Z filter
  
  tree_z_info <- evaluate_z_filter(tree_z)
  
  tree_sub[, z_pass := tree_z_info$z_pass]
  tree_sub[, z_pass_phenotypes := tree_z_info$z_pass_phenotypes]
  tree_sub[, z_pass_pattern := tree_z_info$z_pass_pattern]
  
  keep_tree <- tree_z_info$z_pass
  
  tree_sub <- tree_sub[keep_tree]
  tree_raw <- tree_raw[keep_tree, , drop = FALSE]
  tree_z <- tree_z[keep_tree, , drop = FALSE]
  
  if (nrow(tree_z) < 2)
    stop(
      "Fewer than 2 VMRs remain for sample-tree construction in ",
      label,
      "."
    )
  
  
  # Normalized Euclidean sample distances
  
  d_samples <- dist(t(tree_z))
  d_samples_norm <- d_samples / sqrt(nrow(tree_z))
  
  
  # Complete-linkage sample clustering
  
  hc_cols <- hclust(
    d_samples_norm,
    method = "complete"
  )
  
  cat(
    label,
    ": FIXED sample tree based on ",
    nrow(tree_z),
    " whole-dataset VMRs using normalized Euclidean distance\n",
    sep = ""
  )
  
  
  # Export VMRs used for sample tree
  
  fwrite(
    tree_sub,
    file.path(
      outdir,
      paste0(
        prefix,
        ".sampleTree.VMRs.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  
  
  # Standalone sample dendrogram
  
  pdf(
    file.path(
      outdir,
      paste0(
        prefix,
        ".sampleTree.normalizedEuclidean.pdf"
      )
    ),
    width = 7,
    height = 5
  )
  
  par(mar = c(8, 5, 3, 2))
  
  plot(
    hc_cols,
    hang = -1,
    main = paste0(
      label,
      "\nSample clustering — full filtered VMR dataset"
    ),
    xlab = "",
    sub = "",
    ylab = "Normalized Euclidean distance (Z-score units)",
    cex = 0.9
  )
  
  dev.off()
  
  # Export fixed sample/column order ------------------------------------------
  
  fixed_sample_order <- data.table(
    dendrogram_order = seq_along(hc_cols$order),
    sample_col = hc_cols$labels[hc_cols$order]
  )
  
  fixed_sample_order[
    ,
    sample := sub("^meth_", "", sample_col)
  ]
  
  fixed_sample_order[
    ,
    phenotype := sample_phenotype[sample]
  ]
  
  fwrite(
    fixed_sample_order,
    file.path(
      outdir,
      paste0(
        prefix,
        ".sampleTree.fixedColumnOrder.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  
  return(hc_cols)
}


# Construct fixed reference sample trees ---------------------------------------

hc_cols_auto <- build_sample_tree(
  "Autosome",
  "Autosomes",
  "AUTOSOME.fullDataset"
)

hc_cols_x <- build_sample_tree(
  "Sex_Ch",
  "X chromosome (Parae_12)",
  "XCHROM.fullDataset"
)


# Heatmap ----------------------------------------------------------------------

colors <- inferno(101)
breaks <- seq(0, 1, length.out = 102)


plot_heatmap <- function(chrom_type, label, prefix, hc_cols) {
  
  sub <- dt_plot[
    chromosome_type == chrom_type
  ]
  
  setorder(sub, chr, start, end)
  
  if (nrow(sub) < 2) {
    warning(
      "Skipping ",
      label,
      ": fewer than 2 category-selected VMRs."
    )
    return(invisible(NULL))
  }
  
  
  # Raw methylation matrix
  
  mat_raw <- as.matrix(sub[, ..sample_cols])
  storage.mode(mat_raw) <- "double"
  
  rownames(mat_raw) <- sub$cons_id
  colnames(mat_raw) <- sample_labels
  
  if (anyNA(mat_raw) || any(!is.finite(mat_raw)))
    stop(
      "Missing/non-finite methylation values found in ",
      label,
      "."
    )
  
  
  # Row Z-score
  
  mat_z <- t(scale(t(mat_raw)))
  
  if (anyNA(mat_z) || any(!is.finite(mat_z)))
    stop(
      "Non-finite values produced during row scaling in ",
      label,
      "."
    )
  
  
  # Optional phenotype-consistent Z-score filter
  
  z_info <- evaluate_z_filter(mat_z)
  
  sub[, z_pass := z_info$z_pass]
  sub[, z_pass_phenotypes := z_info$z_pass_phenotypes]
  sub[, z_pass_pattern := z_info$z_pass_pattern]
  
  keep_z <- z_info$z_pass
  
  sub <- sub[keep_z]
  mat_raw <- mat_raw[keep_z, , drop = FALSE]
  mat_z <- mat_z[keep_z, , drop = FALSE]
  
  cat(
    label,
    ": ",
    nrow(mat_raw),
    " VMRs retained for heatmap\n",
    sep = ""
  )
  
  if (nrow(mat_raw) < 2) {
    warning(
      "Skipping ",
      label,
      ": fewer than 2 VMRs remain after Z-score filtering."
    )
    return(invisible(NULL))
  }
  
  
  # ---------------------------------------------------------------------------
  # FIXED WHOLE-DATASET COLUMN TREE
  #
  # Sample clustering is NOT recalculated here.
  # hc_cols was built from dt_all in build_sample_tree(), before category
  # subsetting. VMR_CATEGORIES therefore changes ROWS only.
  # ---------------------------------------------------------------------------
  
  
  # Organize rows: category -> indicator phenotype -> signed r.g ---------------
  #
  # This is ROW ordering only.
  #
  # hc_cols was already calculated from dt_all by build_sample_tree().
  # Nothing in this section can alter the fixed whole-dataset column tree.
  
  category_order <- integer(0)
  within_category_order <- integer(nrow(sub))
  
  for (category in VMR_CATEGORIES) {
    
    idx <- which(
      sub$biological_category == category
    )
    
    if (!length(idx))
      next
    
    phenotype_rank <- match(
      sub$indicator_phenotype[idx],
      PHENOTYPE_ORDER
    )
    
    if (anyNA(phenotype_rank))
      stop(
        "Unknown indicator phenotype(s) in ",
        category,
        ": ",
        paste(
          unique(
            sub$indicator_phenotype[idx][is.na(phenotype_rank)]
          ),
          collapse = ", "
        )
      )
    
    # Nested ordering:
    #   indicator phenotype -> decreasing signed r.g -> genomic tie-breaker
    ordered_idx <- idx[
      order(
        phenotype_rank,
        -sub$r_g[idx],
        sub$chr[idx],
        sub$start[idx],
        sub$end[idx],
        na.last = TRUE
      )
    ]
    
    within_category_order[
      ordered_idx
    ] <- seq_along(ordered_idx)
    
    category_order <- c(
      category_order,
      ordered_idx
    )
  }
  
  
  if (length(category_order) != nrow(sub))
    stop(
      "Not all retained VMRs were assigned to a biological-category block in ",
      label,
      "."
    )
  
  
  # Save within-category ordering
  
  sub[
    ,
    within_category_order := within_category_order
  ]
  
  
  # Reorder table and matrices
  
  sub <- sub[category_order]
  mat_raw <- mat_raw[category_order, , drop = FALSE]
  mat_z <- mat_z[category_order, , drop = FALSE]
  
  sub[, heatmap_row_order := seq_len(.N)]
  
  sub[
    ,
    within_indicator_order := seq_len(.N),
    by = .(
      biological_category,
      indicator_phenotype
    )
  ]
  
  
  # Gaps between biological categories
  
  category_runs <- rle(
    sub$biological_category
  )
  
  gaps_row <- cumsum(
    category_runs$lengths
  )
  
  if (length(gaps_row) > 1) {
    gaps_row <- gaps_row[-length(gaps_row)]
  } else {
    gaps_row <- NULL
  }
  
  
  # Row annotation
  
  annotation_row <- data.frame(
    Category = factor(
      sub$biological_category,
      levels = VMR_CATEGORIES
    ),
    row.names = sub$cons_id
  )
  
  
  # Export final ordered VMR table
  
  fwrite(
    sub,
    file.path(
      outdir,
      paste0(
        prefix,
        ".VMRs.final.category_indicatorPhenotype_rGordered.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  
  
  # Export compact row-order table
  
  row_order_table <- sub[
    ,
    .(
      heatmap_row_order,
      within_category_order,
      within_indicator_order,
      chr,
      start,
      end,
      cons_id,
      biological_category,
      computed_contrast_combo,
      biological_intersection_order,
      sample_range,
      
      indicator_phenotype,
      r_g,
      abs_r_g,
      direction,
      p_value,
      p_adj_BH,
      
      z_pass_phenotypes,
      z_pass_pattern
    )
  ]
  
  fwrite(
    row_order_table,
    file.path(
      outdir,
      paste0(
        prefix,
        ".rowOrder.category_indicatorPhenotype_rG.tsv"
      )
    ),
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  
  
  # Heatmap
  
  pdf(
    file.path(
      outdir,
      paste0(
        prefix,
        ".heatmap.pdf"
      )
    ),
    width = 6.5,
    height = 8
  )
  
  if (!setequal(colnames(mat_raw), hc_cols$labels))
    stop(
      "Heatmap sample columns do not match the samples in the fixed ",
      "whole-dataset tree for ",
      label,
      "."
    )
  
  pheatmap(
    mat_raw,
    cluster_rows = FALSE,
    cluster_cols = hc_cols,
    gaps_row = gaps_row,
    annotation_row = annotation_row,
    color = colors,
    breaks = breaks,
    na_col = "black",
    border_color = NA,
    show_rownames = FALSE,
    fontsize = 10,
    fontsize_col = 8,
    angle_col = 90,
    treeheight_row = 0,
    treeheight_col = 30,
    main = paste0(
      label,
      "\n",
      paste(
        VMR_CATEGORIES,
        collapse = " + "
      )
    )
  )
  
  dev.off()
}


# Run --------------------------------------------------------------------------

plot_heatmap(
  "Autosome",
  "Autosomes",
  paste0(
    "AUTOSOME.",
    tag
  ),
  hc_cols_auto
)

plot_heatmap(
  "Sex_Ch",
  "X chromosome (Parae_12)",
  paste0(
    "XCHROM.",
    tag
  ),
  hc_cols_x
)


# Output location ---------------------------------------------------------------

cat(
  "\nWrote outputs to: ",
  normalizePath(
    outdir,
    winslash = "/",
    mustWork = FALSE
  ),
  "\n",
  sep = ""
)