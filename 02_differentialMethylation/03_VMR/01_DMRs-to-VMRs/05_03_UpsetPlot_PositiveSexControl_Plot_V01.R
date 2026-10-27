#!/usr/bin/env Rscript
# DATE: 20261004
# AUTHOR: MZF
# SCRIPT: 05_03_UpsetPlot_PositiveSexControl_Plot_V01.R
# VERSION: 02
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Build UpSet plots showing how pooled shuffled VMRs are shared across the
#   six shuffled control contrasts.
#
#   The shuffled contrasts are globally mapped to the Sex
#   biological contrasts used in the real-data analysis:
#
#     Shuffled sex-control contrasts:
#           - f_vs_sM_01 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_02 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_03 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_04 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_05 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_06 ->  f_vs_i, f_vs_p, f_vs_y
#
#
#   The same structural contrast-sharing categories used for the real VMRs
#   are applied to the shuffled controls and displayed in the same order:
#
#     1) Sex-associated
#     2) Sex + morph
#     3) Cross-phenotype
#     4) Morph-specific
#     5) Pairwise-specific
#
#   IMPORTANT:
#   These category names are retained for direct comparability with the
#   real-data analysis. For shuffled VMRs, they describe shuffled analogues
#   of the biological contrast-sharing patterns and should NOT be interpreted
#   as evidence of biological sex or morph association.
#
#   Within each category, shuffled contrasts follow the input-table order:
#
#     f_vs_sM_01
#     f_vs_sM_02
#     f_vs_sM_03
#     f_vs_sM_04
#     f_vs_sM_05
#     f_vs_sM_06
#
#   Save one PDF per genomic partition:
#     1) Whole genome
#     2) Autosomes
#     3) Sex chromosome
#
#   IMPORTANT:
#   Only the 16 contrast combinations belonging to the five predefined
#   categories are shown in the UpSet plots.
#
#   VMRs belonging to other contrast combinations are retained in the input
#   data but are not displayed. A summary table reports how many VMRs belong
#   to each predefined category versus "Other contrast combinations".
#
# INPUT:
#   vmr_shuffled_sex.contrastSharing.binary.tsv
#
# OUTPUT:
#   vmr_shuffled_sex.contrastSharing.WG.upset.pdf
#   vmr_shuffled_sex.contrastSharing.AUTO.upset.pdf
#   vmr_shuffled_sex.contrastSharing.SEXCH.upset.pdf
#   vmr_shuffled_sex.contrastSharing.biologicalCategories.summary.tsv
#   vmr_shuffled_sex.contrastSharing.setSizes.summary.txt
#
# USAGE:
#   Rscript 05_03_UpsetPlot_PositiveSexControl_Plot_V01.R
#
#   or
#
#   Rscript 05_03_UpsetPlot_PositiveSexControl_Plot_V01.R input.tsv output_prefix
#
# REQUIRED R PACKAGES:
#   data.table
#   ggplot2
#   ComplexUpset
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%


# Libraries --------------------------------------------------------------------

# devtools::install_github('krassowski/complex-upset')

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexUpset)
})


# Command-line arguments -------------------------------------------------------

args <- commandArgs(
  trailingOnly = TRUE
)


# Defaults ---------------------------------------------------------------------

input_file <- "vmr_shuffled_sex.contrastSharing.binary.tsv"

outdir <- "../plots/complexupset/"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

output_pref <- file.path(
  outdir,
  "vmr_shuffled_sex.contrastSharing"
)

if (length(args) >= 1) {
  input_file <- args[1]
}

if (length(args) >= 2) {
  output_pref <- args[2]
}


# Output files -----------------------------------------------------------------

out_pdf_wg <- paste0(
  output_pref,
  ".WG.upset.pdf"
)

out_pdf_auto <- paste0(
  output_pref,
  ".AUTO.upset.pdf"
)

out_pdf_sex <- paste0(
  output_pref,
  ".SEXCH.upset.pdf"
)

out_category_summary <- paste0(
  output_pref,
  ".biologicalCategories.summary.tsv"
)

out_set_size_summary <- paste0(
  output_pref,
  ".setSizes.summary.txt"
)


# Read input -------------------------------------------------------------------

dt <- fread(
  input_file
)

message(
  "Loaded input file: ",
  input_file
)

message(
  "Rows: ",
  nrow(dt),
  " | Cols: ",
  ncol(dt)
)


# Check required metadata columns ---------------------------------------------

required_cols <- c(
  "cons_id",
  "chromosome_type",
  "chr"
)

missing_cols <- setdiff(
  required_cols,
  names(dt)
)

if (
  length(missing_cols) > 0
) {
  
  stop(
    "Missing required columns: ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}


# Define contrast columns ------------------------------------------------------
#
# IMPORTANT:
# This order is used throughout the script and directly mirrors the order
# of the corresponding biological contrasts in the real-data analysis.
#
#     Shuffled sex-control contrasts:
#           - f_vs_sM_01 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_02 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_03 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_04 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_05 ->  f_vs_i, f_vs_p, f_vs_y
#           - f_vs_sM_06 ->  f_vs_i, f_vs_p, f_vs_y
#
# The _01/_02/_03 suffixes are retained because each shuffled comparison
# represents an independent set in the UpSet analysis.

contrast_cols <- c(
  "f_vs_sM_01",
  "f_vs_sM_02",
  "f_vs_sM_03",
  "f_vs_sM_04",
  "f_vs_sM_05",
  "f_vs_sM_06"
)


# Check that all six shuffled contrasts exist ---------------------------------

missing_contrasts <- setdiff(
  contrast_cols,
  names(dt)
)

if (
  length(missing_contrasts) > 0
) {
  
  stop(
    "Missing required shuffled contrast columns: ",
    paste(
      missing_contrasts,
      collapse = ", "
    )
  )
}


# Validate binary contrast columns --------------------------------------------

for (
  cc in contrast_cols
) {
  
  values_i <- unique(
    dt[[cc]]
  )
  
  values_i <- values_i[
    !is.na(values_i)
  ]
  
  if (
    !all(
      values_i %in% c(
        0,
        1
      )
    )
  ) {
    
    stop(
      "Contrast column ",
      cc,
      " contains values other than 0/1."
    )
  }
  
  if (
    anyNA(
      dt[[cc]]
    )
  ) {
    
    stop(
      "Contrast column ",
      cc,
      " contains NA values."
    )
  }
}

message(
  "Contrast order: ",
  paste(
    contrast_cols,
    collapse = ", "
  )
)


# Define shuffled analogues of biological intersection categories -------------
#
# The list below simultaneously defines:
#
#   1) which intersections are shown
#   2) the order in which they are shown
#
# The category labels are intentionally identical to those used for the
# real-data analysis so category frequencies can be compared directly.
#
# COLUMN BLOCKS IN THE FINAL UPSET:
#
#   Column  1       = Sex-associated
#   Columns 2-4     = Sex + morph
#   Columns 5-7     = Cross-phenotype
#   Columns 8-10    = Morph-specific
#   Columns 11-16   = Pairwise-specific
#


# ------------------------------------------------------------------------------
# 1. Sex-associated VMRs
#
# Shuffled analogue of:
#   all three female-male comparisons
#   no male-male comparison
#
# Therefore:
#   all three sA contrasts
#   no sM contrast
# ------------------------------------------------------------------------------

sex_associated <- list(
  
  c(
    "f_vs_sM_01",
    "f_vs_sM_02",
    "f_vs_sM_03",
    "f_vs_sM_04",
    "f_vs_sM_05",
    "f_vs_sM_06"
  )
  
)


# # ------------------------------------------------------------------------------
# # 2. Sex + morph VMRs
# #
# # Shuffled analogue of:
# #   all three female-male comparisons plus exactly two male-male comparisons.
# #
# # Therefore:
# #   all three sA contrasts plus exactly two sM contrasts.
# # ------------------------------------------------------------------------------
# 
# sex_plus_morph <- list(
#   
#   c(
#     "sA_vs_sA_01",
#     "sA_vs_sA_02",
#     "sA_vs_sA_03",
#     "sM_vs_sM_01",
#     "sM_vs_sM_02",
#     "sM_vs_sM_01",
#     "sM_vs_sM_02"
#   ),
#   
#   c(
#     "sA_vs_sA_01",
#     "sA_vs_sA_02",
#     "sA_vs_sA_03",
#     "sM_vs_sM_01",
#     "sM_vs_sM_03"
#   ),
#   
#   c(
#     "sA_vs_sA_01",
#     "sA_vs_sA_02",
#     "sA_vs_sA_03",
#     "sM_vs_sM_02",
#     "sM_vs_sM_03"
#   )
#   
# )


# ------------------------------------------------------------------------------
# 3. Cross-phenotype VMRs
#
# Shuffled analogues of the three four-contrast patterns defining candidate
# two-versus-two phenotype splits in the real-data analysis.
#
# Original pattern 1:
#   F + I  |  P + Y
#
#   f_vs_p + f_vs_y + i_vs_y + p_vs_i
#
# Shuffled analogue:
#   sA_vs_sA_02 + sA_vs_sA_03 + sM_vs_sM_01 + sM_vs_sM_02
#
#
# Original pattern 2:
#   F + P  |  I + Y
#
#   f_vs_i + f_vs_y + p_vs_i + y_vs_p
#
# Shuffled analogue:
#   sA_vs_sA_01 + sA_vs_sA_03 + sM_vs_sM_02 + sM_vs_sM_03
#
#
# Original pattern 3:
#   F + Y  |  I + P
#
#   f_vs_i + f_vs_p + i_vs_y + y_vs_p
#
# Shuffled analogue:
#   sA_vs_sA_01 + sA_vs_sA_02 + sM_vs_sM_01 + sM_vs_sM_03
# ------------------------------------------------------------------------------

# cross_phenotype <- list(
#   
#   # Shuffled analogue of F + I versus P + Y
#   c(
#     "sA_vs_sA_02",
#     "sA_vs_sA_03",
#     "sM_vs_sM_01",
#     "sM_vs_sM_02"
#   ),
#   
#   # Shuffled analogue of F + P versus I + Y
#   c(
#     "sA_vs_sA_01",
#     "sA_vs_sA_03",
#     "sM_vs_sM_02",
#     "sM_vs_sM_03"
#   ),
#   
#   # Shuffled analogue of F + Y versus I + P
#   c(
#     "sA_vs_sA_01",
#     "sA_vs_sA_02",
#     "sM_vs_sM_01",
#     "sM_vs_sM_03"
#   )
#   
# )


# ------------------------------------------------------------------------------
# 4. Morph-specific VMRs
#
# Shuffled analogue of:
#   one female-male contrast plus the two male-male contrasts involving
#   the same focal male morph.
#
# Pattern 1:
#   Immaculata-specific analogue
#
# Pattern 2:
#   Parae-specific analogue
#
# Pattern 3:
#   Yellow-specific analogue
# ------------------------------------------------------------------------------

# morph_specific <- list(
#   
#   # Immaculata-specific analogue
#   c(
#     "sA_vs_sA_01",
#     "sM_vs_sM_01",
#     "sM_vs_sM_02"
#   ),
#   
#   # Parae-specific analogue
#   c(
#     "sA_vs_sA_02",
#     "sM_vs_sM_02",
#     "sM_vs_sM_03"
#   ),
#   
#   # Yellow-specific analogue
#   c(
#     "sA_vs_sA_03",
#     "sM_vs_sM_01",
#     "sM_vs_sM_03"
#   )
#   
# )


# ------------------------------------------------------------------------------
# 5. Pairwise-specific VMRs
#
# Exactly one shuffled control contrast.
# ------------------------------------------------------------------------------

pairwise_specific <- list(
  c("f_vs_sM_01"),
  c("f_vs_sM_02"),
  c("f_vs_sM_03"),
  c("f_vs_sM_04"),
  c("f_vs_sM_05"),
  c("f_vs_sM_06")
)


# Combine intersections in biological order -----------------------------------

intersection_order <- c(
  sex_associated,
  # sex_plus_morph,
  # cross_phenotype,
  # morph_specific,
  pairwise_specific
)

if (
  length(intersection_order) != 7
) {
  
  stop(
    "Expected 16 shuffled analogue intersections, but found ",
    length(intersection_order),
    "."
  )
}

message(
  "Number of shuffled analogue intersections to plot: ",
  length(intersection_order)
)


# Build lookup table describing the 16 intersections ---------------------------

intersection_key <- data.table(
  
  intersection_order = seq_along(
    intersection_order
  ),
  
  biological_category = c(
    rep(
      "Sex-associated",
      length(sex_associated)
    ),
    # rep(
    #   "Sex + morph",
    #   length(sex_plus_morph)
    # ),
    # rep(
    #   "Cross-phenotype",
    #   length(cross_phenotype)
    # ),
    # rep(
    #   "Morph-specific",
    #   length(morph_specific)
    # ),
    rep(
      "Pairwise-specific",
      length(pairwise_specific)
    )
  ),
  
  contrast_combo = vapply(
    intersection_order,
    function(x) {
      paste(
        x,
        collapse = "&"
      )
    },
    character(1)
  )
  
)

message("")
message(
  "Shuffled analogue intersection order:"
)

for (
  i in seq_len(
    nrow(intersection_key)
  )
) {
  
  message(
    sprintf(
      "%2d  %-20s  %s",
      intersection_key$intersection_order[i],
      intersection_key$biological_category[i],
      intersection_key$contrast_combo[i]
    )
  )
}


# Compute contrast combination directly from binary columns --------------------
#
# This avoids relying on the existing contrast_combo column and guarantees
# that combinations are written according to contrast_cols order.

dt[
  ,
  computed_contrast_combo := apply(
    .SD,
    1,
    function(x) {
      
      x <- as.integer(
        x
      )
      
      present <- contrast_cols[
        x == 1L
      ]
      
      paste(
        present,
        collapse = "&"
      )
    }
  ),
  .SDcols = contrast_cols
]


# Assign each VMR to one of the five categories --------------------------------
#
# The biological category names are intentionally retained from the real-data
# analysis. Here they identify shuffled analogues of those exact sharing
# patterns.

combo_to_category <- setNames(
  intersection_key$biological_category,
  intersection_key$contrast_combo
)

dt[
  ,
  biological_category :=
    unname(
      combo_to_category[
        computed_contrast_combo
      ]
    )
]

dt[
  is.na(
    biological_category
  ),
  biological_category :=
    "Other contrast combinations"
]

category_levels <- c(
  "Sex-associated",
  # "Sex + morph",
  # "Cross-phenotype",
  # "Morph-specific",
  "Pairwise-specific",
  "Other contrast combinations"
)


# Split datasets by genomic partition -----------------------------------------

dt_wg <- copy(
  dt
)

dt_auto <- dt[
  chromosome_type %in% c(
    "Autosome",
    "AUTO",
    "autosome"
  )
]

if (
  nrow(dt_auto) == 0
) {
  
  dt_auto <- dt[
    grepl(
      "auto",
      chromosome_type,
      ignore.case = TRUE
    )
  ]
}

dt_sex <- dt[
  chromosome_type %in% c(
    "Sex_Ch",
    "SEXCH",
    "sex_ch",
    "sex",
    "SexChr"
  ) |
    chr == "Parae_12"
]

if (
  nrow(dt_sex) == 0
) {
  
  dt_sex <- dt[
    grepl(
      "sex",
      chromosome_type,
      ignore.case = TRUE
    ) |
      chr == "Parae_12"
  ]
}

message("")

message(
  "dt_wg rows: ",
  nrow(dt_wg)
)

message(
  "dt_auto rows: ",
  nrow(dt_auto)
)

message(
  "dt_sex rows: ",
  nrow(dt_sex)
)


# Summarize biological categories ---------------------------------------------

summarize_categories <- function(
    data_i,
    partition_i
) {
  
  counts_i <- data_i[
    ,
    .N,
    by = biological_category
  ]
  
  template_i <- data.table(
    biological_category = category_levels
  )
  
  summary_i <- merge(
    template_i,
    counts_i,
    by = "biological_category",
    all.x = TRUE,
    sort = FALSE
  )
  
  summary_i[
    is.na(N),
    N := 0L
  ]
  
  summary_i[
    ,
    partition := partition_i
  ]
  
  summary_i[
    ,
    category_order := match(
      biological_category,
      category_levels
    )
  ]
  
  summary_i[
    ,
    total_vmrs := nrow(
      data_i
    )
  ]
  
  summary_i[
    ,
    percentage :=
      if (
        nrow(data_i) > 0
      ) {
        100 * N / nrow(data_i)
      } else {
        0
      }
  ]
  
  setorder(
    summary_i,
    category_order
  )
  
  summary_i[
    ,
    c(
      "partition",
      "category_order",
      "biological_category",
      "N",
      "percentage",
      "total_vmrs"
    ),
    with = FALSE
  ]
}


category_summary <- rbindlist(
  list(
    
    summarize_categories(
      dt_wg,
      "Whole_genome"
    ),
    
    summarize_categories(
      dt_auto,
      "Autosomes"
    ),
    
    summarize_categories(
      dt_sex,
      "Sex_chromosome"
    )
    
  )
)

fwrite(
  category_summary,
  out_category_summary,
  sep = "\t",
  quote = FALSE
)

message("")

message(
  "Saved category summary: ",
  out_category_summary
)


# Print category summary -------------------------------------------------------

message("")

message(
  "VMR shuffled-analogue biological-category summary:"
)

print(
  category_summary
)


# Prepare matrices for ComplexUpset --------------------------------------------
#
# ComplexUpset expects binary/logical columns defining set membership.
# Convert the original 0/1 columns to logical values.

prepare_upset_matrix <- function(
    data_i
) {
  
  mat_i <- as.data.frame(
    data_i[
      ,
      ..contrast_cols
    ]
  )
  
  for (
    cc in contrast_cols
  ) {
    
    mat_i[[cc]] <- as.logical(
      mat_i[[cc]]
    )
  }
  
  rownames(
    mat_i
  ) <- data_i$cons_id
  
  mat_i
}


mat_wg <- prepare_upset_matrix(
  dt_wg
)

mat_auto <- prepare_upset_matrix(
  dt_auto
)

mat_sex <- prepare_upset_matrix(
  dt_sex
)

message("")

message(
  "mat_wg dims: ",
  paste(
    dim(mat_wg),
    collapse = " x "
  )
)

message(
  "mat_auto dims: ",
  paste(
    dim(mat_auto),
    collapse = " x "
  )
)

message(
  "mat_sex dims: ",
  paste(
    dim(mat_sex),
    collapse = " x "
  )
)


# UpSet plotting function ------------------------------------------------------
#
# Set-size bars are not displayed.
# Instead, the total number of VMRs belonging to each shuffled contrast is
# appended directly to the corresponding contrast label in the intersection
# matrix.
#
# Set sizes are calculated from all VMRs in the genomic partition, including
# VMRs belonging to contrast combinations not displayed among the 16 selected
# shuffled analogues of the biological intersections.
#
# Key options:
#
#   intersect = contrast_cols
#       Keeps the six shuffled set rows in the desired contrast order.
#
#   intersections = intersection_order
#       Restricts the plot to the same 16 structural combinations used in
#       the real-data analysis.
#
#   sort_intersections = FALSE
#       Prevents sorting by frequency and preserves the predefined order.
#
#   sort_sets = FALSE
#       Prevents shuffled contrast rows from being sorted by total size.
#
#   mode = "exclusive_intersection"
#       Counts a VMR only when its exact six-set membership pattern matches
#       the intersection being displayed.
#
#   keep_empty_groups = TRUE
#       Retains requested intersections even when their shuffled count is zero.

plot_upset <- function(
    mat,
    pdf_file,
    partition_label
) {
  
  if (
    nrow(mat) == 0
  ) {
    
    message(
      "Skipping ",
      partition_label,
      " plot: no rows"
    )
    
    return(
      invisible(NULL)
    )
  }
  
  message("")
  
  message(
    "Opening PDF: ",
    pdf_file
  )
  
  
  # Calculate total VMR set size for each shuffled contrast --------------------
  
  set_sizes_n <- colSums(
    mat[
      ,
      contrast_cols,
      drop = FALSE
    ]
  )
  
  set_labels <- setNames(
    paste0(
      contrast_cols,
      " (n = ",
      set_sizes_n[contrast_cols],
      ")"
    ),
    contrast_cols
  )
  
  
  # Plot -----------------------------------------------------------------------
  
  p <- ComplexUpset::upset(
    
    mat,
    
    # Sets
    intersect = contrast_cols,
    encode_sets = FALSE,
    sort_sets = FALSE,
    
    # Intersections
    intersections = intersection_order,
    sort_intersections = FALSE,
    mode = "exclusive_intersection",
    
    # Keep manually requested groups
    keep_empty_groups = TRUE,
    warn_when_dropping_groups = FALSE,
    
    # Remove set-size bar plot
    set_sizes = FALSE,
    
    # Intersection matrix with set-size counts in labels
    matrix = (
      ComplexUpset::intersection_matrix() +
        scale_y_discrete(
          labels = set_labels
        )
    ),
    
    # Main intersection-size bars
    base_annotations = list(
      
      "Intersection size" = (
        
        ComplexUpset::intersection_size(
          counts = TRUE,
          mapping = aes(
            fill = I("black")
          )
        ) +
          
          ylab(
            paste0(
              "Shared VMR\n",
              partition_label,
              " (n = ",
              nrow(mat),
              ")"
            )
          ) +
          
          theme(
            axis.title.y = element_text(
              size = 12
            ),
            axis.text.y = element_text(
              size = 10
            )
          )
        
      )
      
    ),
    
    name = "Contrast-sharing pattern"
    
  )
  
  pdf(
    pdf_file,
    width = 6,
    height = 4,
    useDingbats = FALSE
  )
  
  print(
    p
  )
  
  dev.off()
  
  message(
    "Saved PDF: ",
    pdf_file
  )
}


# Run plots --------------------------------------------------------------------

# Whole genome plot
plot_upset(
  mat = mat_wg,
  pdf_file = out_pdf_wg,
  partition_label = "Whole genome"
)

# Autosome plot
plot_upset(
  mat = mat_auto,
  pdf_file = out_pdf_auto,
  partition_label = "Autosomes"
)

# Sex chromosome plot
plot_upset(
  mat = mat_sex,
  pdf_file = out_pdf_sex,
  partition_label = "Sex chromosome"
)


# Final summary ----------------------------------------------------------------

n_selected_wg <- dt_wg[
  biological_category != "Other contrast combinations",
  .N
]

n_other_wg <- dt_wg[
  biological_category == "Other contrast combinations",
  .N
]

cat("\n")
cat("======================================\n")
cat("Shuffled VMR contrast-sharing UpSet\n")
cat("======================================\n")

cat(
  "Input file: ",
  normalizePath(
    input_file,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)

cat("\nShuffled contrast order:\n")

cat(
  "  ",
  paste(
    contrast_cols,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat("\nReal-to-shuffled contrast mapping:\n")
cat("  f_vs_i  -> sA_vs_sA_01\n")
cat("  f_vs_p  -> sA_vs_sA_02\n")
cat("  f_vs_y  -> sA_vs_sA_03\n")
cat("  i_vs_y  -> sM_vs_sM_01\n")
cat("  p_vs_i  -> sM_vs_sM_02\n")
cat("  y_vs_p  -> sM_vs_sM_03\n")

cat("\nBiological category order:\n")
cat("  1. Sex-associated\n")
cat("  2. Sex + morph\n")
cat("  3. Cross-phenotype\n")
cat("  4. Morph-specific\n")
cat("  5. Pairwise-specific\n")

cat("\n")
cat(
  "NOTE: category names denote shuffled analogues of the real-data\n"
)
cat(
  "      contrast-sharing patterns, not biological classification of\n"
)
cat(
  "      the shuffled VMRs.\n"
)

cat("\nIntersection columns:\n")
cat("  Column 1:      Sex-associated\n")
cat("  Columns 2-4:   Sex + morph\n")
cat("  Columns 5-7:   Cross-phenotype\n")
cat("  Columns 8-10:  Morph-specific\n")
cat("  Columns 11-16: Pairwise-specific\n")

cat(
  "\nRows in whole genome:    ",
  nrow(dt_wg),
  "\n",
  sep = ""
)

cat(
  "Rows in autosomes:       ",
  nrow(dt_auto),
  "\n",
  sep = ""
)

cat(
  "Rows in sex chromosome:  ",
  nrow(dt_sex),
  "\n",
  sep = ""
)

cat(
  "\nWhole-genome VMRs represented by the five categories: ",
  n_selected_wg,
  " / ",
  nrow(dt_wg),
  "\n",
  sep = ""
)

cat(
  "Whole-genome VMRs in other contrast combinations:      ",
  n_other_wg,
  " / ",
  nrow(dt_wg),
  "\n",
  sep = ""
)

cat("\nWrote:\n")

if (
  file.exists(
    out_pdf_wg
  )
) {
  
  cat(
    "  ",
    normalizePath(
      out_pdf_wg,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n",
    sep = ""
  )
}

if (
  file.exists(
    out_pdf_auto
  )
) {
  
  cat(
    "  ",
    normalizePath(
      out_pdf_auto,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n",
    sep = ""
  )
}

if (
  file.exists(
    out_pdf_sex
  )
) {
  
  cat(
    "  ",
    normalizePath(
      out_pdf_sex,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n",
    sep = ""
  )
}

if (
  file.exists(
    out_category_summary
  )
) {
  
  cat(
    "  ",
    normalizePath(
      out_category_summary,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n",
    sep = ""
  )
}


# Save contrast set sizes -------------------------------------------------------

get_set_sizes <- function(
    mat,
    partition_i
) {
  
  data.table(
    partition = partition_i,
    contrast = contrast_cols,
    set_size = as.integer(
      colSums(
        mat[
          ,
          contrast_cols,
          drop = FALSE
        ]
      )
    )
  )
}


set_size_summary <- rbindlist(
  list(
    
    get_set_sizes(
      mat_wg,
      "Whole_genome"
    ),
    
    get_set_sizes(
      mat_auto,
      "Autosomes"
    ),
    
    get_set_sizes(
      mat_sex,
      "Sex_chromosome"
    )
    
  )
)

fwrite(
  set_size_summary,
  out_set_size_summary,
  sep = "\t",
  quote = FALSE
)

message("")

message(
  "Saved set-size summary: ",
  out_set_size_summary
)

if (
  file.exists(
    out_set_size_summary
  )
) {
  
  cat(
    "  ",
    normalizePath(
      out_set_size_summary,
      winslash = "/",
      mustWork = FALSE
    ),
    "\n",
    sep = ""
  )
}

cat("======================================\n")