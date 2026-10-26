#!/usr/bin/env Rscript

# DATE: 20260824
# AUTHOR: MZF
# SCRIPT: 03_complexUpset_02_plot_V02.R
# VERSION: 02
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Build UpSet plots showing how VMRs are shared across phenotype contrasts.
#
#   Intersections are restricted to biologically interpretable contrast-sharing
#   patterns and displayed in the following order:
#
#     1) Sex-associated
#     2) Sex + morph
#     3) Cross-phenotype
#     4) Morph-specific
#     5) Pairwise-specific
#
#   VMRs belonging to other exact contrast combinations are retained and
#   classified as:
#
#     Other contrast combinations
#
#   In addition to the UpSet plots and category summaries, this version writes
#   a PER-VMR annotation table containing the biological category assigned to
#   every VMR. This table can be joined downstream to methylation tables using
#   cons_id.
#
# INPUT:
#   vmr_real.contrastSharing.binary.tsv
#
# OUTPUTS:
#   vmr_real.contrastSharing.WG.upset.pdf
#   vmr_real.contrastSharing.AUTO.upset.pdf
#   vmr_real.contrastSharing.SEXCH.upset.pdf
#
#   vmr_real.contrastSharing.biologicalCategories.summary.tsv
#   vmr_real.contrastSharing.biologicalCategories.perVMR.tsv
#   vmr_real.contrastSharing.setSizes.summary.txt
#
# PER-VMR OUTPUT:
#   chr
#   start
#   end
#   cons_id
#   dmr_ids
#   chromosome_type
#   n_contrasts
#   contrast_list
#   computed_contrast_combo
#   biological_category
#   biological_category_order
#   biological_intersection_order
#   selected_biological_intersection
#   f_vs_i
#   f_vs_p
#   f_vs_y
#   i_vs_y
#   p_vs_i
#   y_vs_p
#
# USAGE:
#   Rscript 03_complexUpset_02_plot_V02.R
#
# or:
#
#   Rscript 03_complexUpset_02_plot_V02.R input.tsv output_prefix
#
# REQUIRED R PACKAGES:
#   data.table
#   ggplot2
#   ComplexUpset
#
# V02:
#   - writes per-VMR biological-category annotations
#   - records exact biological-intersection order
#   - records whether a VMR belongs to one of the 16 selected intersections
#   - fixes output-prefix initialization for set-size summary
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%


# Libraries --------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexUpset)
})


# Command-line arguments -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)


# Defaults ---------------------------------------------------------------------

input_file <- "vmr_real.contrastSharing.binary.tsv"

outdir <- "../plots/complexupset/"
output_pref <- file.path(outdir, "vmr_real.contrastSharing")

if (length(args) >= 1) input_file <- args[1]
if (length(args) >= 2) output_pref <- args[2]

dir.create(dirname(output_pref), showWarnings = FALSE, recursive = TRUE)


# Output files -----------------------------------------------------------------

out_pdf_wg <- paste0(output_pref, ".WG.upset.pdf")
out_pdf_auto <- paste0(output_pref, ".AUTO.upset.pdf")
out_pdf_sex <- paste0(output_pref, ".SEXCH.upset.pdf")

out_category_summary <- paste0(output_pref,".biologicalCategories.summary.tsv")
out_vmr_annotation <- paste0(output_pref,".biologicalCategories.perVMR.tsv")
out_set_size_summary <- paste0(output_pref,".setSizes.summary.txt")


# Read input -------------------------------------------------------------------

dt <- fread(input_file)

message("Loaded input file: ", input_file)
message("Rows: ", nrow(dt), " | Cols: ", ncol(dt))


# Check required metadata columns ---------------------------------------------

required_cols <- c(
  "chr",
  "start",
  "end",
  "cons_id",
  "chromosome_type"
)

missing_cols <- setdiff(required_cols, names(dt))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (anyDuplicated(dt$cons_id) > 0) {
  stop("cons_id is not unique in the input table.")
}


# Define contrast columns ------------------------------------------------------

contrast_cols <- c(
  "f_vs_i",
  "f_vs_p",
  "f_vs_y",
  "i_vs_y",
  "p_vs_i",
  "y_vs_p"
)

missing_contrasts <- setdiff(contrast_cols, names(dt))

if (length(missing_contrasts) > 0) {
  stop(
    "Missing required contrast columns: ",
    paste(missing_contrasts, collapse = ", ")
  )
}


# Validate binary contrast columns --------------------------------------------

for (cc in contrast_cols) {
  
  values_i <- unique(dt[[cc]])
  values_i <- values_i[!is.na(values_i)]
  
  if (!all(values_i %in% c(0, 1))) {
    stop("Contrast column ", cc, " contains values other than 0/1.")
  }
  
  if (anyNA(dt[[cc]])) {
    stop("Contrast column ", cc, " contains NA values.")
  }
}

message("Contrast order: ", paste(contrast_cols, collapse = ", "))


# Biological intersection categories ------------------------------------------

# 1. Sex-associated ------------------------------------------------------------

sex_associated <- list(
  c(
    "f_vs_i",
    "f_vs_p",
    "f_vs_y"
  )
)


# 2. Sex + morph ---------------------------------------------------------------

sex_plus_morph <- list(
  
  c(
    "f_vs_i",
    "f_vs_p",
    "f_vs_y",
    "i_vs_y",
    "p_vs_i"
  ),
  
  c(
    "f_vs_i",
    "f_vs_p",
    "f_vs_y",
    "i_vs_y",
    "y_vs_p"
  ),
  
  c(
    "f_vs_i",
    "f_vs_p",
    "f_vs_y",
    "p_vs_i",
    "y_vs_p"
  )
)


# 3. Cross-phenotype -----------------------------------------------------------

cross_phenotype <- list(
  
  # Female + Immaculata versus Parae + Yellow
  c(
    "f_vs_p",
    "f_vs_y",
    "i_vs_y",
    "p_vs_i"
  ),
  
  # Female + Parae versus Immaculata + Yellow
  c(
    "f_vs_i",
    "f_vs_y",
    "p_vs_i",
    "y_vs_p"
  ),
  
  # Female + Yellow versus Immaculata + Parae
  c(
    "f_vs_i",
    "f_vs_p",
    "i_vs_y",
    "y_vs_p"
  )
)


# 4. Morph-specific ------------------------------------------------------------

morph_specific <- list(
  
  # Immaculata-specific
  c(
    "f_vs_i",
    "i_vs_y",
    "p_vs_i"
  ),
  
  # Parae-specific
  c(
    "f_vs_p",
    "p_vs_i",
    "y_vs_p"
  ),
  
  # Yellow-specific
  c(
    "f_vs_y",
    "i_vs_y",
    "y_vs_p"
  )
)


# 5. Pairwise-specific ---------------------------------------------------------

pairwise_specific <- list(
  c("f_vs_i"),
  c("f_vs_p"),
  c("f_vs_y"),
  c("i_vs_y"),
  c("p_vs_i"),
  c("y_vs_p")
)


# Combine intersections in biological order ----------------------------------

intersection_order <- c(
  sex_associated,
  sex_plus_morph,
  cross_phenotype,
  morph_specific,
  pairwise_specific
)

if (length(intersection_order) != 16) {
  stop(
    "Expected 16 biological intersections, but found ",
    length(intersection_order),
    "."
  )
}

message("Number of biological intersections to plot: ", length(intersection_order))


# Build intersection lookup table ---------------------------------------------

intersection_key <- data.table(
  
  biological_intersection_order = seq_along(intersection_order),
  
  biological_category = c(
    rep("Sex-associated", length(sex_associated)),
    rep("Sex + morph", length(sex_plus_morph)),
    rep("Cross-phenotype", length(cross_phenotype)),
    rep("Morph-specific", length(morph_specific)),
    rep("Pairwise-specific", length(pairwise_specific))
  ),
  
  contrast_combo = vapply(
    intersection_order,
    function(x) paste(x, collapse = "&"),
    character(1)
  )
)

message("")
message("Biological intersection order:")

for (i in seq_len(nrow(intersection_key))) {
  
  message(
    sprintf(
      "%2d  %-20s  %s",
      intersection_key$biological_intersection_order[i],
      intersection_key$biological_category[i],
      intersection_key$contrast_combo[i]
    )
  )
}


# Compute contrast combination directly from binary columns -------------------

dt[
  ,
  computed_contrast_combo := apply(
    .SD,
    1,
    function(x) {
      
      present <- contrast_cols[as.integer(x) == 1L]
      
      if (length(present) == 0) return("none")
      
      paste(present, collapse = "&")
    }
  ),
  .SDcols = contrast_cols
]


# Assign biological category ---------------------------------------------------

category_levels <- c(
  "Sex-associated",
  "Sex + morph",
  "Cross-phenotype",
  "Morph-specific",
  "Pairwise-specific",
  "Other contrast combinations"
)

combo_to_category <- setNames(
  intersection_key$biological_category,
  intersection_key$contrast_combo
)

combo_to_intersection_order <- setNames(
  intersection_key$biological_intersection_order,
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
  is.na(biological_category),
  biological_category := "Other contrast combinations"
]

dt[
  ,
  biological_category_order :=
    match(
      biological_category,
      category_levels
    )
]

dt[
  ,
  biological_intersection_order :=
    as.integer(
      unname(
        combo_to_intersection_order[
          computed_contrast_combo
        ]
      )
    )
]

dt[
  ,
  selected_biological_intersection :=
    biological_category != "Other contrast combinations"
]


# Write per-VMR biological annotation table -----------------------------------

annotation_cols <- c(
  "chr",
  "start",
  "end",
  "cons_id",
  "dmr_ids",
  "chromosome_type",
  "n_contrasts",
  "contrast_list",
  "computed_contrast_combo",
  "biological_category",
  "biological_category_order",
  "biological_intersection_order",
  "selected_biological_intersection",
  contrast_cols
)

annotation_cols <- intersect(
  annotation_cols,
  names(dt)
)

vmr_annotation <- dt[
  ,
  ..annotation_cols
]

setorder(
  vmr_annotation,
  chr,
  start,
  end
)

if (anyDuplicated(vmr_annotation$cons_id) > 0) {
  stop("Duplicated cons_id detected in per-VMR annotation table.")
}

fwrite(
  vmr_annotation,
  out_vmr_annotation,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

message("")
message("Saved per-VMR annotation table: ", out_vmr_annotation)
message("Annotated VMRs: ", nrow(vmr_annotation))


# Split datasets by genomic partition -----------------------------------------

dt_wg <- copy(dt)

dt_auto <- dt[
  chromosome_type %in% c(
    "Autosome",
    "AUTO",
    "autosome"
  )
]

if (nrow(dt_auto) == 0) {
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

if (nrow(dt_sex) == 0) {
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
message("dt_wg rows: ", nrow(dt_wg))
message("dt_auto rows: ", nrow(dt_auto))
message("dt_sex rows: ", nrow(dt_sex))


# Summarize biological categories ---------------------------------------------

summarize_categories <- function(data_i, partition_i) {
  
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
    total_vmrs := nrow(data_i)
  ]
  
  summary_i[
    ,
    percentage :=
      if (nrow(data_i) > 0) {
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
    summarize_categories(dt_wg, "Whole_genome"),
    summarize_categories(dt_auto, "Autosomes"),
    summarize_categories(dt_sex, "Sex_chromosome")
  )
)

fwrite(
  category_summary,
  out_category_summary,
  sep = "\t",
  quote = FALSE
)

message("")
message("Saved category summary: ", out_category_summary)

message("")
message("VMR biological-category summary:")
print(category_summary)


# Prepare matrices for ComplexUpset -------------------------------------------

prepare_upset_matrix <- function(data_i) {
  
  mat_i <- as.data.frame(
    data_i[
      ,
      ..contrast_cols
    ]
  )
  
  for (cc in contrast_cols) {
    mat_i[[cc]] <- as.logical(mat_i[[cc]])
  }
  
  rownames(mat_i) <- data_i$cons_id
  
  mat_i
}

mat_wg <- prepare_upset_matrix(dt_wg)
mat_auto <- prepare_upset_matrix(dt_auto)
mat_sex <- prepare_upset_matrix(dt_sex)

message("")
message("mat_wg dims: ", paste(dim(mat_wg), collapse = " x "))
message("mat_auto dims: ", paste(dim(mat_auto), collapse = " x "))
message("mat_sex dims: ", paste(dim(mat_sex), collapse = " x "))


# UpSet plotting function ------------------------------------------------------

plot_upset <- function(mat, pdf_file, partition_label) {
  
  if (nrow(mat) == 0) {
    message("Skipping ", partition_label, " plot: no rows")
    return(invisible(NULL))
  }
  
  message("")
  message("Opening PDF: ", pdf_file)
  
  # Total number of VMRs in each contrast
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
    # Intersection matrix
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
            axis.title.y = element_text(size = 12),
            axis.text.y = element_text(size = 10)
          )
      )
    ),
    
    name = "Contrast-sharing pattern"
  )
  
  pdf( pdf_file,width = 6, height = 4,useDingbats = FALSE)
  print(p)
  dev.off()
  
  message("Saved PDF: ", pdf_file)
}


# Run plots --------------------------------------------------------------------

plot_upset(
  mat = mat_wg,
  pdf_file = out_pdf_wg,
  partition_label = "Whole genome"
)

plot_upset(
  mat = mat_auto,
  pdf_file = out_pdf_auto,
  partition_label = "Autosomes"
)

plot_upset(
  mat = mat_sex,
  pdf_file = out_pdf_sex,
  partition_label = "Sex chromosome"
)


# Save contrast set sizes ------------------------------------------------------

get_set_sizes <- function(mat, partition_i) {
  
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
    get_set_sizes(mat_wg, "Whole_genome"),
    get_set_sizes(mat_auto, "Autosomes"),
    get_set_sizes(mat_sex, "Sex_chromosome")
  )
)

fwrite(
  set_size_summary,
  out_set_size_summary,
  sep = "\t",
  quote = FALSE
)

message("")
message("Saved set-size summary: ", out_set_size_summary)


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
cat("VMR contrast-sharing UpSet\n")
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

cat("\nContrast order:\n")
cat("  ", paste(contrast_cols, collapse = ", "), "\n", sep = "")

cat("\nBiological category order:\n")
cat("  1. Sex-associated\n")
cat("  2. Sex + morph\n")
cat("  3. Cross-phenotype\n")
cat("  4. Morph-specific\n")
cat("  5. Pairwise-specific\n")
cat("  6. Other contrast combinations\n")

cat("\nIntersection columns:\n")
cat("  Column 1:      Sex-associated\n")
cat("  Columns 2-4:   Sex + morph\n")
cat("  Columns 5-7:   Cross-phenotype\n")
cat("  Columns 8-10:  Morph-specific\n")
cat("  Columns 11-16: Pairwise-specific\n")

cat("\nRows in whole genome:    ", nrow(dt_wg), "\n", sep = "")
cat("Rows in autosomes:       ", nrow(dt_auto), "\n", sep = "")
cat("Rows in sex chromosome:  ", nrow(dt_sex), "\n", sep = "")

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

cat("\nPer-VMR annotations written: ", nrow(vmr_annotation), "\n", sep = "")

cat("\nWrote:\n")

output_files <- c(
  out_pdf_wg,
  out_pdf_auto,
  out_pdf_sex,
  out_category_summary,
  out_vmr_annotation,
  out_set_size_summary
)

for (ff in output_files) {
  
  if (file.exists(ff)) {
    cat(
      "  ",
      normalizePath(
        ff,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n",
      sep = ""
    )
  }
}

cat("======================================\n")