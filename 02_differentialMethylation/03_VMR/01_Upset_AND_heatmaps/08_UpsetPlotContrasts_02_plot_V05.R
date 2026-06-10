#!/usr/bin/env Rscript
# DATE: 20260401
# AUTHOR: MZF & ChatGPT
# SCRIPT: 03_plot_consensusDMR_contrast_upset_PDFonly.R
# VERSION: 05
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Build UpSet plots showing how consensus DMRs are shared across contrasts.
#   Save one PDF per genomic partition:
#     1) Whole genome
#     2) Autosomes
#     3) Sex chromosome
#
# INPUT:
#   dmr_real_consensus.contrastSharing.binary.tsv
#
# OUTPUT:
#   dmr_real_consensus.contrastSharing.WG.upset.pdf
#   dmr_real_consensus.contrastSharing.AUTO.upset.pdf
#   dmr_real_consensus.contrastSharing.SEXCH.upset.pdf
#
# USAGE:
#   Rscript 03_plot_consensusDMR_contrast_upset_PDFonly.R
#
#   or
#
#   Rscript 03_plot_consensusDMR_contrast_upset_PDFonly.R input.tsv output_prefix
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

# suppressPackageStartupMessages({
library(data.table)
library(UpSetR)
library(grid)
# })

args <- commandArgs(trailingOnly = TRUE)

# -----------------------------
# Defaults
# -----------------------------
input_file <- "dmr_real_consensus.contrastSharing.binary.tsv"
outdir      <- "../plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
output_pref <- file.path(outdir, "dmr_real_consensus.contrastSharing")

if (length(args) >= 1) input_file <- args[1]
if (length(args) >= 2) output_pref <- args[2]

out_pdf_wg <- paste0(output_pref, ".WG.upset.pdf")
out_pdf_auto <- paste0(output_pref, ".AUTO.upset.pdf")
out_pdf_sex <- paste0(output_pref, ".SEXCH.upset.pdf")

# -----------------------------
# Read input
# -----------------------------
dt <- fread(input_file)

message("Loaded input file: ", input_file)
message("Rows: ", nrow(dt), " | Cols: ", ncol(dt))

required_cols <- c("cons_id", "chromosome_type", "chr")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# Detect contrast columns
# -----------------------------
meta_cols <- c(
  "chr", "start", "end", "cons_id", "dmr_ids", "chromosome_type",
  "n_contrasts", "contrast_list", "contrast_combo"
)

contrast_cols <- names(dt)[
  sapply(dt, function(x) {
    ux <- unique(x)
    ux <- ux[!is.na(ux)]
    length(ux) > 0 && all(ux %in% c(0, 1))
  }) &
    !(names(dt) %in% meta_cols)
]

if (length(contrast_cols) < 2) {
  stop("Could not detect at least two binary contrast columns.")
}

preferred_order <- c("f_vs_i", "f_vs_p", "f_vs_y", "i_vs_y", "p_vs_i", "y_vs_p")
contrast_cols <- c(
  intersect(preferred_order, contrast_cols),
  setdiff(sort(contrast_cols), preferred_order)
)

message("Detected contrast columns: ", paste(contrast_cols, collapse = ", "))

# -----------------------------
# Split datasets
# -----------------------------
dt_wg <- copy(dt)

dt_auto <- dt[chromosome_type %in% c("Autosome", "AUTO", "autosome")]
if (nrow(dt_auto) == 0) {
  dt_auto <- dt[grepl("auto", chromosome_type, ignore.case = TRUE)]
}

dt_sex <- dt[
  chromosome_type %in% c("Sex_Ch", "SEXCH", "sex_ch", "sex", "SexChr") |
    chr == "Parae_12"
]
if (nrow(dt_sex) == 0) {
  dt_sex <- dt[grepl("sex", chromosome_type, ignore.case = TRUE) | chr == "Parae_12"]
}

message("dt_wg rows: ", nrow(dt_wg))
message("dt_auto rows: ", nrow(dt_auto))
message("dt_sex rows: ", nrow(dt_sex))

# -----------------------------
# Prepare matrices
# -----------------------------
mat_wg <- as.data.frame(dt_wg[, ..contrast_cols])
for (cc in names(mat_wg)) mat_wg[[cc]] <- as.integer(mat_wg[[cc]])
rownames(mat_wg) <- dt_wg$cons_id

mat_auto <- as.data.frame(dt_auto[, ..contrast_cols])
for (cc in names(mat_auto)) mat_auto[[cc]] <- as.integer(mat_auto[[cc]])
rownames(mat_auto) <- dt_auto$cons_id

mat_sex <- as.data.frame(dt_sex[, ..contrast_cols])
for (cc in names(mat_sex)) mat_sex[[cc]] <- as.integer(mat_sex[[cc]])
rownames(mat_sex) <- dt_sex$cons_id

message("mat_wg dims: ", paste(dim(mat_wg), collapse = " x "))
message("mat_auto dims: ", paste(dim(mat_auto), collapse = " x "))
message("mat_sex dims: ", paste(dim(mat_sex), collapse = " x "))

# -----------------------------
# Whole genome plot
# -----------------------------
if (nrow(mat_wg) > 0) {
  message("Opening PDF: ", out_pdf_wg)
  pdf(out_pdf_wg, width = 10, height = 7)
  grid.newpage()

  print(UpSetR::upset(
    mat_wg,
    sets = contrast_cols,
    nsets = length(contrast_cols),
    # intersections # Specific intersections to include in plot entered as a list of lists. Ex: list(list("Set name1", "Set name2"), list("Set name1", "Set name3")). If data is entered into this parameter the only data shown on the UpSet plot will be the specific intersections listed.
    nintersects = NA,
    keep.order = TRUE,
    order.by = "freq",
    mainbar.y.label = paste0(
      "Shared consensus DMRs\nWhole genome (n = ", nrow(dt_wg), ")"
    ),
    sets.x.label = "Consensus DMRs per contrast",
    mb.ratio = c(0.65, 0.35),
    text.scale = c(1.4, 1.2, 1.2, 1.0, 1.2, 0.7),
    main.bar.color = "black",
    sets.bar.color = "grey40",
    matrix.color = "black",
    point.size = 3,
    line.size = 1
  ))

  dev.off()
  message("Saved PDF: ", out_pdf_wg)
} else {
  message("Skipping WG plot: no rows")
}

# -----------------------------
# Autosome plot
# -----------------------------
if (nrow(mat_auto) > 0) {
  message("Opening PDF: ", out_pdf_auto)
  pdf(out_pdf_auto, width = 10, height = 7)
  grid.newpage()

  print(UpSetR::upset(
    mat_auto,
    sets = contrast_cols,
    nsets = length(contrast_cols),
    nintersects = NA,
    keep.order = TRUE,
    order.by = "freq",
    mainbar.y.label = paste0(
      "Shared consensus DMRs\nAutosomes (n = ", nrow(dt_auto), ")"
    ),
    sets.x.label = "Consensus DMRs per contrast",
    mb.ratio = c(0.65, 0.35),
    text.scale = c(1.4, 1.2, 1.2, 1.0, 1.2, 0.7),
    main.bar.color = "black",
    sets.bar.color = "grey40",
    matrix.color = "black",
    point.size = 3,
    line.size = 1
  ))

  dev.off()
  message("Saved PDF: ", out_pdf_auto)
} else {
  message("Skipping AUTO plot: no rows")
}

# -----------------------------
# Sex chromosome plot
# -----------------------------
if (nrow(mat_sex) > 0) {
  message("Opening PDF: ", out_pdf_sex)
  pdf(out_pdf_sex, width = 10, height = 7)
  grid.newpage()

  print(UpSetR::upset(
    mat_sex,
    sets = contrast_cols,
    nsets = length(contrast_cols),
    nintersects = NA,
    keep.order = TRUE,
    order.by = "freq",
    mainbar.y.label = paste0(
      "Shared consensus DMRs\nSex chromosome (n = ", nrow(dt_sex), ")"
    ),
    sets.x.label = "Consensus DMRs per contrast",
    mb.ratio = c(0.65, 0.35),
    text.scale = c(1.4, 1.2, 1.2, 1.0, 1.2, 1.2),
    main.bar.color = "black",
    sets.bar.color = "grey40",
    matrix.color = "black",
    point.size = 3,
    line.size = 1
  ))

  dev.off()
  message("Saved PDF: ", out_pdf_sex)
} else {
  message("Skipping SEX plot: no rows")
}

# -----------------------------
# Final summary
# -----------------------------
cat("======================================\n")
cat("Consensus DMR contrast-sharing UpSet\n")
cat("======================================\n")
cat("Input file: ", normalizePath(input_file, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("Detected contrast columns:\n")
cat("  ", paste(contrast_cols, collapse = ", "), "\n", sep = "")
cat("\n")
cat("Rows in whole genome:    ", nrow(dt_wg), "\n", sep = "")
cat("Rows in autosomes:       ", nrow(dt_auto), "\n", sep = "")
cat("Rows in sex chromosome:  ", nrow(dt_sex), "\n", sep = "")
cat("\n")
cat("Wrote:\n")
if (file.exists(out_pdf_wg)) cat("  ", normalizePath(out_pdf_wg, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (file.exists(out_pdf_auto)) cat("  ", normalizePath(out_pdf_auto, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (file.exists(out_pdf_sex)) cat("  ", normalizePath(out_pdf_sex, winslash = "/", mustWork = FALSE), "\n", sep = "")
