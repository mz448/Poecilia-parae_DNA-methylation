#!/usr/bin/env Rscript
# DATE: 20260221
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   (1) Load the mapped consensus DMR table (DMR presence via NA/non-NA).
#   (2) Filter to the subset of consensus DMRs that passed the methylation-diff
#       filter used for heatmaps (from 05_heatmap_V07 outputs).
#   (3) Plot overlap using:
#       - VennDiagram (classic; not area-scaled)
#       - eulerr (area-scaled Euler)
#       - UpSet plots (area-independent; robust for multiple sets)
#
# Notes:
#   - Presence rule (DMR-sharing): 1 if meth_* is NOT NA (morph participates
#     in at least one contributing DMR in that consensus region), else 0.
#   - Generates separate outputs for AUTOSOME vs CH12/Sex_Ch.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
  library(VennDiagram)
  library(grid)
  library(eulerr)
  library(UpSetR)
})

# ------------------------
# Inputs / settings
# ------------------------
PLOTS <- "../plots_threshold_0.0"
OUTDIR <- paste0(PLOTS,"/venn_methylation_filtered")
INFILE_MAPPED <- "dmr_real_consensus.morph_mapped.tsv"

# This is produced by 05_heatmap_V07.R (filtered observed table)
# (If you changed outdir/threshold, point this to the correct file.)
INFILE_FILTERED_CONS <- paste0(PLOTS,"/dmr_real_consensus.morph_observed.filtered.diffGE_0.tsv")

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)


# ------------------------
# Safe PDF helper (prevents empty/corrupt pdfs if errors occur)
# ------------------------
safe_pdf <- function(path, width=7.5, height=7.5, expr) {
  pdf(path, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  tryCatch(
    expr,
    error = function(e) {
      message("ERROR while plotting: ", basename(path))
      message(conditionMessage(e))
      stop(e)
    }
  )
}
# ------------------------
# Load mapped consensus table
# ------------------------
dt <- fread(INFILE_MAPPED)

meth_cols <- c("meth_female","meth_immaculata","meth_parae","meth_yellow")
stopifnot(all(meth_cols %in% names(dt)))

# region_id MUST be consensus ID when available
if (!"cons_id" %in% names(dt)) stop("Input mapped file must contain cons_id: ", INFILE_MAPPED)
dt[, region_id := as.character(cons_id)]

# QC: one row per consensus ID
dup <- dt[, .N, by = region_id][N > 1]
if (nrow(dup) > 0) {
  print(head(dup, 20))
  stop("Duplicate cons_id (region_id) in mapped table; expected 1 row per consensus region.")
}

# chromosome_type (derive if missing)
if (!"chromosome_type" %in% names(dt)) {
  dt[, chromosome_type := ifelse(chr == "Parae_12", "Sex_Ch", "Autosome")]
}

# ------------------------
# Load filtered consensus IDs from 05_heatmap_V07
# ------------------------
filt <- fread(INFILE_FILTERED_CONS)
if (!"cons_id" %in% names(filt)) stop("Filtered file must contain cons_id: ", INFILE_FILTERED_CONS)
keep_ids <- unique(as.character(filt$cons_id))

# Filter mapped table to filtered consensus set
dt_f <- dt[region_id %in% keep_ids]

# QC
qc_counts <- rbind(
  dt[,   .(n=.N), by=chromosome_type][, stage := "before_filter"],
  dt_f[, .(n=.N), by=chromosome_type][, stage := "after_filter"]
)
setcolorder(qc_counts, c("stage","chromosome_type","n"))
fwrite(qc_counts, file.path(OUTDIR, "QC_counts_before_after_filter.tsv"), sep = "\t")

# ------------------------
# Binary presence matrix (DMR-sharing)
# ------------------------
dt_f[, female     := as.integer(!is.na(meth_female))]
dt_f[, immaculata := as.integer(!is.na(meth_immaculata))]
dt_f[, parae      := as.integer(!is.na(meth_parae))]
dt_f[, yellow     := as.integer(!is.na(meth_yellow))]

bin <- dt_f[, .(region_id, chr, start, end, chromosome_type, female, immaculata, parae, yellow)]
fwrite(bin, file.path(OUTDIR, "binary_DMRsharing_matrix.filtered.tsv"), sep = "\t")

# ------------------------
# Helpers
# ------------------------
save_euler_pdf <- function(fit, filename, fill_cols, main = NULL) {
  pdf(file.path(OUTDIR, filename), width = 7.5, height = 7.5)
  p <- plot(
    fit,
    quantities = TRUE,
    labels = TRUE,
    fills = list(fill = fill_cols, alpha = 0.55),
    edges = list(col = "grey30")
  )
  grid::grid.newpage()
  grid::grid.draw(p)
  if (!is.null(main)) grid::grid.text(main, y = grid::unit(0.97, "npc"))
  dev.off()
}

make_sets <- function(bin_sub) {
  set_F <- unique(bin_sub[female == 1, region_id])
  set_I <- unique(bin_sub[immaculata == 1, region_id])
  set_P <- unique(bin_sub[parae == 1, region_id])
  set_Y <- unique(bin_sub[yellow == 1, region_id])
  set_MALES <- unique(c(set_I, set_P, set_Y))
  list(F=set_F, I=set_I, P=set_P, Y=set_Y, MALES=set_MALES)
}

plot_venn_euler <- function(sets, suffix) {
  
  # ---- VennDiagram A/B/C ----
  venn4 <- venn.diagram(
    x = list(Female = sets$F, Immaculata = sets$I, Parae = sets$P, Yellow = sets$Y),
    filename = NULL,
    fill = c("#149954","#888888","#9e9ac8","#fd8d3c"),
    alpha = 0.55, cex = 0.9, cat.cex = 0.9, margin = 0.10
  )
  pdf(file.path(OUTDIR, paste0("Venn_A_all4morphs_", suffix, ".pdf")), 7.5, 7.5)
  grid.newpage(); grid.draw(venn4); dev.off()
  
  venn2 <- venn.diagram(
    x = list(Female = sets$F, All_males_union = sets$MALES),
    filename = NULL,
    fill = c("#149954","#4D4D4D"),
    alpha = 0.55, cex = 1.0, cat.cex = 1.0, margin = 0.10
  )
  pdf(file.path(OUTDIR, paste0("Venn_B_female_vs_allMalesUnion_", suffix, ".pdf")), 7.5, 7.5)
  grid.newpage(); grid.draw(venn2); dev.off()
  
  venn3 <- venn.diagram(
    x = list(Immaculata = sets$I, Parae = sets$P, Yellow = sets$Y),
    filename = NULL,
    fill = c("#888888","#9e9ac8","#fd8d3c"),
    alpha = 0.55, cex = 0.95, cat.cex = 0.95, margin = 0.10
  )
  pdf(file.path(OUTDIR, paste0("Venn_C_malesOnly_", suffix, ".pdf")), 7.5, 7.5)
  grid.newpage(); grid.draw(venn3); dev.off()
  
  # ---- eulerr A/B/C (scaled) ----
  fit4 <- euler(list(Female = sets$F, Immaculata = sets$I, Parae = sets$P, Yellow = sets$Y))
  save_euler_pdf(fit4,
                 paste0("Euler_A_all4morphs_scaled_", suffix, ".pdf"),
                 c("#149954","#888888","#9e9ac8","#fd8d3c"),
                 main = paste0("Euler scaled — all 4 morphs (", suffix, ")"))
  
  fit2 <- euler(list(Female = sets$F, All_males_union = sets$MALES))
  save_euler_pdf(fit2,
                 paste0("Euler_B_female_vs_allMalesUnion_scaled_", suffix, ".pdf"),
                 c("#149954","#4D4D4D"),
                 main = paste0("Euler scaled — Female vs All males (", suffix, ")"))
  
  fit3 <- euler(list(Immaculata = sets$I, Parae = sets$P, Yellow = sets$Y))
  save_euler_pdf(fit3,
                 paste0("Euler_C_malesOnly_scaled_", suffix, ".pdf"),
                 c("#888888","#9e9ac8","#fd8d3c"),
                 main = paste0("Euler scaled — males only (", suffix, ")"))
}

# plot_upset <- function(bin_sub, suffix) {
#   # UpSetR expects data frame with 0/1 columns
#   x <- as.data.frame(bin_sub[, .(female, immaculata, parae, yellow)])
#   pdf(file.path(OUTDIR, paste0("UpSet_all4morphs_", suffix, ".pdf")), width = 7.5, height = 5.5)
#   upset(x,
#         sets = c("female","immaculata","parae","yellow"),
#         nsets = 4,
#         nintersects = 20,
#         order.by = "freq",
#         keep.order = TRUE)
#   dev.off()
#   
#   # Males only
#   x_m <- as.data.frame(bin_sub[, .(immaculata, parae, yellow)])
#   pdf(file.path(OUTDIR, paste0("UpSet_malesOnly_", suffix, ".pdf")), width = 7.0, height = 5.0)
#   upset(x_m,
#         sets = c("immaculata","parae","yellow"),
#         nsets = 3,
#         nintersects = 15,
#         order.by = "freq",
#         keep.order = TRUE)
#   dev.off()
# }


plot_upset <- function(bin_sub, suffix) {
  
  if (nrow(bin_sub) == 0) {
    message("Skipping UpSet for ", suffix, " (0 rows).")
    return(invisible(NULL))
  }
  
  # all 4
  x <- as.data.frame(bin_sub[, .(female, immaculata, parae, yellow)])
  safe_pdf(file.path(OUTDIR, paste0("UpSet_all4morphs_", suffix, ".pdf")), 7.5, 5.5, {
    # grid::grid.newpage()
    UpSetR::upset(x,
                  sets = c("female","immaculata","parae","yellow"),
                  nsets = 4,
                  # nintersects = 20,
                  order.by = "freq",
                  keep.order = TRUE
    )
  })

  
  # # males only
  # x_m <- as.data.frame(bin_sub[, .(immaculata, parae, yellow)])
  # safe_pdf(file.path(OUTDIR, paste0("UpSet_malesOnly_", suffix, ".pdf")), 7.0, 5.0, {
  #   grid::grid.newpage()
  #   UpSetR::upset(x_m,
  #                 sets = c("immaculata","parae","yellow"),
  #                 nsets = 3,
  #                 # nintersects = 15,
  #                 order.by = "freq",
  #                 keep.order = TRUE
  #   )
  # })
}

# ------------------------
# Split by chromosome_type and plot
# ------------------------
bin_auto <- bin[chromosome_type == "Autosome"]
bin_sex  <- bin[chromosome_type == "Sex_Ch"]

sets_auto <- make_sets(bin_auto)
sets_sex  <- make_sets(bin_sex)

plot_venn_euler(sets_auto, "AUTOSOME.filtered")
plot_venn_euler(sets_sex,  "CH12.filtered")

plot_upset(bin_auto, "AUTOSOME.filtered")
plot_upset(bin_sex,  "CH12.filtered")

# ------------------------
# Optional: write overlap counts
# ------------------------
counts_from_sets <- function(sets, panel) {
  data.table(
    panel = panel,
    n_F = length(sets$F),
    n_I = length(sets$I),
    n_P = length(sets$P),
    n_Y = length(sets$Y),
    n_all4 = length(Reduce(intersect, list(sets$F, sets$I, sets$P, sets$Y))),
    n_F_and_anyMale = length(intersect(sets$F, sets$MALES)),
    n_all3_males = length(Reduce(intersect, list(sets$I, sets$P, sets$Y)))
  )
}

counts_panel <- rbindlist(list(
  counts_from_sets(sets_auto, "AUTOSOME.filtered"),
  counts_from_sets(sets_sex,  "CH12.filtered")
))
fwrite(counts_panel, file.path(OUTDIR, "overlap_counts_byPanel.filtered.tsv"), sep = "\t")