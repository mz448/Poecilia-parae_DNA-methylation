#!/usr/bin/env Rscript
# DATE:       2026-05-11
# AUTHOR:     MZF
# SCRIPT:     10_update_clusterNumbers_fromPlot_V02.R
# VERSION:    02
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# IMPORTANT NOTE:  OPTIONAL!!!
#             This Script is optional. Use only if The clusters where manually 
#             relabeled during figure postprocessing 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#             Add cluster numbers matching the plot labels to VMR cluster files.
#
#             The original cluster files contain cluster IDs from the original
#             cluster lists. The correspondence table maps those original IDs
#             to the cluster numbers used in the manuscript plots.
#
# INPUTS:
#             1. cluster_correspondence_plot-files.tsv
#             2. clusters_AUTOSOME.sampleHeatmap.diffGE_0.4.tsv
#             3. clusters_CH12.sampleHeatmap.diffGE_0.4.tsv
#
# OUTPUTS:
#             1. clusters_AUTOSOME.sampleHeatmap.diffGE_0.4.withPlotClusters.tsv
#             2. clusters_CH12.sampleHeatmap.diffGE_0.4.withPlotClusters.tsv
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# User settings
# ==============================================================================

outdir      <- "VMR_tracks"

correspondence_file <- file.path(outdir,"cluster_correspondence_plot-files.tsv")

autosome_cluster_file <- "clusters_AUTOSOME.sampleHeatmap.diffGE_0.4.tsv"
ch12_cluster_file     <- "clusters_CH12.sampleHeatmap.diffGE_0.4.tsv"

output_autosome_file <- file.path(outdir,"clusters_AUTOSOME.sampleHeatmap.diffGE_0.4.withPlotClusters.tsv")
output_ch12_file     <- file.path(outdir,"clusters_CH12.sampleHeatmap.diffGE_0.4.withPlotClusters.tsv")

# ==============================================================================
# Load correspondence table
# ==============================================================================

corr <- fread(correspondence_file)

required_corr_cols <- c("partition", "Cluster_inplot", "Cluster_inList")
missing_corr_cols <- setdiff(required_corr_cols, names(corr))

if (length(missing_corr_cols) > 0) {
  stop(
    "Missing required columns from correspondence file: ",
    paste(missing_corr_cols, collapse = ", ")
  )
}

corr[, partition := as.character(partition)]
corr[, Cluster_inplot := as.integer(Cluster_inplot)]
corr[, Cluster_inList := as.integer(Cluster_inList)]

corr_clean <- corr[
  ,
  .(
    partition,
    Cluster_inList,
    Cluster_in_plots = Cluster_inplot
  )
]

# ==============================================================================
# Function to add plot cluster numbers
# ==============================================================================

add_plot_clusters <- function(cluster_file, partition_name, output_file) {
  
  dt <- fread(cluster_file)
  
  required_dt_cols <- c("Cluster")
  missing_dt_cols <- setdiff(required_dt_cols, names(dt))
  
  if (length(missing_dt_cols) > 0) {
    stop(
      "Missing required columns from cluster file ",
      cluster_file,
      ": ",
      paste(missing_dt_cols, collapse = ", ")
    )
  }
  
  dt[, Cluster := as.integer(Cluster)]
  dt[, partition := partition_name]
  
  map_dt <- corr_clean[partition == partition_name]
  
  if (nrow(map_dt) == 0) {
    stop("No correspondence rows found for partition: ", partition_name)
  }
  
  dt <- merge(
    dt,
    map_dt,
    by.x = c("partition", "Cluster"),
    by.y = c("partition", "Cluster_inList"),
    all.x = TRUE,
    sort = FALSE
  )
  
  missing_map <- dt[is.na(Cluster_in_plots), unique(Cluster)]
  
  if (length(missing_map) > 0) {
    stop(
      "Some clusters in ",
      cluster_file,
      " did not have a correspondence in the mapping file: ",
      paste(missing_map, collapse = ", ")
    )
  }
  
  dt[, partition := NULL]
  
  base_cols <- names(dt)
  base_cols <- base_cols[base_cols != "Cluster_in_plots"]
  
  cluster_pos <- which(base_cols == "Cluster")
  
  final_cols <- append(
    base_cols,
    values = "Cluster_in_plots",
    after = cluster_pos
  )
  
  setcolorder(dt, final_cols)
  
  fwrite(
    dt,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  
  message("Saved: ", output_file)
  message("Cluster mapping used for ", partition_name, ":")
  
  print(
    unique(dt[, .(Cluster, Cluster_in_plots)])[order(Cluster)]
  )
}

# ==============================================================================
# Add plot cluster numbers and save new files
# ==============================================================================

add_plot_clusters(
  cluster_file = autosome_cluster_file,
  partition_name = "AUTO",
  output_file = output_autosome_file
)

add_plot_clusters(
  cluster_file = ch12_cluster_file,
  partition_name = "CH12",
  output_file = output_ch12_file
)