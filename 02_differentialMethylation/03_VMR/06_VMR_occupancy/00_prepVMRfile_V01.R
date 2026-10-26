#!/usr/bin/env bash
# DATE:     20260315
# AUTHOR:   MZF
# UPDATED:  20260829
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#           Merge VMR files from AUTO and CH12 into WG
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
THRESHOLD=0

awk 'BEGIN {OFS="\t"; print "chr","start","end","VMR_Cluster","score","strand"}
     FNR > 1 {
         print $1, $2, $3, $4"|Cluster_"$6, $5, "+"
     }' \
clusters_AUTOSOME.sampleHeatmap.diffGE_$THRESHOLD.tsv \
clusters_CH12.sampleHeatmap.diffGE_$THRESHOLD.tsv \
> clusters_WG.sampleheatmap.diffGE_$THRESHOLD.bed