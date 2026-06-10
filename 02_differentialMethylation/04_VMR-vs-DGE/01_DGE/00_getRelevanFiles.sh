#!/bin/bash

# DATE:         20250916
# AUTHOR:       MZF
# VERSION:      01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:         Get reference files from 
#                         - Features
#                         - DMRs
#                         - Differential Gene expression Table
#                         - Gene expression CPM (just in case)
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

# Features
ln -s  /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/04_methylation-overlapping/01_geneticFeatures/data/allFeatures.feature.bed .

# DMRs
ln -s /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_differential_methylation/06_cleanParameters_Methylasso_DMR/data/methylasso_DMRs/dmr_all.tsv .

# Differential Gene expression
ln -s /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/06_dmg-vs-deg/00_DEG-Real-AND-Shuffled_02_CPMfilter/data/Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv .

# gene counts 
ln -s /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/06_dmg-vs-deg/00_DEG-Real-AND-Shuffled_02_CPMfilter/data/Combined_Counts.tsv .

# groups mapping
ln -s /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/06_dmg-vs-deg/00_DEG-Real-AND-Shuffled_02_CPMfilter/data/methyl.tsv .

# samples mapping
ln -s /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/06_dmg-vs-deg/00_DEG-Real-AND-Shuffled_02_CPMfilter/data/rnaseq_to_methyl_map.tsv .