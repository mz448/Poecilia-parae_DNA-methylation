#!/bin/bash

# DATE:         20260801
# AUTHOR:       MZF
# VERSION:      01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:         Get reference files from
#                         - Features
#                         - DMRs
#                         - Differential Gene expression Table
#                         - Gene expression CPM (just in case)
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

SOURCE_FOLDER="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/04_GE/data/00_GE-and-DGE"

# Features
ln -s  ${SOURCE_FOLDER}/allFeatures.feature.bed .

# DMRs
ln -s /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/02_DMR/data/01_methylasso/methylasso_DMRs/dmr_all.tsv .

# Differential Gene expression
ln -s ${SOURCE_FOLDER}/Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv .

# gene counts
ln -s ${SOURCE_FOLDER}/Combined_Counts.tsv .

# groups mapping
ln -s ${SOURCE_FOLDER}/methyl.tsv .

# samples mapping
ln -s ${SOURCE_FOLDER}/rnaseq_to_methyl_map.tsv .

# Additional Files:
# VMRs & distance to genes
cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/plots/0/VMRs_diff0.closest_mRNA.tsv ./VMR_diff0.closest_mRNA.tsv
# OVMRS overlapping with Genes
cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/plots/0/dist_VMR-to-gene/VMR_diff0.overlap_to_mRNA.tsv .
# Only intergenic VMRs with their distance
cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/plots/0/dist_VMR-to-gene/VMR_diff0.intergenic.tsv .

cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/plots/0/vmr.morph_and_sample_observed.filtered.diffGE_0.tsv .