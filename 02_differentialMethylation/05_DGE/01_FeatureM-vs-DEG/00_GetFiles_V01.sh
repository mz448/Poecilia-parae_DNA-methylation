#!/bin/bash

# DATE: 20260809
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:         Get reference files from
#                         - Methylation at Genetic features
#                         - VMRs overlaping with Genes
#                         - Differential Gene expression Table
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# METHYILATION AT FEATURES
cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/04_GE/data/00_GE-and-DGE/gene_methylation_means/final_tables/meth.combined.tsv .

# DEG_FILE
cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/04_GE/data/00_GE-and-DGE/Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv .

# VMR_FILE
cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/05_DGE/data/VMR_diff0.overlap_to_mRNA.tsv .