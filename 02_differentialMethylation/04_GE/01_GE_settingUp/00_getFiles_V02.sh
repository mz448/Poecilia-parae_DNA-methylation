#!bin/bash
# DATE:       2026-0831
# AUTHOR:     MZF 
# SCRIPT:     01_getFiles_V02.R
# VERSION:    02
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:       To call Gene Expression files form the Flutamide Dataset.
#             These Files where generated and processed by Ben and Ehren. 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

# Data Folder (Already on this project's Data Folder)
GEDATA='/local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/gene_expression'

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Flutamide Study Samples
cp $GEDATA/Sample_Info.tsv Sample_Info.tsv
# Gene counts from all samples already extracted by Ehren
cp $GEDATA/Combined_Counts.tsv Combined_Counts.tsv
# Names of Genes in the ppar Fem genome
cp $GEDATA/gene_names.tsv gene_names.tsv 

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
## Data Processed By Ehren
## Differential Gene expression
# cp $GEDATA/Muscle_combined_DEGs.tsv .

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Files I made manually to make both datasets  (meth & flurtamideGE) compatible.
# 1) methyl.tsv <-  Name and grouping of real and shuffled samples from methylations analysis
cp $GEDATA/methyl.tsv .
# 2) rnaseq_to_methyl_map.tsv <-  Map of artificial correspondence between methylation and RNA seq samples 
cp $GEDATA/rnaseq_to_methyl_map.tsv .
# 3) Sample tsv necessary for Making the DEseq analysis compatible with Methylations.
cp $GEDATA/sample_info_with_comparisons.tsv .
# 4) Order of contrasts
cp $GEDATA/orderOfcontrasts.tsv .