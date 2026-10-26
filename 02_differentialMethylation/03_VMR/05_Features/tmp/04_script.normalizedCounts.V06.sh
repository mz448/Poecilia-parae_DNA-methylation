#!/bin/bash

# ==============================================================================
# SCRIPT:     04_script.normalizedCounts.V06.sh
# DATE:       20250707
# AUTHOR:     MZF & ChatGPT
# VERSION:    06
# ==============================================================================
# PURPOSE:
# This script calculates **normalized counts** of genes corresponding to specific 
# genetic features that overlap with **Differentially Methylated Regions (DMRs)**. 
# Normalization is performed by dividing the raw count in a group-specific input 
# file by a reference count for the same gene/feature. 
#
# The script also appends all normalized counts into a consolidated summary file 
# for easier downstream analysis.
#
# ==============================================================================
# INPUTS:
# The script requires the following input files for each **feature** and **group**:
#
# 1. Reference Files:
#    - Format:    perChr.countOfGenes.reference.<feature>.PparFemVer2024.transcriptome.6col.bed
#    - Content:   Contains reference counts per gene/chromosome for each feature.
#
# 2. Group-Specific Files:
#    - Format:    perChr.countPerGene.DMRsOvelapping.<group>.<feature>.PparFemVer2024.transcriptome.6col.bed.tmp
#    - Content:   Raw counts of genes overlapping DMRs for a group and feature.
#
# Requirements:
# - Input files must follow the above naming conventions.
# - Files must be tab-delimited with at least two columns:
#     • Column 1: Gene or chromosome identifier.
#     • Column 2: Count value.
#
# ==============================================================================
# OUTPUTS:
#
# 1. Individual Normalized Count Files:
#    - Format:     normalized.<group>.<feature>.PparFemVer2024.transcriptome.6col.bed
#    - Columns:
#        • Chr              → Gene/Chromosome ID
#        • count            → Raw count
#        • referenceCount   → Value from reference file
#        • normalizedCount  → count / referenceCount
#        • group            → Group name
#        • feature          → Feature name
#
# 2. Summary File:
#    - Format:     summary.normalized.perChr.countPerGene.DMRsOvelapping.allGroups.allFeatures.PparFemVer2024.transcriptome.6col.bed
#    - Description: Aggregated normalized data across all groups and features.
#    - Same column structure as individual files.
#
# 3. Optional Deletion:
#    - Individual normalized files are deleted after consolidation unless the 
#      final `rm` line is commented out.
#
# ==============================================================================
# USAGE:
# 1. Ensure all required input files are in the current directory.
# 2. Run the script:
#       bash 04_script.normalizedCounts.V06.sh
# 3. Output files will be saved in a subdirectory: `normalizedCounts/`
#
# ==============================================================================
# VERSION HISTORY:
# ------------------------------------------------------------------------------
# 02 (20250113)
# - Initial implementation of normalized counts
# - Included features: five_prime_UTR, gene, three_prime_UTR, u2000, d2000
#
# 04 (20250115)
# - Expanded feature list:
#   CDS, d2000, exon, five_prime_UTR_exon, five_prime_UTR, geneBody_exon, 
#   geneBody_intron, geneBody, gene, intron, mRNA, three_prime_UTR_exon, 
#   three_prime_UTR, u2000, UTR_exon, UTR
#
# 05 (20250126)
# - Fully working version
# 
# 06 (20250126)
# - Finalized and documented version
# - Adds structured header documentation for reproducibility and clarity
# ==============================================================================

# ------------------------------- SCRIPT START ----------------------------------

# Features
features=("CDS" "d2000" "exon" "five_prime_UTR_exon" "five_prime_UTR" "geneBody_exon" "geneBody_intron" "geneBody" "gene" "intron" "mRNA" "three_prime_UTR_exon" "three_prime_UTR" "u2000" "UTR_exon" "UTR")

# Groups
groups=("f_vs_i" "f_vs_p" "f_vs_y" "i_vs_y" "p_vs_i" "y_vs_p")

# Directory for outputs
output_dir="normalizedCounts"
mkdir -p "$output_dir"

# Iterate over features
for feature in "${features[@]}"; do
    reference_file="perChr.countOfGenes.reference.${feature}.PparFemVer2024.transcriptome.6col.bed"
    if [[ ! -f "$reference_file" ]]; then
        echo "Reference file $reference_file not found! Skipping..."
        continue
    fi

    for group in "${groups[@]}"; do
        group_file="perChr.countPerGene.DMRsOvelapping.${group}.${feature}.PparFemVer2024.transcriptome.6col.bed.tmp"
        if [[ ! -f "$group_file" ]]; then
            echo "Group file $group_file not found! Skipping..."
            continue
        fi

        output_file="${output_dir}/normalized.${group}.${feature}.PparFemVer2024.transcriptome.6col.bed"

        awk -v GROUP="$group" -v FEATURE="$feature" \
            'NR==FNR {ref[$1]=$2; next} $1 in ref {print $1 "\t" $2 "\t" ref[$1] "\t" $2/ref[$1] "\t" GROUP "\t" FEATURE}' \
            "$reference_file" "$group_file" > "$output_file"

        echo "Normalized counts written to $output_file"
    done
done

cd normalizedCounts/

for FILE in normalized.*; do
    awk '{print $0}' $FILE >> summary.normalized.perChr.countPerGene.DMRsOvelapping.allGroups.allFeatures.PparFemVer2024.transcriptome.6col.bed
done

# Comment the next line to keep the individual files of normalized counts
rm normalized.*

# Print DONE
echo "--------------------------------------------------------------------------------------------------------------"
echo "Run report:"
echo "1) Normalization done for all files (check messages for any errors)"
echo "2) Summary file saved in $output_dir"
