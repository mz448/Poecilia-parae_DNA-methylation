#!/bin/bash
# Link files with DMRs
ln -s /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_differential_methylation/02_DSS/default_DMR_combined.bed /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/04_methylation-overlapping/98_gff3_2_features_BedFiles

# sort the reference file and crop the header
awk  '{print $0}' default_DMR_combined.bed | tail -n +2 | sort -k1,1V -k2,2n > sorted.default_DMR_combined.bed

# Check groups 
awk  '{print $10}' sorted.default_DMR_combined.bed | sort | uniq 

# Filter in individual files per group
for GROUP in f_vs_i f_vs_p f_vs_y i_vs_y p_vs_i y_vs_p;
do
	awk -v "GROUP=$GROUP" -F'[\t]' '$10 ==GROUP {print $0}' sorted.default_DMR_combined.bed > $GROUP.sorted.default_DMR_combined.bed
	done
