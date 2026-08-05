#!/bin/bash

# Link files with DMRs
cp /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_differential_methylation/02_DSS/01_DSS-DifferentialMethylation/Routput/DATA/DMRs.compiled.bed .
# sort the reference file and crop the header
awk  '{print $0}' DMRs.compiled.bed | tail -n +2 | sort -k1,1V -k2,2n > DMRs.bed

# Check groups 
awk  '{print $10}' DMRs.bed | sort | uniq 

# Filter in individual files per group
for GROUP in f_vs_i f_vs_p f_vs_y i_vs_y_01 p_vs_i_01 y_vs_p_01;
do
	awk -v "GROUP=$GROUP" -F'[\t]' '$10 ==GROUP {print $0}' DMRs.bed > $GROUP.real.DMRs.bed
done


# Filter in individual files per group
for GROUP in rf_vs_sM_01 rf_vs_sM_02 rf_vs_sM_03 rf_vs_sM_04 rf_vs_sM_05 rf_vs_sM_06;
do
	awk -v "GROUP=$GROUP" -F'[\t]' '$10 ==GROUP {print $0}' DMRs.bed > $GROUP.shSex.DMRs.bed
done

# Filter in individual files per group
for GROUP in sM_vs_sM_01 sM_vs_sM_02 sM_vs_sM_03;
do
	awk -v "GROUP=$GROUP" -F'[\t]' '$10 ==GROUP {print $0}' DMRs.bed > $GROUP.shMorph.DMRs.bed
done

# Filter in individual files per group
for GROUP in sA_vs_sA_01 sA_vs_sA_02 sA_vs_sA_03;
do
	awk -v "GROUP=$GROUP" -F'[\t]' '$10 ==GROUP {print $0}' DMRs.bed > $GROUP.shAll.DMRs.bed
done
