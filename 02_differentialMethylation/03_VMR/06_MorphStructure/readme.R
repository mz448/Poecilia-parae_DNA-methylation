DATE:   20260803
AUTHOR: MZF

GOAL:
    To determine whether methylation patterns across VMRs are structured by 
    biological phenotype, sex, or male-morph identity, rather than arising 
    primarily from individual sample variability and group assignments.
CONCEPTUALLY:
    Do samples from the same biological identity have more similar VMR methylation 
    profiles than samples randomly selected?
SCRIPTS:
1)  01_compile_conditionSpecific_VMRs_V06.sh
    Filters previously indexed DMRs by condition and generates VMRs
    using bedtools while preserving the existing stable DMR IDs.
2)  02_calculateMethylation_V08.sh
    Compute per-VMR mean methylation and mean coverageFrom existing formatted
    symmetric CpG report files and a condition-specific VMR file
3)  03_PERMANOVA_V02.R
    Test whether the effect of sample identity (in a VMR dataset) is 
    relevant for the methylation structure 
METHODS:
    For each VMR set, we calculated mean methylation values for each 
    individual VMR in all individual samples. Then Euclidean distances 
    were calculated between sample methylation profiles, and PERMANOVA 
    was used to test three factors: overall phenotype among the four 
    groups (~phenotype), sex across all samples (~sex), and male morph 
    among males only (~male-morph). VMR datasets were analyzed separately 
    and compared using their (R^2), pseudo-(F), permutation p-values, 
    and dispersion results. PERMDISP was used alongside the sex and 
    male-morph tests to assess whether significant PERMANOVA results 
    reflected differences in group centroids or unequal within-group 
    dispersion. 