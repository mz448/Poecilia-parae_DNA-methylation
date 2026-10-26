#!/usr/bin/env bash
# DATE: 20260923
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#      Compile Reference Files to perform the VMR analysis 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# NOTE: Run from 03_VMR/data 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

# All DMRs 
cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/02_DMR/data/01_methylasso/methylasso_DMRs/dmr_all.tsv .

# DNA methylation counts files (symetric) 
cp /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/01_DMC/01_DNMTools/data/Sym.forSym.ppar*.CpG_report.txt.meth .