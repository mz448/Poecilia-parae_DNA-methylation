#!/usr/bin/env bash
# DATE:       2025-09-31
# AUTHOR:     MZF
# SCRIPT:     01_prepDMRs_02_filterByAutosomes_V01.sh
# VERSION:    01
# GOAL:       Filter for DMRs in Autosomes only
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# USAGE:
#             
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# NOTES:
#             
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
INTERSECT="intersect_features_vs_dmrs/intersect_features_vs_dmrs.tsv"
# Filter Buy autosome Only
awk -F'\t' '$1 =="feat_chr" || $1=="Parae_01" || $1=="Parae_02" || $1=="Parae_03" || $1=="Parae_04" || $1=="Parae_05" || $1=="Parae_06" || $1=="Parae_07" || $1=="Parae_08" || $1=="Parae_09" || $1=="Parae_10" || $1=="Parae_11" || $1=="Parae_13" || $1=="Parae_14" || $1=="Parae_15" || $1=="Parae_16" || $1=="Parae_17" || $1=="Parae_18" || $1=="Parae_19" || $1=="Parae_20" || $1=="Parae_21" || $1=="Parae_22" || $1=="Parae_23" || $1=="LambdaNEB" || $1=="P_parae_Mitochondria" || $1=="pUC19" {print $0}' "$INTERSECT" > intersect_features_vs_dmrs/AUTO_intersect_features_vs_dmrs.tsv
# Check Autosome Only
awk '{print $1}' intersect_features_vs_dmrs/AUTO_intersect_features_vs_dmrs.tsv | sort | uniq
