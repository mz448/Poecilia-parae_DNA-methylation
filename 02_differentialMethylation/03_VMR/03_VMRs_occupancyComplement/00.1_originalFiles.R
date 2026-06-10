#!/bin/bash
# DATE: 20260214     
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#       List of original files:
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%


  # metrics 
  - promoter_CpG_metrics_u400-TSS-d200.tsv
  # promoter classes (use class_2)
  - promoterClasses_u400-TSS-d200.tsv
  # expression data
  - Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv
  # CGIs
  - CGI-Pparae.cleaned.sorted.txt

  # Diferential expression data
  Muscle_real-and-shuffled_combined_DEGs_CPMfilter_labeled.tsv
  # promoter metrics and positions
  promoter_CpG_metrics_u400-TSS-d200.tsv
  # promoter clasifications
  promoterClasses_u400-TSS-d200.tsv
  # cDMRs filtered to threshold of min difference + methylaiton values per sample
  dmr_real_consensus.morph_and_sample_observed.tsv
  # Consensus BED
  dmr_real_consensus.sorted.bed