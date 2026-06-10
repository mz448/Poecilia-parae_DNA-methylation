#!/bin/bash
# Copy bismark files that have been proccessed by DNMToosl Sym function
ln -s /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_differential_methylation/01_DNMTools/Sym.forSym.ppar* . 

for file in Sym.forSym.ppar*mem*.CpG_report.txt.meth; do
  # Extract just the sample ID portion
  sample_id=$(echo "$file" | sed 's/^Sym\.forSym\.//; s/\.CpG_report\.txt\.meth$//')
  
  # Define output name
  outfile="sym.formatted.${sample_id}.CpG_report.txt"
  
  # Format and & filter for relevant chromosomes only & save 
  awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $5, $6}' "$file" | awk '$1=="Parae_01" || $1=="Parae_02" || $1=="Parae_03" || $1=="Parae_04" || $1=="Parae_05" || $1=="Parae_06" || $1=="Parae_07" || $1=="Parae_08" || $1=="Parae_09" || $1=="Parae_10" || $1=="Parae_11" || $1=="Parae_12" || $1=="Parae_13" || $1=="Parae_14" || $1=="Parae_15" || $1=="Parae_16" || $1=="Parae_17" || $1=="Parae_18" || $1=="Parae_19" || $1=="Parae_20" || $1=="Parae_21" || $1=="Parae_22" || $1=="Parae_23" || $1=="LambdaNEB" || $1=="P_parae_Mitochondria" || $1=="pUC19" {print $0}' > "$outfile"
done


