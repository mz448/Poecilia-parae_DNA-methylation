#!/bin/bash
# DATE:       20250907
# AUTHOR:     MZF
# SCRIPT:     00_getRelevantFiles.sh
# ~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%
# GOAL:       Copy & reformat Bismark CpG_report.meth files processed by DNMTools Sym
#             into methylasso-compatible format, filtering to nuclear Chromosomes only.
# ~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%~~~~~%~~~~~~%~~~~~%

SOURCE="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/01_DMC/01_DNMTools/data"
LOGDIR="logs"
LOGFILE="${LOGDIR}/00_getRelevantFiles.log"

# Make sure log dir exists
mkdir -p "$LOGDIR"
echo "==== Preparing Files: $(date) ====" > "$LOGFILE"
echo "SOURCE FOLDER: $SOURCE" >> "$LOGFILE"

for FILE in "$SOURCE"/Sym.forSym.ppar*mem*.CpG_report.txt.meth; do
  # Clean basename
  base=$(basename "$FILE")                      # Sym.forSym.pparfmem001.CpG_report.txt.meth
  base=${base#Sym.forSym.ppar}                  # pparfmem001.CpG_report.txt.meth
  base=${base%.CpG_report.txt.meth}             # pparfmem001

  # Define output name
  outfile="methylasso.${base}.txt"

  # Format for methylasso (chr start end mC covC) and filter to real chromosomes
  awk 'BEGIN{OFS="\t"} {print $1, $2, $2+1, $5, $6}' "$FILE" |
    awk '$1 ~ /^Parae_[0-9]+$/ && $1 <= "Parae_23"' > "$outfile"

  # Log result
  echo "Sample -> $outfile [✓]" >> "$LOGFILE"
done

echo "==== Finished run: $(date) ====" >> "$LOGFILE"


