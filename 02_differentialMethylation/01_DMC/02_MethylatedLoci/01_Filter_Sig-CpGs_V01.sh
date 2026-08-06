#!/bin/bash
# DATE:   20260806
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:   Filters and calculates %mCpG per chromosome from DNMTools output. 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# INPUT:  1) Posterior methylation probability calculated with hmr function in DNMTools 
#         2) A karyotype of CpGs with the same chromosome names as in the samples
# OUTPUT: 2 column files with normalized (%) values of CpG methylation
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail

OUTDIR="methylatedLoci"
mkdir -p "$OUTDIR"

# Keep loci with posterior probability >= 0.9 (90%)
for SAMPLE in postProb_meth.*.tsv; do
BASENAME=$(basename "$SAMPLE" .tsv)

awk '$5 >= 0.9' "$SAMPLE" \
> "${OUTDIR}/${BASENAME}.sig_0.9.tsv"
echo "${SAMPLE} done"
done


# Copy the Reference Chromosome lenght into the folder

cp /local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/karyotype.CpGs.PparFemMitoLambdaPuc19Ver2024.txt${OUTDIR}/.


# Sumarize and normalize counts per chromosome
SCRIPTSFOLDER="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/01_DMC/01_DNMTools/scripts/03_MethylatedLoci"
for SAMPLE in ${OUTDIR}/postProb_meth.*.tsv; do
BASENAME=$(basename "$SAMPLE" .tsv)
echo "${BASENAME} Counting"
  python ${SCRIPTSFOLDER}/99_counts_per_chromosome_03.py ${OUTDIR}/karyotype.CpGs.PparFemMitoLambdaPuc19Ver2024.txt ${SAMPLE} ${OUTDIR}/count.${BASENAME}.tsv
  python ${SCRIPTSFOLDER}/99_normalized_counts_per_chromosome_01.py ${OUTDIR}/karyotype.CpGs.PparFemMitoLambdaPuc19Ver2024.txt ${OUTDIR}/count.${BASENAME}.tsv ${OUTDIR}/norm.${BASENAME}.tsv

done

echo "Done !"