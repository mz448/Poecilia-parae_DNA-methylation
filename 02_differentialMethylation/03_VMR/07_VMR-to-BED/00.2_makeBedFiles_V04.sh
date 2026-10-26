#!/usr/bin/env bash
# DATE: 20260314
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Make BED files from consensus DMR sets and recover the matching
#   rows from the morph_and_sample_observed table.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail

DIFF_THRESH="0"
OUTDIR="../plots_${DIFF_THRESH}"
mkdir ${OUTDIR}

AUTO_TSV="clusters_AUTOSOME.sampleHeatmap.diffGE_${DIFF_THRESH}.tsv"
CH12_TSV="clusters_CH12.sampleHeatmap.diffGE_${DIFF_THRESH}.tsv"

AUTO_BED="${OUTDIR}/clusters_AUTOSOME.sampleHeatmap.diffGE_${DIFF_THRESH}.bed"
CH12_BED="${OUTDIR}/clusters_CH12.sampleHeatmap.diffGE_${DIFF_THRESH}.bed"

STRICT_BED="cDMRs_diff${DIFF_THRESH}.bed" 

OBSERVED_TSV="dmr_real_consensus.morph_and_sample_observed.tsv"
STRICT_OBSERVED_TSV="${OUTDIR}/consensus-DMRs_diff${DIFF_THRESH}.morph_and_sample_observed.tsv"
STRICT_AUTO_TSV="${OUTDIR}/consensus-DMRs_diff${DIFF_THRESH}.AUTOSOME.morph_and_sample_observed.tsv"
STRICT_CH12_TSV="${OUTDIR}/consensus-DMRs_diff${DIFF_THRESH}.CH12.morph_and_sample_observed.tsv"

# ------------------------------------------------------------
# 1) Convert cluster TSVs to BED for IGV
# ------------------------------------------------------------
awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1,$2,$3,$4"|Cluster_"$6,$5,"."}' \
  "${AUTO_TSV}" > "${AUTO_BED}"

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $1,$2,$3,$4"|Cluster_"$6,$5,"."}' \
  "${CH12_TSV}" > "${CH12_BED}"

# ------------------------------------------------------------
# 2) Make a single BED for convenience
# ------------------------------------------------------------
cat "${AUTO_BED}" "${CH12_BED}" > "${STRICT_BED}"

# ------------------------------------------------------------
# 3) Filter the morph_and_sample_observed table using the strict BED
# ------------------------------------------------------------
awk '
BEGIN{FS=OFS="\t"}
NR==FNR {
  split($4, a, "|")
  keep[a[1]] = 1
  next
}
FNR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="cons_id") cons_col=i
  }
  print
  next
}
($cons_col in keep)
' "${STRICT_BED}" "${OBSERVED_TSV}" > "${STRICT_OBSERVED_TSV}"

# ------------------------------------------------------------
# 4) Split filtered table by chromosome_type
# ------------------------------------------------------------
awk '
BEGIN{FS=OFS="\t"}
FNR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="chromosome_type") chrom_col=i
  }
  print
  next
}
$chrom_col=="Autosome"
' "${STRICT_OBSERVED_TSV}" > "${STRICT_AUTO_TSV}"

awk '
BEGIN{FS=OFS="\t"}
FNR==1 {
  for (i=1; i<=NF; i++) {
    if ($i=="chromosome_type") chrom_col=i
  }
  print
  next
}
$chrom_col=="Sex_Ch"
' "${STRICT_OBSERVED_TSV}" > "${STRICT_CH12_TSV}"

echo "DONE"
echo "BED files:"
echo "  ${AUTO_BED}"
echo "  ${CH12_BED}"
echo "  ${STRICT_BED}"
echo
echo "Filtered tables:"
echo "  ${STRICT_OBSERVED_TSV}"
echo "  ${STRICT_AUTO_TSV}"
echo "  ${STRICT_CH12_TSV}"