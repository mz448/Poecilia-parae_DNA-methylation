#!/usr/bin/env bash
# DATE: 20261026
# AUTHOR: MZF
# 05_index_positiveSexControl_V01.sh
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#       Build stable identifiers and index DMR positions, then generate pooled
#       VMRs from the Female vs. shuffled_morph dataset using bedtools.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# NOTES:
#       1) This script is specific to pooled shuffled VMRs generated from:
#
#          - shuffled_sex:
#                       f_vs_sM_01
#                       f_vs_sM_02
#                       f_vs_sM_03
#                       f_vs_sM_04
#                       f_vs_sM_05
#                       f_vs_sM_06
#
#       2) DMRs from shuffled_morph, shuffled_all and all real_* conditions are 
#          excluded.
#
#       3) Stable DMR IDs are assigned BEFORE filtering the input dataset.
#          Therefore, DMR IDs correspond to the globally sorted dmr_all.tsv
#          dataset rather than being renumbered within the shuffled subset.
#
#       4) All six shuffled comparisons are pooled BEFORE bedtools merge.
#          Consequently, overlapping DMRs from different shuffled conditions
#          and/or comparisons can contribute to the same VMR.
#
# Purpose:
#   (A) Create a "substrate" TSV containing ONLY shuffled_sex
#       rows, with a stable DMR ID per row. IDs are assigned AFTER sorting the
#       complete dataset by condition, chromosome, coordinates, and comparison,
#       but BEFORE filtering to the selected shuffled conditions.
#   (B) Create an indexed BED for bedtools (chr, start, end, dmr_id).
#   (C) Create pooled VMR regions using bedtools merge, collapsing contributing
#       DMR IDs (union).
#   (D) Add consensus VMR IDs (VMR#######) to the merged regions.
#
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Usage:
#         bash ../scripts/01_DMRs-to-VMRs/01_compileShuffledVMRs_V04.sh dmr_all.tsv .
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
set -euo pipefail


IN_TSV="${1:-}"
OUTDIR="${2:-}"

if [[ -z "${IN_TSV}" || -z "${OUTDIR}" ]]; then
echo "ERROR: missing arguments."
echo "Usage: bash 05_index_positiveSexControl_V01.sh dmr_all.tsv outdir"
exit 1
fi

mkdir -p "${OUTDIR}"

# Outputs
SUBSTRATE_TSV="${OUTDIR}/dmr_shuffled_sex.withIDs.tsv"
INDEXED_BED="${OUTDIR}/dmr_shuffled_sex.positions.bed"
CONSENSUS_BED="${OUTDIR}/vmr_shuffled_sex.bed"

# Temp files
TMP_ALL_SORTED_WITHID="${OUTDIR}/dmr_all.sorted.withID.tsv"
TMP_SHUFFLED_WITHID="${OUTDIR}/.tmp.shuffled_sex.withID.tsv"
TMP_CONS_MERGE="${OUTDIR}/.tmp.shuffled_sex.consensus.merge.bed"

# -----------------------------
# (A) Build substrate TSV:
#     - discard header
#     - sort the complete dataset deterministically
#     - add stable DMR IDs in that global sorted order
#     - then filter to retain ONLY:
#           shuffled_all
#           shuffled_morph
#     - restore header and save as substrate
# -----------------------------

HEADER="$(head -n 1 "${IN_TSV}")"

# Sort keys:
#  - condition (col 13) alphabetically
#  - then chr (col 1), start (col 2), end (col 3), comparison (col 14)
#    for determinism
#
# IMPORTANT:
# DMR IDs are assigned to the complete dmr_all.tsv dataset before filtering.
tail -n +2 "${IN_TSV}" \
| sort -t $'\t' -k13,13 -k1,1 -k2,2n -k3,3n -k14,14 \
| awk -F $'\t' -v OFS=$'\t' '{id=sprintf("DMR%07d", NR); print $0, id}' \
> "${TMP_ALL_SORTED_WITHID}"

# Filter to shuffled_sex conditions ONLY
# (condition is still column 13 here)
awk -F $'\t' \
'$13 == "shuffled_sex"' \
"${TMP_ALL_SORTED_WITHID}" \
> "${TMP_SHUFFLED_WITHID}"

# Restore header and write substrate (adds "dmr_id" as last column)
{
  printf "%s\tdmr_id\n" "${HEADER}"
  cat "${TMP_SHUFFLED_WITHID}"
} > "${SUBSTRATE_TSV}"

# -----------------------------
# (B) Create indexed BED for bedtools
# -----------------------------
# BED columns: chr start end dmr_id (dmr_id is last column, $NF)
#
# DMRs from all six selected shuffled comparisons are pooled here before
# merging. The condition/comparison annotations remain recoverable from the
# substrate TSV through the stable DMR IDs.
awk -F $'\t' -v OFS=$'\t' 'NR>1 {print $1,$2,$3,$NF}' "${SUBSTRATE_TSV}" \
| sort -k1,1 -k2,2n -k3,3n \
> "${INDEXED_BED}"

# -----------------------------
# (C) Merge into pooled VMR regions (union)
#     + collapse contributing DMR IDs
# -----------------------------
# All DMRs from shuffled_sex are merged together.
# Therefore, DMRs from different shuffled conditions/comparisons can
# contribute to the same VMR when their genomic coordinates overlap/touch.
bedtools merge \
-i "${INDEXED_BED}" \
-d 0 \
-c 4 \
-o collapse \
> "${TMP_CONS_MERGE}"
# TMP_CONS_MERGE columns: chr start end collapsed_dmr_ids

# -----------------------------
# (D) Add VMR IDs (VMR#######)
# -----------------------------

awk -F $'\t' -v OFS=$'\t' \
'{cid=sprintf("VMR%07d", NR); print $1,$2,$3,cid,$4}' \
"${TMP_CONS_MERGE}" > "${CONSENSUS_BED}"
# CONSENSUS_BED columns:
# chr start end vmr_id collapsed_dmr_ids

# -----------------------------
# QC checks
# -----------------------------
echo "=== QC: file sanity ==="
echo "Substrate TSV:   ${SUBSTRATE_TSV}"
echo "Indexed BED:     ${INDEXED_BED}"
echo "Pooled VMR BED:  ${CONSENSUS_BED}"
echo

echo "=== QC: conditions present in substrate ==="
echo "Expected: shuffled_sex ONLY"
cut -f13 "${SUBSTRATE_TSV}" | tail -n +2 | sort | uniq -c
echo

echo "=== QC: comparisons present in substrate ==="
echo "Expected exactly six comparisons:"
echo "  f_vs_sM_01"
echo "  f_vs_sM_02"
echo "  f_vs_sM_03"
echo "  f_vs_sM_04"
echo "  f_vs_sM_05"
echo "  f_vs_sM_06"
echo
cut -f14 "${SUBSTRATE_TSV}" | tail -n +2 | sort | uniq -c
echo

# Explicitly verify that ONLY the six expected comparisons are present.
EXPECTED_COMPARISONS="${OUTDIR}/.qc.expected_comparisons.txt"
OBSERVED_COMPARISONS="${OUTDIR}/.qc.observed_comparisons.txt"

printf "%s\n" \
"f_vs_sM_01" \
"f_vs_sM_02" \
"f_vs_sM_03" \
"f_vs_sM_04" \
"f_vs_sM_05" \
"f_vs_sM_06" \
| sort > "${EXPECTED_COMPARISONS}"

cut -f14 "${SUBSTRATE_TSV}" \
| tail -n +2 \
| sort -u \
> "${OBSERVED_COMPARISONS}"

COMPARISON_DIFF=$(comm -3 \
                  "${EXPECTED_COMPARISONS}" \
                  "${OBSERVED_COMPARISONS}" \
                  | wc -l | tr -d ' ')

echo "Comparison-set differences (should be 0): ${COMPARISON_DIFF}"

if [[ "${COMPARISON_DIFF}" -ne 0 ]]; then
echo "ERROR: observed comparisons do not match the six expected comparisons."
echo
echo "Differences:"
comm -3 "${EXPECTED_COMPARISONS}" "${OBSERVED_COMPARISONS}"
exit 1
fi
echo

echo "=== QC: counts ==="
N_SHUFFLED_TSV=$(( $(wc -l < "${SUBSTRATE_TSV}") - 1 ))
N_SHUFFLED_BED=$(wc -l < "${INDEXED_BED}")
N_VMRS=$(wc -l < "${CONSENSUS_BED}")

echo "Selected shuffled DMR rows in substrate (excluding header): ${N_SHUFFLED_TSV}"
echo "Selected shuffled DMR rows in indexed BED:                 ${N_SHUFFLED_BED}"
echo "Pooled VMR regions:                                         ${N_VMRS}"
echo

echo "=== QC: DMR ID uniqueness (should be 0 duplicates) ==="
DUP_IDS=$(cut -f4 "${INDEXED_BED}" | sort | uniq -d | wc -l | tr -d ' ')
echo "Duplicate DMR IDs: ${DUP_IDS}"

if [[ "${DUP_IDS}" -ne 0 ]]; then
echo "First few duplicate IDs:"
cut -f4 "${INDEXED_BED}" | sort | uniq -d | head
exit 1
fi
echo

echo "=== QC: all DMR IDs recovered in pooled VMR collapse lists ==="

# Count total IDs represented across all VMR rows (sum of splits)
N_IDS_IN_VMRS=$(awk -F $'\t' \
                '{n+=split($5,a,",")} END{print n}' \
                "${CONSENSUS_BED}")

echo "Total DMR IDs referenced across VMR rows: ${N_IDS_IN_VMRS}"
echo "Expected (should equal # selected DMRs):  ${N_SHUFFLED_BED}"

if [[ "${N_IDS_IN_VMRS}" -ne "${N_SHUFFLED_BED}" ]]; then
echo "ERROR: total collapsed ID count != number of selected shuffled DMRs."
echo "This suggests missing IDs or unexpected formatting of collapse lists."
exit 1
fi
echo

echo "=== QC: exact ID set match between indexed BED and pooled VMR collapse ==="

cut -f4 "${INDEXED_BED}" \
| sort \
> "${OUTDIR}/.qc.ids_from_bed.txt"

awk -F $'\t' \
'{split($5,a,","); for(i in a) print a[i]}' \
"${CONSENSUS_BED}" \
| sort \
> "${OUTDIR}/.qc.ids_from_vmr.txt"

# comm output should be empty if sets match exactly
DIFF_LINES=$(comm -3 \
             "${OUTDIR}/.qc.ids_from_bed.txt" \
             "${OUTDIR}/.qc.ids_from_vmr.txt" \
             | wc -l | tr -d ' ')

echo "ID set differences (lines in comm -3; should be 0): ${DIFF_LINES}"

if [[ "${DIFF_LINES}" -ne 0 ]]; then
echo "First few differences (unexpected):"
comm -3 \
"${OUTDIR}/.qc.ids_from_bed.txt" \
"${OUTDIR}/.qc.ids_from_vmr.txt" \
| head
exit 1
fi
echo

echo "=== DONE ==="

# Cleanup temp/QC files
rm -f \
"${TMP_SHUFFLED_WITHID}" \
"${TMP_CONS_MERGE}" \
"${EXPECTED_COMPARISONS}" \
"${OBSERVED_COMPARISONS}" \
"${OUTDIR}/.qc.ids_from_bed.txt" \
"${OUTDIR}/.qc.ids_from_vmr.txt"