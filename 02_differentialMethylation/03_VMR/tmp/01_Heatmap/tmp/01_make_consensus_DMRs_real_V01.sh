#!/usr/bin/env bash
# DATE: 20260220
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#       Build indexed DMRs + consensus regions using bedtools
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
#
# Purpose:
#   (A) Create a "substrate" TSV that contains ONLY real_* rows, with a stable DMR ID per row,
#       where IDs are assigned AFTER sorting the full dataset by condition (alphabetical; real_* first).
#   (B) Create an indexed BED for bedtools (chr, start, end, dmr_id).
#   (C) Create consensus (union) DMR regions using bedtools merge, collapsing contributing DMR IDs.
#   (D) Add consensus IDs (CONS#######) to the merged regions.
#
# Requirements:
#   - bedtools on PATH
#   - GNU coreutils (sort, awk, comm)
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Usage:
#         bash ../scripts/01_make_consensus_DMRs_real_V01.sh dmr_all.tsv outdir
#         bash ../scripts/01_make_consensus_DMRs_real_V01.sh dmr_all.tsv .
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
set -euo pipefail


IN_TSV="${1:-}"
OUTDIR="${2:-}"

if [[ -z "${IN_TSV}" || -z "${OUTDIR}" ]]; then
  echo "ERROR: missing arguments."
  echo "Usage: bash 01_make_consensus_real_dmrs.sh dmr_all.tsv outdir"
  exit 1
fi

mkdir -p "${OUTDIR}"

# Outputs
SUBSTRATE_TSV="${OUTDIR}/dmr_real.withIDs.tsv"
INDEXED_BED="${OUTDIR}/dmr_real.indexed.bed"
# CONSENSUS_BED="${OUTDIR}/dmr_real_consensus.bed"
CONSENSUS_BED="${OUTDIR}/vmr_real.bed"

# Temp files
TMP_ALL_SORTED_WITHID="${OUTDIR}/dmr_all.sorted.withID.tsv"
TMP_REAL_WITHID="${OUTDIR}/.tmp.real.withID.tsv"
TMP_CONS_MERGE="${OUTDIR}/.tmp.consensus.merge.bed"

# -----------------------------
# (A) Build substrate TSV:
#     - discard header
#     - sort by condition alphabetically (real_* first)
#     - add stable DMR IDs in that sorted order
#     - then filter to keep only real_* rows
#     - restore header and save as substrate
# -----------------------------

HEADER="$(head -n 1 "${IN_TSV}")"

# Sort keys:
#  - condition (col 13) alphabetically  => real_* comes before shuffled_*
#  - then chr (col 1), start (col 2), end (col 3), comparison (col 14) for determinism
tail -n +2 "${IN_TSV}" \
  | sort -t $'\t' -k13,13 -k1,1 -k2,2n -k3,3n -k14,14 \
  | awk -F $'\t' -v OFS=$'\t' '{id=sprintf("DMR%07d", NR); print $0, id}' \
  > "${TMP_ALL_SORTED_WITHID}"

# Filter to real_* conditions ONLY (condition is still column 13 here)
awk -F $'\t' '$13 ~ /^real_/' "${TMP_ALL_SORTED_WITHID}" > "${TMP_REAL_WITHID}"

# Restore header and write substrate (adds "dmr_id" as last column)
{
  printf "%s\tdmr_id\n" "${HEADER}"
  cat "${TMP_REAL_WITHID}"
} > "${SUBSTRATE_TSV}"

# -----------------------------
# (B) Create indexed BED for bedtools
# -----------------------------
# BED columns: chr start end dmr_id (dmr_id is last column, $NF)
awk -F $'\t' -v OFS=$'\t' 'NR>1 {print $1,$2,$3,$NF}' "${SUBSTRATE_TSV}" \
  | sort -k1,1 -k2,2n -k3,3n \
  > "${INDEXED_BED}"

# -----------------------------
# (C) Merge into consensus DMR regions (union) + collapse contributing DMR IDs
# -----------------------------
bedtools merge \
  -i "${INDEXED_BED}" \
  -d 0 \
  -c 4 \
  -o collapse \
  > "${TMP_CONS_MERGE}"
# TMP_CONS_MERGE columns: chr start end collapsed_dmr_ids

# -----------------------------
# (D) Add consensus IDs (CONS#######)
# -----------------------------
awk -F $'\t' -v OFS=$'\t' '{cid=sprintf("CONS%07d", NR); print $1,$2,$3,cid,$4}' \
  "${TMP_CONS_MERGE}" > "${CONSENSUS_BED}"
# CONSENSUS_BED columns: chr start end cons_id collapsed_dmr_ids

# -----------------------------
# QC checks
# -----------------------------
echo "=== QC: file sanity ==="
echo "Substrate TSV:   ${SUBSTRATE_TSV}"
echo "Indexed BED:     ${INDEXED_BED}"
echo "Consensus BED:   ${CONSENSUS_BED}"
echo

echo "=== QC: conditions present in substrate (should be only real_*) ==="
cut -f13 "${SUBSTRATE_TSV}" | tail -n +2 | sort | uniq -c
echo

echo "=== QC: counts ==="
N_REAL_TSV=$(( $(wc -l < "${SUBSTRATE_TSV}") - 1 ))
N_REAL_BED=$(wc -l < "${INDEXED_BED}")
N_CONS=$(wc -l < "${CONSENSUS_BED}")
echo "Real DMR rows in substrate (excluding header): ${N_REAL_TSV}"
echo "Real DMR rows in indexed BED:                 ${N_REAL_BED}"
echo "Consensus merged regions:                     ${N_CONS}"
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

echo "=== QC: all DMR IDs recovered in consensus collapse lists ==="
# Count total IDs represented across all consensus rows (sum of splits)
N_IDS_IN_CONS=$(awk -F $'\t' '{n+=split($5,a,",")} END{print n}' "${CONSENSUS_BED}")
echo "Total DMR IDs referenced across consensus rows: ${N_IDS_IN_CONS}"
echo "Expected (should equal # real DMRs):            ${N_REAL_BED}"
if [[ "${N_IDS_IN_CONS}" -ne "${N_REAL_BED}" ]]; then
  echo "ERROR: total collapsed ID count != number of real DMRs."
  echo "This suggests missing IDs or unexpected formatting of collapse lists."
  exit 1
fi
echo

echo "=== QC: exact ID set match between indexed BED and consensus collapse ==="
cut -f4 "${INDEXED_BED}" | sort > "${OUTDIR}/.qc.ids_from_bed.txt"
awk -F $'\t' '{split($5,a,","); for(i in a) print a[i]}' "${CONSENSUS_BED}" | sort > "${OUTDIR}/.qc.ids_from_consensus.txt"

# comm output should be empty if sets match exactly
DIFF_LINES=$(comm -3 "${OUTDIR}/.qc.ids_from_bed.txt" "${OUTDIR}/.qc.ids_from_consensus.txt" | wc -l | tr -d ' ')
echo "ID set differences (lines in comm -3 (should be empty if sets match exactly)): ${DIFF_LINES}"
if [[ "${DIFF_LINES}" -ne 0 ]]; then
  echo "First few differences (unexpected):"
  comm -3 "${OUTDIR}/.qc.ids_from_bed.txt" "${OUTDIR}/.qc.ids_from_consensus.txt" | head
  exit 1
fi
echo

echo "=== DONE ==="

# Cleanup temp files (alternative)
rm -f "${TMP_REAL_WITHID}" "${TMP_CONS_MERGE}" \
      "${OUTDIR}/.qc.ids_from_bed.txt" "${OUTDIR}/.qc.ids_from_consensus.txt"
