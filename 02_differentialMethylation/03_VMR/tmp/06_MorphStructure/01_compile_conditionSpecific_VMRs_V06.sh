#!/usr/bin/env bash
# DATE: 20260220
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#       Filter previously indexed DMRs by condition and generate VMRs
#       using bedtools while preserving the existing stable DMR IDs.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
#
# Purpose:
#   (A) Filter a previously sorted and indexed DMR file to retain only rows
#       matching CONDITION while preserving the existing DMR IDs.
#   (B) Create an indexed BED for bedtools:
#       chr, start, end, dmr_id.
#   (C) Create VMRs using bedtools merge, collapsing contributing DMR IDs.
#   (D) Add stable VMR identifiers in the form VMR#######.
#
# Input:
#   dmr_all.sorted.withID.tsv
#
#   This file:
#     - has already been sorted;
#     - already contains stable DMR IDs in its final column;
#     - does not contain a header.
#
# Requirements:
#   - bedtools on PATH
#   - GNU coreutils: sort, awk, comm
#
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Usage:
#   bash ../scripts/01_compile_conditionSpecific_VMRs_V06 dmr_all.sorted.withID.tsv outdir
#   bash ../scripts/01_compile_conditionSpecific_VMRs_V06 dmr_all.sorted.withID.tsv .
#
# Defaults:
#   bash ../scripts/01_compile_conditionSpecific_VMRs_V06
#
# V04: Updated object names to VMR terminology.
# V05: Added CONDITION-based exact filtering to compile VMRs separately from
#      real_morph or shuffled_morph DMRs.
# V06: Uses the previously generated dmr_all.sorted.withID.tsv as input,
#      preserving the existing stable DMR IDs in condition-specific VMR sets.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail


# -----------------------------
# Global parameters
# -----------------------------

IN_TSV="${1:-dmr_all.sorted.withID.tsv}"
OUTDIR="${2:-dataset_Structure}"

if [[ ! -f "${IN_TSV}" ]]; then
  echo "ERROR: input file not found: ${IN_TSV}"
  echo "Usage: bash 01_compileVMRs_V06.sh dmr_all.sorted.withID.tsv outdir"
  exit 1
fi

mkdir -p "${OUTDIR}"


# -----------------------------
# Select one condition
# -----------------------------

CONDITION="real_morph"
# CONDITION="shuffled_morph"
# CONDITION="real_sex"
# CONDITION="shuffled_sex"
# CONDITION="shuffled_all"


# -----------------------------
# Output files
# -----------------------------

SUBSTRATE_TSV="${OUTDIR}/dmr_${CONDITION}.withIDs.tsv"
INDEXED_BED="${OUTDIR}/dmr_${CONDITION}.positions.bed"
VMR_BED="${OUTDIR}/vmr_${CONDITION}.bed"


# -----------------------------
# Temporary files
# -----------------------------

TMP_CONDITION_WITHID="${OUTDIR}/.tmp.${CONDITION}.withID.tsv"
TMP_VMR_MERGE="${OUTDIR}/.tmp.${CONDITION}.vmr.merge.bed"


# -----------------------------
# (A) Build condition-specific substrate TSV
#
# The input file:
#   - is already sorted;
#   - already contains stable DMR IDs;
#   - does not contain a header.
#
# Therefore, rows are filtered directly without
# sorting or assigning new DMR IDs.
# -----------------------------

HEADER=$'chr\tstart\tend\tnum.cpgs1\tnum.cpgs2\tcov.score\tmeth1\tmeth2\tdiff\tpvalue\tFDR\tannotation\tcondition\tcomparison\tchromosome_type\tdmr_id'

# Keep only rows whose condition column exactly matches CONDITION.
awk -F $'\t' -v COND="${CONDITION}" '$13 == COND' \
  "${IN_TSV}" \
  > "${TMP_CONDITION_WITHID}"

# Restore the header while preserving the existing DMR IDs.
{
  printf "%s\n" "${HEADER}"
  cat "${TMP_CONDITION_WITHID}"
} > "${SUBSTRATE_TSV}"


# -----------------------------
# (B) Create indexed BED
#
# BED columns:
#   chr, start, end, dmr_id
#
# The existing dmr_id is the final column.
# -----------------------------

awk -F $'\t' -v OFS=$'\t' \
  'NR > 1 {print $1, $2, $3, $NF}' \
  "${SUBSTRATE_TSV}" \
  | sort -k1,1 -k2,2n -k3,3n \
  > "${INDEXED_BED}"


# -----------------------------
# (C) Merge DMRs into VMRs
#
# Collapse the IDs of all DMRs contributing
# to each merged VMR.
# -----------------------------

bedtools merge \
  -i "${INDEXED_BED}" \
  -d 0 \
  -c 4 \
  -o collapse \
  > "${TMP_VMR_MERGE}"

# TMP_VMR_MERGE columns:
#   chr, start, end, collapsed_dmr_ids


# -----------------------------
# (D) Add VMR identifiers
# -----------------------------

awk -F $'\t' -v OFS=$'\t' \
  '{
     vmr_id = sprintf("VMR%07d", NR)
     print $1, $2, $3, vmr_id, $4
   }' \
  "${TMP_VMR_MERGE}" \
  > "${VMR_BED}"

# VMR_BED columns:
#   chr, start, end, vmr_id, collapsed_dmr_ids


# -----------------------------
# QC checks
# -----------------------------

echo "=== QC: file sanity ==="
echo "Input DMR file: ${IN_TSV}"
echo "Condition:      ${CONDITION}"
echo "Substrate TSV:  ${SUBSTRATE_TSV}"
echo "Indexed BED:    ${INDEXED_BED}"
echo "VMR BED:        ${VMR_BED}"
echo


echo "=== QC: conditions present in substrate ==="
echo "Expected condition: ${CONDITION}"
cut -f13 "${SUBSTRATE_TSV}" \
  | tail -n +2 \
  | sort \
  | uniq -c
echo


echo "=== QC: counts ==="

N_DMR_TSV=$(( $(wc -l < "${SUBSTRATE_TSV}") - 1 ))
N_DMR_BED=$(wc -l < "${INDEXED_BED}")
N_VMR=$(wc -l < "${VMR_BED}")

echo "DMR rows in substrate, excluding header: ${N_DMR_TSV}"
echo "DMR rows in indexed BED:                 ${N_DMR_BED}"
echo "Merged VMR regions:                      ${N_VMR}"
echo


echo "=== QC: DMR ID uniqueness ==="

DUP_IDS=$(
  cut -f4 "${INDEXED_BED}" \
    | sort \
    | uniq -d \
    | wc -l \
    | tr -d ' '
)

echo "Duplicate DMR IDs, expected 0: ${DUP_IDS}"

if [[ "${DUP_IDS}" -ne 0 ]]; then
  echo "ERROR: duplicate DMR IDs detected."
  echo "First few duplicate IDs:"

  cut -f4 "${INDEXED_BED}" \
    | sort \
    | uniq -d \
    | head

  exit 1
fi

echo


echo "=== QC: all DMR IDs recovered in VMR collapse lists ==="

N_IDS_IN_VMRS=$(
  awk -F $'\t' \
    '{
       n += split($5, ids, ",")
     }
     END {
       print n + 0
     }' \
    "${VMR_BED}"
)

echo "Total DMR IDs referenced across VMRs: ${N_IDS_IN_VMRS}"
echo "Expected number of DMR IDs:           ${N_DMR_BED}"

if [[ "${N_IDS_IN_VMRS}" -ne "${N_DMR_BED}" ]]; then
  echo "ERROR: collapsed DMR ID count does not equal the number of input DMRs."
  echo "This suggests missing IDs or unexpected collapse-list formatting."
  exit 1
fi

echo


echo "=== QC: exact DMR ID set match ==="

QC_IDS_FROM_BED="${OUTDIR}/.qc.${CONDITION}.ids_from_bed.txt"
QC_IDS_FROM_VMRS="${OUTDIR}/.qc.${CONDITION}.ids_from_vmrs.txt"

cut -f4 "${INDEXED_BED}" \
  | sort \
  > "${QC_IDS_FROM_BED}"

awk -F $'\t' \
  '{
     split($5, ids, ",")
     for (i in ids) {
       print ids[i]
     }
   }' \
  "${VMR_BED}" \
  | sort \
  > "${QC_IDS_FROM_VMRS}"

DIFF_LINES=$(
  comm -3 "${QC_IDS_FROM_BED}" "${QC_IDS_FROM_VMRS}" \
    | wc -l \
    | tr -d ' '
)

echo "ID set differences, expected 0: ${DIFF_LINES}"

if [[ "${DIFF_LINES}" -ne 0 ]]; then
  echo "ERROR: indexed BED and VMR collapse lists contain different DMR ID sets."
  echo "First few differences:"

  comm -3 "${QC_IDS_FROM_BED}" "${QC_IDS_FROM_VMRS}" \
    | head

  exit 1
fi

echo


echo "=== DONE ==="


# -----------------------------
# Cleanup temporary files
# -----------------------------

rm -f \
  "${TMP_CONDITION_WITHID}" \
  "${TMP_VMR_MERGE}" \
  "${QC_IDS_FROM_BED}" \
  "${QC_IDS_FROM_VMRS}"