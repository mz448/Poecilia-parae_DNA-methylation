#!/usr/bin/env bash
# DATE: 20260315
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Calculate CpG coverage for VMRs using each sym CpG report as
#   the reference list of CpGs for each sample.
#
#   For each  DMR and each sample, report:
#     - total CpGs in the region
#     - covered CpGs in the region
#     - percent covered
#
#   Then summarize coverage separately for:
#     - Autosomes
#     - Sex chromosome (Parae_12)
#
# DETAILS OF COLUMNS:
#   total_CpGs: total number of CpGs in the reference sym file for that chromosome class
#   n_regions: number of  DMRs in that chromosome class
#   total_CpGs_inRegions: total number of CpGs falling inside the  DMRs
#   pct_overlap: percent of all CpGs in that chromosome class that fall inside  DMRs
#   covered_CpGs_inSample: number of CpGs inside  DMRs with sample coverage > 0
#   pct_covered_inSample: percent of CpGs inside  DMRs that are covered in the sample
#
# NOTES:
#   - Run this script from the data/ directory
#   - Uses only columns 1-3 (chr, start, end) from the sym file to define CpGs
#   - "Covered CpGs" are CpGs with coverage > 0 in column 5 of the sym file
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail
shopt -s nullglob



DIFF_THRESH="0"
OUTDIR="../plots/${DIFF_THRESH}/coverage"
DMR_BED="../plots/${DIFF_THRESH}/clusters_WG.sampleheatmap.diffGE_${DIFF_THRESH}.bed"

# DMR_BED="VMRs_diff${DIFF_THRESH}.bed"
SYM_GLOB="sym.formatted.ppar*mem*.CpG_report.txt"

PER_REGION_DIR="${OUTDIR}/coverage_VMRs/per_region"
SUMMARY_DIR="${OUTDIR}/coverage_VMRs/summary"
MASTER_SUMMARY="${OUTDIR}/coverage_VMRs/coverage.dmr_real_consensus.summary.ALLsamples.tsv"

CHROM_RE='^Parae_(0[1-9]|1[0-9]|2[0-3])$'

mkdir -p "${OUTDIR}" "${PER_REGION_DIR}" "${SUMMARY_DIR}"

TMPDIR="$(mktemp -d)"
trap 'rm -rf "${TMPDIR}"' EXIT

_SORTED="${TMPDIR}/.sorted.bed"

# ------------------------------------------------------------
# 0) Basic checks
# ------------------------------------------------------------
[[ -f "${DMR_BED}" ]] || { echo "ERROR: Missing  BED: ${DMR_BED}"; exit 1; }

sym_files=( ${SYM_GLOB} )
if [[ ${#sym_files[@]} -eq 0 ]]; then
  echo "ERROR: No sym files matched: ${SYM_GLOB}"
  exit 1
fi

# ------------------------------------------------------------
# 1) Sort  BED once
# ------------------------------------------------------------
LC_ALL=C sort -k1,1 -k2,2n -k3,3n "${DMR_BED}" > "${_SORTED}"

# ------------------------------------------------------------
# 2) Loop over all samples
# ------------------------------------------------------------
master_tmp="${TMPDIR}/master_summary.tsv"
echo -e "sample\tchromosome_type\ttotal_CpGs\tn_regions\ttotal_CpGs_inRegions\tpct_overlap\tcovered_CpGs_inSample\tpct_covered_inSample" > "${master_tmp}"

for REF_SYM in "${sym_files[@]}"; do
  SAMPLE_ID="$(basename "${REF_SYM}" .CpG_report.txt)"
  SAMPLE_ID="${SAMPLE_ID#sym.formatted.}"

  OUT_PER_REGION="${PER_REGION_DIR}/coverage.dmr_real_consensus..${SAMPLE_ID}.perRegion.tsv"
  OUT_SUMMARY="${SUMMARY_DIR}/coverage.dmr_real_consensus..${SAMPLE_ID}.summary.tsv"

  ALL_CPGS="${TMPDIR}/${SAMPLE_ID}.all_cpgs.bed"
  COVERED_CPGS="${TMPDIR}/${SAMPLE_ID}.covered_cpgs.bed"
  ALL_COUNTS="${TMPDIR}/${SAMPLE_ID}..totalCounts.bed"
  COV_COUNTS="${TMPDIR}/${SAMPLE_ID}..coveredCounts.bed"
  GENOME_TOTALS="${TMPDIR}/${SAMPLE_ID}.genome_totals.tsv"

  # ----------------------------------------------------------
  # 2A) Build CpG BEDs from sym file
  # ----------------------------------------------------------
  awk -F $'\t' -v OFS=$'\t' -v re="${CHROM_RE}" '
    ($1 ~ re) { print $1,$2,$3 }
  ' "${REF_SYM}" \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
  > "${ALL_CPGS}"

  awk -F $'\t' -v OFS=$'\t' -v re="${CHROM_RE}" '
    ($1 ~ re) && ($5+0 > 0) { print $1,$2,$3 }
  ' "${REF_SYM}" \
  | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
  > "${COVERED_CPGS}"

  # ----------------------------------------------------------
  # 2B) Count overlaps per  DMR
  # ----------------------------------------------------------
  bedtools intersect -a "${_SORTED}" -b "${ALL_CPGS}"     -c > "${ALL_COUNTS}"
  bedtools intersect -a "${_SORTED}" -b "${COVERED_CPGS}" -c > "${COV_COUNTS}"

  # ----------------------------------------------------------
  # 2C) Build per-region output
  # ----------------------------------------------------------
  {
    echo -e "sample\tchr\tstart\tend\tname\tscore\tstrand\tcons_id\tcluster\tchromosome_type\ttotal_CpGs\tcovered_CpGs\tpct_covered"

    paste "${ALL_COUNTS}" "${COV_COUNTS}" \
    | awk -F $'\t' -v OFS=$'\t' -v sample="${SAMPLE_ID}" '
        {
          chr    = $1
          start  = $2
          end    = $3
          name   = $4
          score  = $5
          strand = $6

          total   = $7
          covered = $14

          split(name, a, "|")
          cons_id = a[1]
          cluster = a[2]

          chrom_type = (chr=="Parae_12" ? "Sex_Ch" : "Autosome")
          pct = (total > 0 ? 100 * covered / total : "NA")

          print sample, chr, start, end, name, score, strand, cons_id, cluster, chrom_type, total, covered, pct
        }
      '
  } > "${OUT_PER_REGION}"

  # ----------------------------------------------------------
  # 2D) Genome-wide total CpGs by chromosome class
  # ----------------------------------------------------------
  {
    echo -e "chromosome_type\ttotal_CpGs"
    awk -F $'\t' -v OFS=$'\t' '
      $1=="Parae_12" {sex++}
      $1 ~ /^Parae_(0[1-9]|1[0-9]|2[0-3])$/ && $1!="Parae_12" {auto++}
      END {
        print "Autosome", auto+0
        print "Sex_Ch", sex+0
      }
    ' "${ALL_CPGS}"
  } > "${GENOME_TOTALS}"

  # ----------------------------------------------------------
  # 2E) Summarize by chromosome type
  # ----------------------------------------------------------
  {
    echo -e "sample\tchromosome_type\ttotal_CpGs\tn_regions\ttotal_CpGs_inRegions\tpct_overlap\tcovered_CpGs_inSample\tpct_covered_inSample"

    awk -F $'\t' -v OFS=$'\t' -v sample="${SAMPLE_ID}" '
      BEGIN {
        order[1] = "Autosome"
        order[2] = "Sex_Ch"
      }

      FNR==NR && NR>1 {
        type = $1
        total_genome[type] = $2
        next
      }

      FNR>1 {
        type = $10
        n[type]++
        total_in_regions[type] += $11
        covered_in_sample[type] += $12
      }

      END {
        for (i=1; i<=2; i++) {
          t = order[i]
          if (t in total_genome || t in n) {
            pct_overlap = (total_genome[t] > 0 ? 100 * total_in_regions[t] / total_genome[t] : "NA")
            pct_covered = (total_in_regions[t] > 0 ? 100 * covered_in_sample[t] / total_in_regions[t] : "NA")

            print sample,
                  t,
                  total_genome[t] + 0,
                  n[t] + 0,
                  total_in_regions[t] + 0,
                  pct_overlap,
                  covered_in_sample[t] + 0,
                  pct_covered
          }
        }
      }
    ' "${GENOME_TOTALS}" "${OUT_PER_REGION}"
  } > "${OUT_SUMMARY}"

  tail -n +2 "${OUT_SUMMARY}" >> "${master_tmp}"

  echo "DONE: ${SAMPLE_ID}"
  echo "  Per-region: ${OUT_PER_REGION}"
  echo "  Summary:    ${OUT_SUMMARY}"
done

# ------------------------------------------------------------
# 3) Write combined summary across all samples
# ------------------------------------------------------------
cp "${master_tmp}" "${MASTER_SUMMARY}"

echo
echo "ALL DONE"
echo "Combined summary:"
echo "  ${MASTER_SUMMARY}"