#!/usr/bin/env bash
#SBATCH --job-name=methmap
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

# DATE: 20260220
# AUTHOR: MZF
# Script: 02_calculateMethylation_V08.sh
# Version 08
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Compute per-VMR mean methylation and mean coverageFrom existing formatted
#   symmetric CpG report files and a condition-specific VMR file:
#
#     1) per sample using bedtools map
#     2) averaged across the 3 samples per morph
#
#   Then generate three output tables:
#
#     A) Sample-level table with all per-sample values:
#        chr  start  end
#        cov.score_<sample1> ... cov.score_<sampleN>
#        meth_<sample1>      ... meth_<sampleN>
#        vmr_id  n_dmrs  dmr_ids  chromosome_type  n_meth_nonNA_samples
#
#     B) Morph-level table with the mean of the 3 samples per morph:
#        chr  start  end
#        cov.score_female  cov.score_immaculata  cov.score_parae  cov.score_yellow
#        meth_female       meth_immaculata       meth_parae       meth_yellow
#        vmr_id  n_dmrs  dmr_ids  chromosome_type  n_meth_nonNA
#
#     C) Combined table with both morph means and all sample values:
#        chr  start  end
#        cov.score_female  cov.score_immaculata  cov.score_parae  cov.score_yellow
#        cov.score_<sample1> ... cov.score_<sampleN>
#        meth_female       meth_immaculata       meth_parae       meth_yellow
#        meth_<sample1>    ... meth_<sampleN>
#        vmr_id  n_dmrs  dmr_ids  chromosome_type
#        n_meth_nonNA  n_meth_nonNA_samples
#
#   The script must be run from the VMR data directory containing:
#
#     sym.formatted.ppar*mem*.CpG_report.txt
#     dataset_Structure/vmr_<condition>.bed
#
# Requirements:
#   - bedtools
#   - awk, sort, paste, wc
#
# Usage:
#   sbatch ../scripts/06_MorphStructure/02_calculateMethylation_V08.sh real_morph
#   sbatch ../scripts/06_MorphStructure/02_calculateMethylation_V08.sh shuffled_morph
#   sbatch ../scripts/06_MorphStructure/02_calculateMethylation_V08.sh real_sex
#   sbatch ../scripts/06_MorphStructure/02_calculateMethylation_V08.sh shuffled_sex
#   sbatch ../scripts/06_MorphStructure/02_calculateMethylation_V08.sh shuffled_all
#
# V05: Adds the initial formatting of the symmetric files.
# V06: Fixes VMR naming, replacing CONS#### and DMR terminology.
# V07: Uses the existing sym.formatted CpG files from the data directory,
#      reads condition-specific VMR files from dataset_Structure, and saves
#      all outputs in a condition-specific folder inside dataset_Structure.
# V08: Fix awk problem, labels master outputs with condition and reuses samples 
#      if already calculated.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail


# -----------------------------
# Inputs
# -----------------------------

CONDITION="${1:-real_morph}"

VMR_BED="dataset_Structure/vmr_${CONDITION}.bed"

CPG_GLOB="sym.formatted.ppar*mem*.CpG_report.txt"

OUTDIR="dataset_Structure/${CONDITION}_vmr_methylation"

MINCOV=1
CHROM_RE='^Parae_(0[1-9]|1[0-9]|2[0-3])$'


# -----------------------------
# Check inputs and create outputs
# -----------------------------

if [[ ! -f "${VMR_BED}" ]]; then
  echo "ERROR: condition-specific VMR file not found:"
  echo "  ${VMR_BED}"
  exit 1
fi

mkdir -p "${OUTDIR}"/{mapped,morph_means,qc}


# ------------------------------------------------------------
# 0) Sort VMR BED for bedtools map -sorted
# ------------------------------------------------------------

LC_ALL=C sort -k1,1 -k2,2n "${VMR_BED}" \
  > "${OUTDIR}/vmr.${CONDITION}.sorted.bed"

VMR_SORTED="${OUTDIR}/vmr.${CONDITION}.sorted.bed"
VMR_N="$(wc -l < "${VMR_SORTED}")"


# ------------------------------------------------------------
# 1) Per-sample mapping
#
# Input formatted CpG files:
#   chr  start  end  meth  cov
#
# VMR BED:
#   chr  start  end  vmr_id  dmr_ids
#
# Output from bedtools map:
#   chr  start  end  vmr_id  dmr_ids  mean_meth  mean_cov
# ------------------------------------------------------------

shopt -s nullglob

files_in=( ${CPG_GLOB} )

if [[ ${#files_in[@]} -eq 0 ]]; then
  echo "ERROR: no existing formatted CpG files matched:"
  echo "  ${CPG_GLOB}"
  exit 1
fi

for INFILE in "${files_in[@]}"; do

  BASE="$(basename "${INFILE}" .CpG_report.txt)"

  OUTMAP="${OUTDIR}/mapped/${BASE}.vmr_meanMeth_meanCov.tsv"
  
  # Reuse a completed mapped file when it has the expected number of VMR rows.
  if [[ -s "${OUTMAP}" ]] &&
     [[ "$(wc -l < "${OUTMAP}")" -eq "${VMR_N}" ]]; then
    echo "Reusing completed mapped file: ${OUTMAP}"
    continue
  fi
  
  awk -F $'\t' \
    -v OFS=$'\t' \
    -v re="${CHROM_RE}" \
    -v mincov="${MINCOV}" '
      ($1 ~ re) && ($5 + 0 >= mincov) {
        print $1, $2, $3, $4, $5
      }
    ' "${INFILE}" \
  | bedtools map \
      -a "${VMR_SORTED}" \
      -b - \
      -c 4,5 \
      -o mean,mean \
      -null NA \
      -sorted \
  > "${OUTMAP}"

  # QC: every mapped file must contain one row per VMR.
  n="$(wc -l < "${OUTMAP}")"

  if [[ "${n}" -ne "${VMR_N}" ]]; then
    echo "ERROR: ${OUTMAP} has ${n} rows; expected ${VMR_N}"
    exit 1
  fi

done


# ------------------------------------------------------------
# 1b) Build per-VMR table with all sample values
#
# Output columns:
#   chr start end
#   cov.score_<sample1> ... cov.score_<sampleN>
#   meth_<sample1> ... meth_<sampleN>
#   vmr_id n_dmrs dmr_ids chromosome_type n_meth_nonNA_samples
# ------------------------------------------------------------

SAMPLE_FINAL="${OUTDIR}/vmr.${CONDITION}.sample_observed-methylation.tsv"

mapped_files=()

while IFS= read -r f; do
  mapped_files+=( "${f}" )
done < <(
  printf '%s\n' \
    "${OUTDIR}"/mapped/sym.formatted.ppar*mem*.vmr_meanMeth_meanCov.tsv \
  | LC_ALL=C sort
)

if [[ ${#mapped_files[@]} -eq 0 ]]; then
  echo "ERROR: no mapped sample files found in:"
  echo "  ${OUTDIR}/mapped/"
  exit 1
fi

tmp_sample_dir="$(mktemp -d)"
sample_ids=()

for f in "${mapped_files[@]}"; do

  sid="$(basename "${f}" .vmr_meanMeth_meanCov.tsv)"
  sid="${sid#sym.formatted.}"

  sample_ids+=( "${sid}" )

done

meta_file="${tmp_sample_dir}/meta.tsv"

cut -f1-5 "${mapped_files[0]}" > "${meta_file}"

cov_parts=()
meth_parts=()

i=0

for f in "${mapped_files[@]}"; do

  i=$((i + 1))

  cov_file="${tmp_sample_dir}/cov_${i}.tsv"
  meth_file="${tmp_sample_dir}/meth_${i}.tsv"

  # Mapped-file columns:
  #   1 chr
  #   2 start
  #   3 end
  #   4 vmr_id
  #   5 dmr_ids
  #   6 mean_meth
  #   7 mean_cov

  cut -f7 "${f}" > "${cov_file}"
  cut -f6 "${f}" > "${meth_file}"

  cov_parts+=( "${cov_file}" )
  meth_parts+=( "${meth_file}" )

done

{
  printf "chr\tstart\tend"

  for sid in "${sample_ids[@]}"; do
    printf "\tcov.score_%s" "${sid}"
  done

  for sid in "${sample_ids[@]}"; do
    printf "\tmeth_%s" "${sid}"
  done

  printf "\tvmr_id\tn_dmrs\tdmr_ids\tchromosome_type"
  printf "\tn_meth_nonNA_samples\n"

  paste \
    "${meta_file}" \
    "${cov_parts[@]}" \
    "${meth_parts[@]}" \
  | awk -F $'\t' \
      -v OFS=$'\t' \
      -v ns="${#sample_ids[@]}" '
      {
        chr     = $1
        start   = $2
        end     = $3
        vmr_id  = $4
        dmr_ids = $5

        n_dmrs = split(dmr_ids, a, ",")

        chrom_type = (chr == "Parae_12" ? "Sex_Ch" : "Autosome")

        n_nonNA = 0

        for (i = 6 + ns; i <= 5 + (2 * ns); i++) {
          if ($i != "NA" && $i != "." && $i != "nan" &&
              $i != "NaN" && $i != "Inf" && $i != "-Inf") {
            n_nonNA++
          }
        }

        printf "%s\t%s\t%s", chr, start, end

        for (i = 6; i <= 5 + ns; i++) {
          printf "\t%s", $i
        }

        for (i = 6 + ns; i <= 5 + (2 * ns); i++) {
          printf "\t%s", $i
        }

        printf "\t%s\t%s\t%s\t%s\t%s\n",
          vmr_id,
          n_dmrs,
          dmr_ids,
          chrom_type,
          n_nonNA
      }
    '

} > "${SAMPLE_FINAL}"

rm -rf "${tmp_sample_dir}"

echo "=== QC: sample-level output ==="
echo "  ${SAMPLE_FINAL}"
wc -l "${SAMPLE_FINAL}" | awk '{print "  rows:", $1}'
echo "  header:"
head -n 1 "${SAMPLE_FINAL}"
echo "  first data row:"
sed -n '2p' "${SAMPLE_FINAL}"


# ------------------------------------------------------------
# 2) Average across the 3 samples per morph
#
# For each morph:
#   f = female
#   i = immaculata
#   p = parae
#   y = yellow
#
# Output columns:
#   chr start end vmr_id dmr_ids
#   meth_s1 meth_s2 meth_s3 meth_mean
#   cov_s1 cov_s2 cov_s3 cov_mean
# ------------------------------------------------------------

morph_mean() {

  local morph="$1"
  local out="${OUTDIR}/morph_means/${morph}.perVMR.tsv"

  local f1=""
  local f2=""
  local f3=""
  local count=0
  local n=0

  for f in \
    "${OUTDIR}/mapped/sym.formatted.ppar${morph}mem"*.vmr_meanMeth_meanCov.tsv
  do
    [[ -e "${f}" ]] || continue

    count=$((count + 1))

    if [[ ${count} -eq 1 ]]; then
      f1="${f}"
    elif [[ ${count} -eq 2 ]]; then
      f2="${f}"
    elif [[ ${count} -eq 3 ]]; then
      f3="${f}"
    else
      echo "ERROR: found more than 3 files for morph '${morph}'."
      echo "Unexpected file: ${f}"
      return 1
    fi
  done

  if [[ ${count} -ne 3 ]]; then
    echo "ERROR: expected 3 mapped files for morph '${morph}', found ${count}."
    echo "Looked for:"
    echo "  ${OUTDIR}/mapped/sym.formatted.ppar${morph}mem*.vmr_meanMeth_meanCov.tsv"
    return 1
  fi

    if ! paste \
    <(cut -f1-5 "${f1}") \
    <(cut -f6 "${f1}") \
    <(cut -f6 "${f2}") \
    <(cut -f6 "${f3}") \
    <(cut -f7 "${f1}") \
    <(cut -f7 "${f2}") \
    <(cut -f7 "${f3}") \
  | awk -F $'\t' -v OFS=$'\t' '
      function valid(x) {
        return x != "NA" && x != "." && x != "nan" && x != "NaN" && x != "Inf" && x != "-Inf"
      }

      function mean3(a, b, c, s, n) {
        s = 0
        n = 0

        if (valid(a)) {
          s = s + a
          n = n + 1
        }

        if (valid(b)) {
          s = s + b
          n = n + 1
        }

        if (valid(c)) {
          s = s + c
          n = n + 1
        }

        if (n > 0) {
          return s / n
        }

        return "NA"
      }

      {
        meth_mean = mean3($6, $7, $8)
        cov_mean = mean3($9, $10, $11)

        print $1, $2, $3, $4, $5, $6, $7, $8, meth_mean, $9, $10, $11, cov_mean
      }
    ' > "${out}"
  then
    echo "ERROR: failed to calculate morph means for '${morph}'."
    rm -f "${out}"
    return 1
  fi

  n="$(wc -l < "${out}")"

  if [[ "${n}" -ne "${VMR_N}" ]]; then
    echo "ERROR: ${out} has ${n} rows; expected ${VMR_N}."
    rm -f "${out}"
    return 1
  fi

  echo "Created: ${out} (${n} rows)"
}


# Generate the four morph-level intermediate tables.
morph_mean f
morph_mean i
morph_mean p
morph_mean y

# Explicitly define their filenames.
F_FILE="${OUTDIR}/morph_means/f.perVMR.tsv"
I_FILE="${OUTDIR}/morph_means/i.perVMR.tsv"
P_FILE="${OUTDIR}/morph_means/p.perVMR.tsv"
Y_FILE="${OUTDIR}/morph_means/y.perVMR.tsv"


# ------------------------------------------------------------
# 3) Build morph-level master table
# ------------------------------------------------------------

FINAL="${OUTDIR}/vmr.${CONDITION}.morph_observed.tsv"

{
  echo -e \
    "chr\tstart\tend"\
"\tcov.score_female\tcov.score_immaculata\tcov.score_parae\tcov.score_yellow"\
"\tmeth_female\tmeth_immaculata\tmeth_parae\tmeth_yellow"\
"\tvmr_id\tn_dmrs\tdmr_ids\tchromosome_type\tn_meth_nonNA"

  paste \
    <(cut -f1-5 "${F_FILE}") \
    <(cut -f13 "${F_FILE}") \
    <(cut -f13 "${I_FILE}") \
    <(cut -f13 "${P_FILE}") \
    <(cut -f13 "${Y_FILE}") \
    <(cut -f9 "${F_FILE}") \
    <(cut -f9 "${I_FILE}") \
    <(cut -f9 "${P_FILE}") \
    <(cut -f9 "${Y_FILE}") \
  | awk -F $'\t' -v OFS=$'\t' '
      {
        chr     = $1
        start   = $2
        end     = $3
        vmr_id  = $4
        dmr_ids = $5

        cov_f = $6
        cov_i = $7
        cov_p = $8
        cov_y = $9

        meth_f = $10
        meth_i = $11
        meth_p = $12
        meth_y = $13

        n_dmrs = split(dmr_ids, a, ",")

        chrom_type = (chr == "Parae_12" ? "Sex_Ch" : "Autosome")

        n_nonNA = 0

        if (meth_f != "NA") n_nonNA++
        if (meth_i != "NA") n_nonNA++
        if (meth_p != "NA") n_nonNA++
        if (meth_y != "NA") n_nonNA++

        print \
          chr, start, end, \
          cov_f, cov_i, cov_p, cov_y, \
          meth_f, meth_i, meth_p, meth_y, \
          vmr_id, n_dmrs, dmr_ids, chrom_type, n_nonNA
      }
    '

} > "${FINAL}"


# ------------------------------------------------------------
# 3b) Build combined morph-mean and sample-level table
# ------------------------------------------------------------

BIG_FINAL="${OUTDIR}/vmr.${CONDITION}.morph_and_sample_observed.tsv"

ns="${#sample_ids[@]}"

sample_cov_start=4
sample_cov_end=$((3 + ns))

sample_meth_start=$((4 + ns))
sample_meth_end=$((3 + (2 * ns)))

sample_nnonna_col=$((8 + (2 * ns)))

{
  printf "chr\tstart\tend"

  printf "\tcov.score_female"
  printf "\tcov.score_immaculata"
  printf "\tcov.score_parae"
  printf "\tcov.score_yellow"

  for sid in "${sample_ids[@]}"; do
    printf "\tcov.score_%s" "${sid}"
  done

  printf "\tmeth_female"
  printf "\tmeth_immaculata"
  printf "\tmeth_parae"
  printf "\tmeth_yellow"

  for sid in "${sample_ids[@]}"; do
    printf "\tmeth_%s" "${sid}"
  done

  printf "\tvmr_id"
  printf "\tn_dmrs"
  printf "\tdmr_ids"
  printf "\tchromosome_type"
  printf "\tn_meth_nonNA"
  printf "\tn_meth_nonNA_samples\n"

  paste \
    <(tail -n +2 "${FINAL}"        | cut -f1-3) \
    <(tail -n +2 "${FINAL}"        | cut -f4-7) \
    <(tail -n +2 "${SAMPLE_FINAL}" | cut -f"${sample_cov_start}-${sample_cov_end}") \
    <(tail -n +2 "${FINAL}"        | cut -f8-11) \
    <(tail -n +2 "${SAMPLE_FINAL}" | cut -f"${sample_meth_start}-${sample_meth_end}") \
    <(tail -n +2 "${FINAL}"        | cut -f12-16) \
    <(tail -n +2 "${SAMPLE_FINAL}" | cut -f"${sample_nnonna_col}")

} > "${BIG_FINAL}"


# ------------------------------------------------------------
# QC
# ------------------------------------------------------------

echo "=== QC: analysis condition ==="
echo "  ${CONDITION}"

echo "=== QC: input VMR file ==="
echo "  ${VMR_BED}"

echo "=== QC: sorted VMR rows ==="
echo "  ${VMR_SORTED} (${VMR_N} rows)"

echo "=== QC: formatted CpG input files ==="
echo "  ${#files_in[@]} files matched"

echo "=== QC: mapped sample files ==="
ls -1 "${OUTDIR}/mapped/"*.vmr_meanMeth_meanCov.tsv \
  | wc -l \
  | awk '{print "  mapped sample files:", $1}'

echo "=== QC: morph-level output ==="
echo "  ${FINAL}"
wc -l "${FINAL}" | awk '{print "  rows:", $1}'
echo "  header:"
head -n 1 "${FINAL}"
echo "  first data row:"
sed -n '2p' "${FINAL}"

echo "=== QC: combined morph + sample output ==="
echo "  ${BIG_FINAL}"
wc -l "${BIG_FINAL}" | awk '{print "  rows:", $1}'
echo "  header:"
head -n 1 "${BIG_FINAL}"
echo "  first data row:"
sed -n '2p' "${BIG_FINAL}"

echo "=== DONE ==="
