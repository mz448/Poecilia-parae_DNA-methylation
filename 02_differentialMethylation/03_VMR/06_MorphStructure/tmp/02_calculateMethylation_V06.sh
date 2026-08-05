#!/usr/bin/env bash
# DATE: 20260220
# AUTHOR: MZF
# Script: 04_calculateMeth_V06.sh
# Version 06
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   From symetric CpG report files and a VMR file, compute 
#   per-consensus-region mean methylation and mean coverage:
#
#     1) per sample (bedtools map)
#     2) averaged across the 3 samples per morph
#
#   Then generate three output tables:
#
#     A) Sample-level table with all per-sample values:
#        chr  start  end
#        cov.score_<sample1> ... cov.score_<sampleN>
#        meth_<sample1>      ... meth_<sampleN>
#        cons_id  n_dmrs  dmr_ids  chromosome_type  n_meth_nonNA_samples
#
#     B) Morph-level table with the mean of the 3 samples per morph:
#        chr  start  end
#        cov.score_female  cov.score_immaculata  cov.score_parae  cov.score_yellow
#        meth_female       meth_immaculata       meth_parae       meth_yellow
#        cons_id  n_dmrs  dmr_ids  chromosome_type  n_meth_nonNA
#
#     C) Combined table with both morph means and all sample values:
#        chr  start  end
#        cov.score_female  cov.score_immaculata  cov.score_parae  cov.score_yellow
#        cov.score_<sample1> ... cov.score_<sampleN>
#        meth_female       meth_immaculata       meth_parae       meth_yellow
#        meth_<sample1>    ... meth_<sampleN>
#        cons_id  n_dmrs  dmr_ids  chromosome_type  n_meth_nonNA  n_meth_nonNA_samples
#
#   This allows downstream analyses to use:
#     - morph-consensus methylation values
#     - per-sample methylation values
#     - a single integrated table containing both
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Requirements:
#   - bedtools
#   - awk, sort, paste, wc
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# V05: Adds the initial formating of the symetric files
# V06: Fixes VMR naming (replaces CONS#### and DMR in column and file names)
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail

# Inputs (edit if needed):
# -----------------------------
# Select one condition
# -----------------------------

CONDITION="real_morph"
# CONDITION="shuffled_morph"
# CONDITION="real_sex"
# CONDITION="shuffled_sex"


CONS_BED="vmr_real.bed"   # chr start end cons_id dmr_ids
SOURCE_DIR="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/01_DMC/01_DNMTools/data" # source of CpG counts
SAMPLES="Sym.forSym.ppar*mem*.CpG_report.txt.meth" # Modify if source file naming changes

CPG_GLOB="sym.formatted.ppar*mem*.CpG_report.txt"
OUTDIR="vmr_methylation"
MINCOV=1                            # ignore CpGs with cov < MINCOV
CHROM_RE='^Parae_(0[1-9]|1[0-9]|2[0-3])$'

mkdir -p "${OUTDIR}"/{mapped,morph_means,qc}




# ------------------------------------------------------------
# Prepare Samples from sym files
# Input:  chr start strand context meth% cov
# Output: chr start end    meth%   cov
# ------------------------------------------------------------

# Copy sample CpG report files into the current project folder
cp "${SOURCE_DIR}"/Sym.forSym.ppar*mem*.CpG_report.txt.meth .

# Reformat and rename each sample
for INPUT in Sym.forSym.ppar*mem*.CpG_report.txt.meth; do

    # Example:
    # Sym.forSym.pparfmem001.CpG_report.txt.meth
    # ->
    # sym.formatted.pparfmem001.CpG_report.txt

    OUTPUT="${INPUT/Sym.forSym./sym.formatted.}"
    OUTPUT="${OUTPUT%.meth}"

    awk 'BEGIN {OFS="\t"} {
        print $1, $2, $2+1, $5, $6
    }' "$INPUT" > "$OUTPUT"

done


# ------------------------------------------------------------
# 0) Sort consensus BED for -sorted bedtools map
# ------------------------------------------------------------
LC_ALL=C sort -k1,1 -k2,2n "${CONS_BED}" > "${OUTDIR}/vmr.sorted.bed"
CONS_SORTED="${OUTDIR}/vmr.sorted.bed"
CONS_N="$(wc -l < "${CONS_SORTED}")"

# ------------------------------------------------------------
# 1) Per-sample mapping: mean(meth) and mean(cov) per consensus DMR
#    (no intermediate 4-col files needed; stream filter -> bedtools map)
# ------------------------------------------------------------
shopt -s nullglob
files_in=( ${CPG_GLOB} )
if [[ ${#files_in[@]} -eq 0 ]]; then
  echo "ERROR: No input CpG files matched: ${CPG_GLOB}"
  exit 1
fi

for INFILE in "${files_in[@]}"; do
  BASE="$(basename "${INFILE}" .CpG_report.txt)"
  OUTMAP="${OUTDIR}/mapped/${BASE}.vmr_meanMeth_meanCov.tsv"

  awk -F $'\t' -v OFS=$'\t' -v re="${CHROM_RE}" -v mincov="${MINCOV}" '
    ($1 ~ re) && ($5+0 >= mincov) { print $1,$2,$3,$4,$5 }
  ' "${INFILE}" \
  | bedtools map \
      -a "${CONS_SORTED}" \
      -b - \
      -c 4,5 \
      -o mean,mean \
      -null NA \
      -sorted \
  > "${OUTMAP}"

  # QC: ensure same number of rows as consensus
  n="$(wc -l < "${OUTMAP}")"
  if [[ "${n}" -ne "${CONS_N}" ]]; then
    echo "ERROR: ${OUTMAP} has ${n} rows; expected ${CONS_N}"
    exit 1
  fi
done

# ------------------------------------------------------------
# 1b) Build per-consensus table with all sample values
#     output columns:
#       chr start end
#       cov.score_<sample1> ... cov.score_<sampleN>
#       meth_<sample1> ... meth_<sampleN>
#       cons_id n_dmrs dmr_ids chromosome_type n_meth_nonNA_samples
# ------------------------------------------------------------
SAMPLE_FINAL="${OUTDIR}/vmr.sample_observed-methylation.tsv"

mapped_files=()
while IFS= read -r f; do
  mapped_files+=( "$f" )
done < <(printf '%s\n' "${OUTDIR}"/mapped/sym.formatted.ppar*mem*.vmr_meanMeth_meanCov.tsv | LC_ALL=C sort)

if [[ ${#mapped_files[@]} -eq 0 ]]; then
  echo "ERROR: no mapped sample files found in ${OUTDIR}/mapped/"
  exit 1
fi

tmp_sample_dir="$(mktemp -d)"
sample_ids=()

for f in "${mapped_files[@]}"; do
  sid="$(basename "$f" .vmr_meanMeth_meanCov.tsv)"
  sid="${sid#sym.formatted.}"
  sample_ids+=( "$sid" )
done

meta_file="${tmp_sample_dir}/meta.tsv"
cut -f1-5 "${mapped_files[0]}" > "${meta_file}"

cov_parts=()
meth_parts=()

i=0
for f in "${mapped_files[@]}"; do
  i=$((i+1))
  cov_file="${tmp_sample_dir}/cov_${i}.tsv"
  meth_file="${tmp_sample_dir}/meth_${i}.tsv"

  # mapped file columns:
  # 1 chr, 2 start, 3 end, 4 cons_id, 5 dmr_ids, 6 mean_meth, 7 mean_cov
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
  printf "\tcons_id\tn_dmrs\tdmr_ids\tchromosome_type\tn_meth_nonNA_samples\n"

  paste "${meta_file}" "${cov_parts[@]}" "${meth_parts[@]}" \
  | awk -F $'\t' -v OFS=$'\t' -v ns="${#sample_ids[@]}" '
      {
        chr     = $1
        start   = $2
        end     = $3
        cons_id = $4
        dmr_ids = $5

        n_dmrs = split(dmr_ids, a, ",")
        chrom_type = (chr=="Parae_12" ? "Sex_Ch" : "Autosome")

        n_nonNA = 0
        for (i = 6 + ns; i <= 5 + (2 * ns); i++) {
          if ($i != "NA" && $i != "." && $i != "nan" && $i != "NaN" && $i != "Inf" && $i != "-Inf") {
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

        printf "\t%s\t%s\t%s\t%s\t%s\n", cons_id, n_dmrs, dmr_ids, chrom_type, n_nonNA
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
#    For each morph (f,i,p,y):
#      output columns:
#        chr start end cons_id dmr_ids  meth_s1 meth_s2 meth_s3  meth_mean  cov_s1 cov_s2 cov_s3  cov_mean
# ------------------------------------------------------------
morph_mean() {
  local morph="$1"  # f i p y
  local out="${OUTDIR}/morph_means/${morph}.perVMR.tsv"

  # Collect files for this morph (expect exactly 3)
  local f1="" f2="" f3=""
  local count=0

  for f in "${OUTDIR}/mapped/sym.formatted.ppar${morph}mem"*.vmr_meanMeth_meanCov.tsv; do
    [[ -e "$f" ]] || continue
    count=$((count+1))
    if   [[ $count -eq 1 ]]; then f1="$f"
    elif [[ $count -eq 2 ]]; then f2="$f"
    elif [[ $count -eq 3 ]]; then f3="$f"
    else
      echo "ERROR: Found >3 files for morph '${morph}'. Unexpected: $f"
      exit 1
    fi
  done

  if [[ $count -ne 3 ]]; then
    echo "ERROR: expected 3 mapped files for morph '${morph}', found ${count}."
    echo "  Looked for: ${OUTDIR}/mapped/sym.formatted.ppar${morph}mem*.vmr_meanMeth_meanCov.tsv"
    exit 1
  fi

  # File columns from bedtools map:
  #   1 chr, 2 start, 3 end, 4 cons_id, 5 dmr_ids, 6 mean_meth, 7 mean_cov
  paste \
    <(cut -f1-5 "$f1") \
    <(cut -f6   "$f1") <(cut -f6 "$f2") <(cut -f6 "$f3") \
    <(cut -f7   "$f1") <(cut -f7 "$f2") <(cut -f7 "$f3") \
  | awk -F $'\t' -v OFS=$'\t' '
      function mean3(a,b,c,   s,n) {
        s=0; n=0;
        if(a!="NA" && a!="." && a!="nan" && a!="NaN" && a!="Inf" && a!="-Inf"){ s+=a; n++; }
        if(b!="NA" && b!="." && b!="nan" && b!="NaN" && b!="Inf" && b!="-Inf"){ s+=b; n++; }
        if(c!="NA" && c!="." && c!="nan" && c!="NaN" && c!="Inf" && c!="-Inf"){ s+=c; n++; }
        return (n>0 ? s/n : "NA");
      }
      {
        meth_mean = mean3($6,$7,$8);
        cov_mean  = mean3($9,$10,$11);
        print $1,$2,$3,$4,$5, $6,$7,$8,meth_mean, $9,$10,$11,cov_mean
      }
    ' > "${out}"

  echo "${out}"
}

F_FILE="$(morph_mean f)"
I_FILE="$(morph_mean i)"
P_FILE="$(morph_mean p)"
Y_FILE="$(morph_mean y)"

# ------------------------------------------------------------
# 3) Build final master table matching your previous format
# ------------------------------------------------------------
FINAL="${OUTDIR}/vmr.morph_observed.tsv"

{
  echo -e "chr\tstart\tend\tcov.score_female\tcov.score_immaculata\tcov.score_parae\tcov.score_yellow\tmeth_female\tmeth_immaculata\tmeth_parae\tmeth_yellow\tcons_id\tn_dmrs\tdmr_ids\tchromosome_type\tn_meth_nonNA"

  paste \
    <(cut -f1-5 "$F_FILE") \
    <(cut -f13 "$F_FILE") \
    <(cut -f13 "$I_FILE") \
    <(cut -f13 "$P_FILE") \
    <(cut -f13 "$Y_FILE") \
    <(cut -f9  "$F_FILE") \
    <(cut -f9  "$I_FILE") \
    <(cut -f9  "$P_FILE") \
    <(cut -f9  "$Y_FILE") \
  | awk -F $'\t' -v OFS=$'\t' '
      {
        chr=$1; start=$2; end=$3; cons_id=$4; dmr_ids=$5;
        cov_f=$6; cov_i=$7; cov_p=$8; cov_y=$9;
        meth_f=$10; meth_i=$11; meth_p=$12; meth_y=$13;

        # n_dmrs = number of comma-separated IDs
        n_dmrs = split(dmr_ids, a, ",");

        # chromosome_type
        chrom_type = (chr=="Parae_12" ? "Sex_Ch" : "Autosome");

        # n_meth_nonNA
        n_nonNA=0;
        if(meth_f!="NA") n_nonNA++;
        if(meth_i!="NA") n_nonNA++;
        if(meth_p!="NA") n_nonNA++;
        if(meth_y!="NA") n_nonNA++;

        print chr,start,end,
              cov_f,cov_i,cov_p,cov_y,
              meth_f,meth_i,meth_p,meth_y,
              cons_id,n_dmrs,dmr_ids,chrom_type,n_nonNA
      }
    '
} > "${FINAL}"

# ------------------------------------------------------------
# 3b) Build combined table with morph means + all sample values
#     output columns:
#       chr start end
#       cov.score_female cov.score_immaculata cov.score_parae cov.score_yellow
#       cov.score_<sample1> ... cov.score_<sampleN>
#       meth_female meth_immaculata meth_parae meth_yellow
#       meth_<sample1> ... meth_<sampleN>
#       cons_id n_dmrs dmr_ids chromosome_type n_meth_nonNA n_meth_nonNA_samples
# ------------------------------------------------------------
BIG_FINAL="${OUTDIR}/vmr.morph_and_sample_observed.tsv"

ns="${#sample_ids[@]}"
sample_cov_start=4
sample_cov_end=$((3 + ns))
sample_meth_start=$((4 + ns))
sample_meth_end=$((3 + 2 * ns))
sample_nnonna_col=$((8 + 2 * ns))

{
  printf "chr\tstart\tend"
  printf "\tcov.score_female\tcov.score_immaculata\tcov.score_parae\tcov.score_yellow"
  for sid in "${sample_ids[@]}"; do
    printf "\tcov.score_%s" "${sid}"
  done
  printf "\tmeth_female\tmeth_immaculata\tmeth_parae\tmeth_yellow"
  for sid in "${sample_ids[@]}"; do
    printf "\tmeth_%s" "${sid}"
  done
  printf "\tcons_id\tn_dmrs\tdmr_ids\tchromosome_type\tn_meth_nonNA\tn_meth_nonNA_samples\n"

  paste \
    <(tail -n +2 "${FINAL}"        | cut -f1-3) \
    <(tail -n +2 "${FINAL}"        | cut -f4-7) \
    <(tail -n +2 "${SAMPLE_FINAL}" | cut -f"${sample_cov_start}-${sample_cov_end}") \
    <(tail -n +2 "${FINAL}"        | cut -f8-11) \
    <(tail -n +2 "${SAMPLE_FINAL}" | cut -f"${sample_meth_start}-${sample_meth_end}") \
    <(tail -n +2 "${FINAL}"        | cut -f12-16) \
    <(tail -n +2 "${SAMPLE_FINAL}" | cut -f"${sample_nnonna_col}")
} > "${BIG_FINAL}"

echo "=== QC: combined morph + sample output ==="
echo "  ${BIG_FINAL}"
wc -l "${BIG_FINAL}" | awk '{print "  rows:", $1}'
echo "  header:"
head -n 1 "${BIG_FINAL}"
echo "  first data row:"
sed -n '2p' "${BIG_FINAL}"

# ------------------------------------------------------------
# QC
# ------------------------------------------------------------
echo "=== QC: consensus rows ==="
echo "  ${CONS_SORTED}  (${CONS_N} rows)"

echo "=== QC: mapped rows (sample) ==="
ls -1 "${OUTDIR}/mapped/"*.vmr_meanMeth_meanCov.tsv | wc -l | awk '{print "  mapped sample files:", $1}'

echo "=== QC: final output ==="
echo "  ${FINAL}"
wc -l "${FINAL}" | awk '{print "  rows:", $1}'
echo "  header:"
head -n 1 "${FINAL}"
echo "  first data row:"
sed -n '2p' "${FINAL}"