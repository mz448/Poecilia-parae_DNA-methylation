#!/usr/bin/env bash
# DATE: 20260305
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Goal:
#   From aligned CpG report files (chr start end meth cov) and genomic features,
#   compute (per-gene, per-feature) the mean methylation and mean coverage
#     - per sample (bedtools map)
#     - averaged across the 3 samples per morph
#
# OUTPUTS generate three output tables:
#
#   A) Sample-level table with all per-sample values:
#      chr  start  end
#      cov.sample.<sample1> ... cov.sample.<sampleN>
#      meth.sample.<sample1> ... meth.sample.<sampleN>
#      gene  feature  chromosome_type
#
#   B) Morph-level table with the mean of the 3 samples per morph:
#      chr  start  end
#      cov.morph.female  cov.morph.immaculata  cov.morph.parae  cov.morph.yellow
#      meth.morph.female meth.morph.immaculata meth.morph.parae meth.morph.yellow
#      gene  feature  chromosome_type
#
#   C) Combined table with both morph means and all sample values:
#      chr  start  end
#      cov.morph.female  cov.morph.immaculata  cov.morph.parae  cov.morph.yellow
#      meth.morph.female meth.morph.immaculata meth.morph.parae meth.morph.yellow
#      cov.sample.<sample1> ... cov.sample.<sampleN>
#      meth.sample.<sample1> ... meth.sample.<sampleN>
#      gene  feature  chromosome_type
#
# IMPORTANT CHANGE (new CONS_BED format):
#   - CONS_BED uses columns 1-5 as the unique key: chr, start, end, gene, feature
#   - There are NO cons_id / dmr_ids columns
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Requirements:
#   - bedtools
#   - awk, sort, paste, wc
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail

# Inputs (edit if needed):
CONS_BED="allFeatures.feature.promoters_u400-TSS-d200.bed"  # cols 1-5 are key: chr start end gene feature
CPG_GLOB="sym.formatted.ppar*mem*.CpG_report.txt"
OUTDIR="methylation_means"
MINCOV=1
CHROM_RE='^Parae_(0[1-9]|1[0-9]|2[0-3])$'

mkdir -p "${OUTDIR}"/{mapped,morph_means,qc,final_tables}

# ------------------------------------------------------------
# 0) Sort CONS_BED for -sorted bedtools map
#    Keep all columns, but we will ONLY rely on cols 1-5 downstream.
# ------------------------------------------------------------
LC_ALL=C sort -k1,1 -k2,2n "${CONS_BED}" > "${OUTDIR}/${CONS_BED}.sorted.bed"
CONS_SORTED="${OUTDIR}/${CONS_BED}.sorted.bed"
CONS_N="$(wc -l < "${CONS_SORTED}")"

# ------------------------------------------------------------
# 1) Per-sample mapping: mean(meth) and mean(cov) per CONS_BED row
#    Input CpG reports are expected as: chr start end meth cov
# ------------------------------------------------------------
shopt -s nullglob
files_in=( ${CPG_GLOB} )
if [[ ${#files_in[@]} -eq 0 ]]; then
  echo "ERROR: No input CpG files matched: ${CPG_GLOB}"
  exit 1
fi

LC_ALL=C IFS=$'\n' files_in=( $(printf '%s\n' "${files_in[@]}" | sort) )
unset IFS

declare -a SAMPLE_NAMES=()
declare -a MAPPED_FILES=()

for INFILE in "${files_in[@]}"; do
  BASE="$(basename "${INFILE}" .CpG_report.txt)"
  OUTMAP="${OUTDIR}/mapped/${BASE}._meanMeth_meanCov.tsv"

  # remove prefix from sample name
  SAMPLE="${BASE#sym.formatted.}"

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

  n="$(wc -l < "${OUTMAP}")"
  if [[ "${n}" -ne "${CONS_N}" ]]; then
    echo "ERROR: ${OUTMAP} has ${n} rows; expected ${CONS_N}"
    exit 1
  fi

  SAMPLE_NAMES+=( "${SAMPLE}" )
  MAPPED_FILES+=( "${OUTMAP}" )
done

echo '✓ Mapped mean methylation in all regions (per sample)'
# ------------------------------------------------------------
# 2) Average across the 3 samples per morph
#    For each morph (f,i,p,y):
#      output columns:
#        chr start end gene feature
#        meth_s1 meth_s2 meth_s3 meth_mean
#        cov_s1  cov_s2  cov_s3  cov_mean
#
# NOTE:
#   bedtools map output columns are:
#     1..N  = columns from CONS_SORTED (we only use 1..5)
#     last2 = mean_meth, mean_cov
# ------------------------------------------------------------
morph_mean() {
  local morph="$1"       # f i p y
  local morph_label="$2" # female immaculata parae yellow
  local out="${OUTDIR}/morph_means/${morph}.perFeature.tsv"

  local f1="" f2="" f3=""
  local count=0

  for f in "${OUTDIR}/mapped/sym.formatted.ppar${morph}mem"*._meanMeth_meanCov.tsv; do
    [[ -e "$f" ]] || continue
    count=$((count+1))
    if   [[ $count -eq 1 ]]; then f1="$f"
    elif [[ $count -eq 2 ]]; then f2="$f"
    elif [[ $count -eq 3 ]]; then f3="$f"
    else
      echo "ERROR: Found >3 files for morph '${morph_label}'. Unexpected: $f"
      exit 1
    fi
  done

  if [[ $count -ne 3 ]]; then
    echo "ERROR: expected 3 mapped files for morph '${morph_label}', found ${count}."
    echo "  Looked for: ${OUTDIR}/mapped/sym.formatted.ppar${morph}mem*._meanMeth_meanCov.tsv"
    exit 1
  fi

  paste \
    <(cut -f1-5 "$f1") \
    <(awk -F $'\t' '{print $(NF-1)}' "$f1") \
    <(awk -F $'\t' '{print $(NF-1)}' "$f2") \
    <(awk -F $'\t' '{print $(NF-1)}' "$f3") \
    <(awk -F $'\t' '{print $NF}'     "$f1") \
    <(awk -F $'\t' '{print $NF}'     "$f2") \
    <(awk -F $'\t' '{print $NF}'     "$f3") \
  | awk -F $'\t' -v OFS=$'\t' '
      function ok(x) {
        return (x!="NA" && x!="." && x!="nan" && x!="NaN" && x!="Inf" && x!="-Inf" && x!="");
      }
      function mean3(a,b,c,   s,n) {
        s=0; n=0;
        if(ok(a)) { s+=a; n++; }
        if(ok(b)) { s+=b; n++; }
        if(ok(c)) { s+=c; n++; }
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

F_FILE="$(morph_mean f female)"
I_FILE="$(morph_mean i immaculata)"
P_FILE="$(morph_mean p parae)"
Y_FILE="$(morph_mean y yellow)"

echo '✓ Calculated DMR methylation average per morph'
# ------------------------------------------------------------
# 3A) Build sample-level table
#     chr start end
#     cov.sample.<sample1> ... cov.sample.<sampleN>
#     meth.sample.<sample1> ... meth.sample.<sampleN>
#     gene feature chromosome_type
# ------------------------------------------------------------
SAMPLE_TABLE="${OUTDIR}/final_tables/meth.sample_level.tsv"
TMPDIR_BLOCK3A="$(mktemp -d)"

# key columns
cut -f1-5 "${MAPPED_FILES[0]}" > "${TMPDIR_BLOCK3A}/key.tsv"

# sample cov columns
cov_files=()
for i in "${!MAPPED_FILES[@]}"; do
  f="${MAPPED_FILES[$i]}"
  out="${TMPDIR_BLOCK3A}/cov_${i}.tsv"
  awk -F $'\t' '{print $NF}' "$f" > "$out"
  cov_files+=( "$out" )
done

# sample meth columns
meth_files=()
for i in "${!MAPPED_FILES[@]}"; do
  f="${MAPPED_FILES[$i]}"
  out="${TMPDIR_BLOCK3A}/meth_${i}.tsv"
  awk -F $'\t' '{print $(NF-1)}' "$f" > "$out"
  meth_files+=( "$out" )
done

{
  printf "chr\tstart\tend"
  for s in "${SAMPLE_NAMES[@]}"; do
    printf "\tcov.sample.%s" "${s}"
  done
  for s in "${SAMPLE_NAMES[@]}"; do
    printf "\tmeth.sample.%s" "${s}"
  done
  printf "\tgene\tfeature\tchromosome_type\n"

  paste "${TMPDIR_BLOCK3A}/key.tsv" "${cov_files[@]}" "${meth_files[@]}" \
  | awk -F $'\t' -v OFS=$'\t' -v ns="${#SAMPLE_NAMES[@]}" '
      {
        chr=$1; start=$2; end=$3; gene=$4; feature=$5;

        if (chr=="Parae_12") chrom_type="Sex_Ch";
        else if (chr ~ /^Parae_/) chrom_type="Autosome";
        else chrom_type="Other";

        out = chr OFS start OFS end;
        for (i=6; i<6+ns; i++) out = out OFS $i;
        for (i=6+ns; i<6+(2*ns); i++) out = out OFS $i;
        out = out OFS gene OFS feature OFS chrom_type;
        print out;
      }
    '
} > "${SAMPLE_TABLE}"

rm -rf "${TMPDIR_BLOCK3A}"

echo '✓ Build sample-level table'
# ------------------------------------------------------------
# 3B) Build morph-level table
#     chr start end
#     cov.morph.female cov.morph.immaculata cov.morph.parae cov.morph.yellow
#     meth.morph.female meth.morph.immaculata meth.morph.parae meth.morph.yellow
#     gene feature chromosome_type
# ------------------------------------------------------------
MORPH_TABLE="${OUTDIR}/final_tables/meth.morph_level.tsv"

{
  echo -e "chr\tstart\tend\tcov.morph.female\tcov.morph.immaculata\tcov.morph.parae\tcov.morph.yellow\tmeth.morph.female\tmeth.morph.immaculata\tmeth.morph.parae\tmeth.morph.yellow\tgene\tfeature\tchromosome_type"

  paste \
    <(cut -f1-5  "$F_FILE") \
    <(cut -f13   "$F_FILE") \
    <(cut -f13   "$I_FILE") \
    <(cut -f13   "$P_FILE") \
    <(cut -f13   "$Y_FILE") \
    <(cut -f9    "$F_FILE") \
    <(cut -f9    "$I_FILE") \
    <(cut -f9    "$P_FILE") \
    <(cut -f9    "$Y_FILE") \
  | awk -F $'\t' -v OFS=$'\t' '
      {
        chr=$1; start=$2; end=$3; gene=$4; feature=$5;
        cov_f=$6; cov_i=$7; cov_p=$8; cov_y=$9;
        meth_f=$10; meth_i=$11; meth_p=$12; meth_y=$13;

        if (chr=="Parae_12") chrom_type="Sex_Ch";
        else if (chr ~ /^Parae_/) chrom_type="Autosome";
        else chrom_type="Other";

        print chr,start,end,
              cov_f,cov_i,cov_p,cov_y,
              meth_f,meth_i,meth_p,meth_y,
              gene,feature,chrom_type
      }
    '
} > "${MORPH_TABLE}"

echo '✓ Build morph-level table'
# ------------------------------------------------------------
# 3C) Build combined table
#     chr start end
#     cov.morph.female cov.morph.immaculata cov.morph.parae cov.morph.yellow
#     meth.morph.female meth.morph.immaculata meth.morph.parae meth.morph.yellow
#     cov.sample.<sample1> ... cov.sample.<sampleN>
#     meth.sample.<sample1> ... meth.sample.<sampleN>
#     gene feature chromosome_type
# ------------------------------------------------------------
COMBINED_TABLE="${OUTDIR}/final_tables/meth.combined.tsv"
TMPDIR_BLOCK3C="$(mktemp -d)"

# base morph block
paste \
  <(cut -f1-5  "$F_FILE") \
  <(cut -f13   "$F_FILE") \
  <(cut -f13   "$I_FILE") \
  <(cut -f13   "$P_FILE") \
  <(cut -f13   "$Y_FILE") \
  <(cut -f9    "$F_FILE") \
  <(cut -f9    "$I_FILE") \
  <(cut -f9    "$P_FILE") \
  <(cut -f9    "$Y_FILE") \
> "${TMPDIR_BLOCK3C}/morph_base.tsv"

# sample cov columns
cov_files=()
for i in "${!MAPPED_FILES[@]}"; do
  f="${MAPPED_FILES[$i]}"
  out="${TMPDIR_BLOCK3C}/cov_${i}.tsv"
  awk -F $'\t' '{print $NF}' "$f" > "$out"
  cov_files+=( "$out" )
done

# sample meth columns
meth_files=()
for i in "${!MAPPED_FILES[@]}"; do
  f="${MAPPED_FILES[$i]}"
  out="${TMPDIR_BLOCK3C}/meth_${i}.tsv"
  awk -F $'\t' '{print $(NF-1)}' "$f" > "$out"
  meth_files+=( "$out" )
done

{
  printf "chr\tstart\tend"
  printf "\tcov.morph.female\tcov.morph.immaculata\tcov.morph.parae\tcov.morph.yellow"
  printf "\tmeth.morph.female\tmeth.morph.immaculata\tmeth.morph.parae\tmeth.morph.yellow"
  for s in "${SAMPLE_NAMES[@]}"; do
    printf "\tcov.sample.%s" "${s}"
  done
  for s in "${SAMPLE_NAMES[@]}"; do
    printf "\tmeth.sample.%s" "${s}"
  done
  printf "\tgene\tfeature\tchromosome_type\n"

  paste "${TMPDIR_BLOCK3C}/morph_base.tsv" "${cov_files[@]}" "${meth_files[@]}" \
  | awk -F $'\t' -v OFS=$'\t' -v ns="${#SAMPLE_NAMES[@]}" '
      {
        chr=$1; start=$2; end=$3; gene=$4; feature=$5;

        if (chr=="Parae_12") chrom_type="Sex_Ch";
        else if (chr ~ /^Parae_/) chrom_type="Autosome";
        else chrom_type="Other";

        out = chr OFS start OFS end;

        for (i=6; i<=9; i++) out = out OFS $i;
        for (i=10; i<=13; i++) out = out OFS $i;
        for (i=14; i<14+ns; i++) out = out OFS $i;
        for (i=14+ns; i<14+(2*ns); i++) out = out OFS $i;

        out = out OFS gene OFS feature OFS chrom_type;
        print out;
      }
    '
} > "${COMBINED_TABLE}"

rm -rf "${TMPDIR_BLOCK3C}"

echo '✓ Saved combined table'
# ------------------------------------------------------------
# QC
# ------------------------------------------------------------
{
  echo "=== QC: consensus rows ==="
  echo "  ${CONS_SORTED}  (${CONS_N} rows)"

  echo "=== QC: mapped rows (sample) ==="
  ls -1 "${OUTDIR}/mapped/"*._meanMeth_meanCov.tsv 2>/dev/null | wc -l | awk '{print "  mapped sample files:", $1}'

  echo "=== QC: morph mean files ==="
  ls -1 "${OUTDIR}/morph_means/"*.perFeature.tsv 2>/dev/null | wc -l | awk '{print "  morph mean files:", $1}'

  echo "=== QC: final outputs ==="
  for f in "${SAMPLE_TABLE}" "${MORPH_TABLE}" "${COMBINED_TABLE}"; do
    echo "  ${f}"
    wc -l "${f}" | awk '{print "    rows (incl header):", $1}'
    echo "    header:"
    head -n 1 "${f}"
    echo "    first data row:"
    sed -n '2p' "${f}"
  done
} > "${OUTDIR}/qc/qc.summary.txt"

echo "DONE"
echo "Sample-level: ${SAMPLE_TABLE}"
echo "Morph-level:  ${MORPH_TABLE}"
echo "Combined:     ${COMBINED_TABLE}"
echo "QC:           ${OUTDIR}/qc/qc.summary.txt"