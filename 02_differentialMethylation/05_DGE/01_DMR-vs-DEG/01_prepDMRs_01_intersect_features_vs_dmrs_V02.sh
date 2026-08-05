#!/usr/bin/env bash
# DATE:       2025-09-16
# AUTHOR:     MZF
# SCRIPT:     01_intersect_features_vs_dmrs_V02.sh
# VERSION:    02
# GOAL:       Intersect compiled features (A) vs compiled DMRs (B) using bedtools
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# USAGE:
#   bash ../scripts/01_intersect_features_vs_dmrs_V02.sh allFeatures.feature.bed dmr_all.tsv intersect_features_vs_dmrs
#   # Optional 4th arg: "loj" to use -loj (left outer join)
#   bash 04_intersect_features_vs_dmrs.sh allFeatures.feature.bed dmr_all.tsv intersect loj
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# NOTES:
# - A (features) is assumed to be BED-like without header:
#     chr  start  end  gene  feature  strand  length
# - B (DMRs) is TSV with a header; first 3 columns are chr,start,end.
# - We sort both files; we don’t use `-sorted` mode (no genome file required).
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <features_bed> <dmr_tsv> <out_name> (e.g intersect_features_vs_dmrs) [loj]"
  exit 1
fi

A_FEATURES="$1"
B_DMRS="$2"
OUT_TSV="$3"
USE_LOJ="${4:-}"

mkdir -p ${OUT_TSV}

if [[ ! -s "$A_FEATURES" ]]; then echo "ERROR: features file not found or empty: $A_FEATURES" >&2; exit 2; fi
if [[ ! -s "$B_DMRS" ]]; then echo "ERROR: DMR file not found or empty: $B_DMRS" >&2; exit 2; fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

A_SORT="$tmpdir/features.sorted.bed"
B_NOHDR="$tmpdir/dmr.noheader.tsv"
B_SORT="$tmpdir/dmr.noheader.sorted.tsv"
HDR_B="$tmpdir/dmr.header.tsv"

# 1) Capture DMR header and strip it
head -n1 "$B_DMRS" > "$HDR_B"
tail -n +2 "$B_DMRS" > "$B_NOHDR"

# 2) Sort both inputs (lexicographic chr, numeric start)
#    Use C locale to keep ordering stable for names like Parae_01 ... Parae_12
LC_ALL=C sort -k1,1 -k2,2n "$A_FEATURES" > "$A_SORT"
LC_ALL=C sort -k1,1 -k2,2n "$B_NOHDR"   > "$B_SORT"

# 3) Build a combined header for output
#    Features header is fixed; DMR header read from file and prefixed with "dmr_"
FEATURES_HDR=$'feat_chr\tfeat_start\tfeat_end\tgene\tfeature\tstrand\tlength'
DMR_HDR_PREFIXED="$(awk -F'\t' 'NR==1{
  for(i=1;i<=NF;i++){
    printf("%s%s", (i>1?OFS:""), "dmr_" $i)
  }
  printf("\n")
}' OFS='\t' "$HDR_B")"

{
  echo -e "${FEATURES_HDR}\t${DMR_HDR_PREFIXED}"

  # 4) Intersect:
  #    default = report only overlaps (-wa -wb)
  #    if 4th arg == "loj", use -loj to keep all A even if no overlap (B fields = -1)
  if [[ "${USE_LOJ}" == "loj" ]]; then
    bedtools intersect -a "$A_SORT" -b "$B_SORT" -loj
  else
    bedtools intersect -a "$A_SORT" -b "$B_SORT" -wa -wb
  fi
} > "${OUT_TSV}/${OUT_TSV}.tsv"

echo "Wrote: $OUT_TSV"
