#!/usr/bin/env bash
# DATE: 20260315
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   For each consensus DMR, identify the closest mRNA feature and report
#   the signed distance relative to gene strand:
#     - upstream of the gene start = negative
#     - overlap with the gene body = 0
#     - downstream of the gene body = positive
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# BIOLOGICAL CONTEXT:
#   This annotates strict consensus DMRs by the closest gene model using mRNA
#   intervals from the genome feature file. Distances are gene-centric and
#   strand-aware:
#     + genes: upstream = lower coordinates, downstream = higher coordinates
#     - genes: upstream = higher coordinates, downstream = lower coordinates
#
#   Because full mRNA intervals are used:
#     - any cDMR overlapping the gene body gets distance = 0
#     - non-overlapping upstream cDMRs get negative distances
#     - non-overlapping downstream cDMRs get positive distances
#
# INPUTS:
#   1) allFeatures.feature.promoters_u400-TSS-d200.bed
#      Expected columns:
#        1 chr
#        2 start
#        3 end
#        4 gene
#        5 feature
#        6 strand
#        7 length
#      Only rows with feature == "mRNA" are used.
#
#   2) cDMRs_diff0.26.bed
#      BED6:
#        chr start end name score strand
#      where:
#        name = CONSxxxxxxx|Cluster_X
#
# OUTPUT:
#   strict-consensus-DMRs_diff0.26.closest_mRNA.tsv
#     dmr_chr, dmr_start, dmr_end, cons_id, cluster, dmr_score,
#     gene_chr, gene_start, gene_end, gene, gene_strand,
#     dist_signed, abs_dist_bp, chromosome_type
#
# METHOD SUMMARY:
#   1) Extract mRNA intervals from the genome feature file
#   2) Build a clean strict scDMR BED
#   3) Use bedtools closest -D a with genes as A and scDMRs as B
#   4) Convert gene->DMR matches into a DMR-centered candidate table
#   5) Retain the closest gene per DMR
#   6) Add chromosome_type
#
# DEPENDENCIES:
#   bedtools, awk, sort (GNU coreutils)
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail

DIFF_THRESH="0"

FEATURE_BED="allFeatures.feature.promoters_u400-TSS-d200.bed"
STRICT_BED="cDMRs_diff0.bed"

OUT_DIR="../plots_${DIFF_THRESH}"
mkdir -p "${OUT_DIR}"
OUT="${OUT_DIR}/cDMRs_diff${DIFF_THRESH}.closest_mRNA.tsv"


tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

A_BED="${tmpdir}/genes.mrna.bed"
B_BED="${tmpdir}/strict_dmrs.bed"
RAW="${tmpdir}/gene_to_dmr.closest.raw.tsv"
DMR_CANDIDATES="${tmpdir}/dmr_to_gene.candidates.tsv"

# ------------------------------------------------------------
# 0) Basic checks
# ------------------------------------------------------------
[[ -f "${FEATURE_BED}" ]] || { echo "ERROR: Missing feature BED: ${FEATURE_BED}"; exit 1; }
[[ -f "${STRICT_BED}"  ]] || { echo "ERROR: Missing strict BED: ${STRICT_BED}"; exit 1; }

# ------------------------------------------------------------
# 1) Build mRNA BED from genome feature file
#    A = genes/mRNA intervals (reference frame for signed distance)
#    Columns emitted:
#      chr start end gene 0 strand
# ------------------------------------------------------------
awk -F $'\t' -v OFS=$'\t' '
  $5 == "mRNA" {
    print $1, $2, $3, $4, 0, $6
  }
' "${FEATURE_BED}" | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "${A_BED}"

# ------------------------------------------------------------
# 2) Build a clean strict scDMR BED
#    B = strict scDMRs
#    Columns emitted:
#      chr start end cons_id dmr_score cluster
# ------------------------------------------------------------
awk -F $'\t' -v OFS=$'\t' '
  {
    split($4, a, "|")
    cons_id = a[1]
    cluster = a[2]
    print $1, $2, $3, cons_id, $5, cluster
  }
' "${STRICT_BED}" | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "${B_BED}"

# ------------------------------------------------------------
# 3) Closest strict scDMR to each gene/mRNA interval
#    -D a  => signed distance relative to A (the gene), which is what we want
#    -t first => choose one hit if tied
#
#    Interpretation of dist_signed:
#      negative = scDMR upstream of gene (relative to gene strand)
#      0        = scDMR overlaps gene body
#      positive = scDMR downstream of gene (relative to gene strand)
# ------------------------------------------------------------
bedtools closest -a "${A_BED}" -b "${B_BED}" -D a -t first > "${RAW}"

# ------------------------------------------------------------
# 4) Convert gene->DMR matches into a DMR-centered candidate table
#    Keep:
#      dmr fields
#      gene fields
#      dist_signed
#      abs_dist_bp
# ------------------------------------------------------------
awk -F $'\t' -v OFS=$'\t' '
  {
    # A (gene/mRNA): 1-6
    # B (DMR):       7-12
    # dist:          13

    if ($7=="." || $8==-1 || $9==-1) next

    dist_signed = $13
    abs_dist = (dist_signed < 0 ? -dist_signed : dist_signed)

    print $7, $8, $9, $10, $12, $11, \
          $1, $2, $3, $4, $6, \
          dist_signed, abs_dist
  }
' "${RAW}" > "${DMR_CANDIDATES}"

# ------------------------------------------------------------
# 5) For each DMR, keep the closest gene
#    Tie-breakers:
#      1) smallest absolute distance
#      2) smallest signed distance
#      3) lexicographic gene name
# ------------------------------------------------------------
{
  echo -e "dmr_chr\tdmr_start\tdmr_end\tcons_id\tcluster\tdmr_score\tgene_chr\tgene_start\tgene_end\tgene\tgene_strand\tdist_signed\tabs_dist_bp\tchromosome_type"
  
  LC_ALL=C sort -t $'\t' \
  -k1,1 -k2,2n -k3,3n -k4,4 \
  -k13,13n -k12,12n -k10,10 \
  "${DMR_CANDIDATES}" \
  | awk -F $'\t' -v OFS=$'\t' '
      {
        key = $1 OFS $2 OFS $3 OFS $4
        if (!(key in seen)) {
          seen[key] = 1
          chrom_type = ($1 == "Parae_12" ? "Sex_Ch" : "Autosome")
          print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,chrom_type
        }
      }
    '
} > "${OUT}"

echo "Wrote: ${OUT}"