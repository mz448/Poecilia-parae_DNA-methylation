#!/usr/bin/env bash
# DATE: 20260315
# AUTHOR: MZF
# UPDATED: 20260609
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Annotate each VMR with its closest mRNA feature.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# BIOLOGICAL CONTEXT:
#   This script assigns one closest gene model to every input VMR.
#   Distances are calculated relative to the mRNA interval and are
#   strand-aware:
#
#     distance = 0  : VMR overlaps the mRNA body
#     distance < 0  : VMR is upstream of the gene
#     distance > 0  : VMR is downstream of the gene
#
#   For genes on the negative strand, upstream and downstream are interpreted
#   relative to transcriptional orientation.
#
# INPUTS:
#   1) allFeatures.feature.promoters_u400-TSS-d200.bed
#      Genome feature BED. Only rows with feature == "mRNA" are used.
#      Expected columns:
#        chr, start, end, gene, feature, strand, length
#
#   2) VMRs_diff0.bed
#      BED6 file of VMRs/VMRs.
#      Expected columns:
#        chr, start, end, CONSxxxxxxx|Cluster_X, score, strand
#
# OUTPUT:
#   VMRs_diff0.closest_mRNA.tsv
#      One row per input VMR/VMR, plus header.
#      Columns:
#        dmr_chr, dmr_start, dmr_end, cons_id, cluster, dmr_score,
#        gene_chr, gene_start, gene_end, gene, gene_strand,
#        dist_signed, abs_dist_bp, chromosome_type
#
# METHOD:
#   1) Extract mRNA intervals from the genome feature BED.
#   2) Clean VMR IDs by separating consensus ID and cluster.
#   3) Use bedtools closest with VMRs as A and mRNAs as B.
#   4) Use -D b so signed distances are relative to gene strand.
#   5) Retain one closest mRNA per VMR using explicit tie-breaking.
#   6) Add chromosome type: Parae_12 = Sex_Ch; all others = Autosome.
#
# IMPORTANT:
#   VMRs/VMRs must be used as the -a file in bedtools closest.
#   Using genes as -a changes the question to "closest VMR per gene" and
#   causes VMRs that are not closest to any gene to be lost.
#
# DEPENDENCIES:
#   bedtools, awk, sort
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail

DIFF_THRESH="0"

FEATURE_BED="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/data/features/allFeatures.feature.bed"
VMR_BED="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/03_VMR/plots/0/clusters_WG.sampleheatmap.diffGE_0.bed"

OUT_DIR="../plots/${DIFF_THRESH}"
mkdir -p "${OUT_DIR}"
OUT="${OUT_DIR}/VMRs_diff${DIFF_THRESH}.closest_mRNA.tsv"


tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

A_BED="${tmpdir}/genes.mrna.bed"
B_BED="${tmpdir}/VMRs.bed"
RAW="${tmpdir}/gene_to_VMR.closest.raw.tsv"
DMR_CANDIDATES="${tmpdir}/VMR_to_gene.candidates.tsv"

# ------------------------------------------------------------
# 0) Basic checks
# ------------------------------------------------------------
[[ -f "${FEATURE_BED}" ]] || { echo "ERROR: Missing feature BED: ${FEATURE_BED}"; exit 1; }
[[ -f "${VMR_BED}"  ]] || { echo "ERROR: Missing strict BED: ${VMR_BED}"; exit 1; }

# ------------------------------------------------------------
# 1) Build mRNA BED from genome feature file
#    A = genes/mRNA intervals (reference frame for signed distance)
#    Columns emitted:
#      chr start end gene 0 strand
# ------------------------------------------------------------
awk -F $'\t' -v OFS=$'\t' '
  $5 == "gene" {
    print $1, $2, $3, $4, 0, $6
  }
' "${FEATURE_BED}" | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "${A_BED}"

# ------------------------------------------------------------
# 2) Build a clean VMR BED
#    B = VMRs
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
' "${VMR_BED}" | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "${B_BED}"

# ------------------------------------------------------------
# 3) Closest mRNA interval to each VMR
#    A = VMRs
#    B = genes/mRNA intervals
#
#    -D b  => signed distance relative to the gene/mRNA strand
#    -t all => keep all ties so our downstream tie-breaker can choose
# ------------------------------------------------------------
bedtools closest -a "${B_BED}" -b "${A_BED}" -D b -t all > "${RAW}"

# ------------------------------------------------------------
# 4) Convert DMR->gene matches into a DMR-centered candidate table
# ------------------------------------------------------------
awk -F $'\t' -v OFS=$'\t' '
  {
    # A (DMR):       1-6
    # B (gene/mRNA): 7-12
    # dist:          13

    if ($7=="." || $8==-1 || $9==-1) {
      print $1, $2, $3, $4, $6, $5, \
            "NA", "NA", "NA", "NA", "NA", \
            "NA", "NA"
      next
    }

    dist_signed = $13
    abs_dist = (dist_signed < 0 ? -dist_signed : dist_signed)

    print $1, $2, $3, $4, $6, $5, \
          $7, $8, $9, $10, $12, \
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