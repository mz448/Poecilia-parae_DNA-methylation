#!/usr/bin/env bash

#SBATCH --job-name=cov_sym
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=12
#SBATCH --mem=36G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err
# NOTE: %x_%j is going to be the name of the job

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# AUTHOR: MZF
# DATE: 20260718
# SCRIPT: 03_Calculate_CpG_Coverage_meth_V01.sh
# VERSION: 01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# DESCRIPTION:
#   Calculate CpG sequencing depth from DNMTools-style *.meth files.
#   Do this in Symetric files! So the real CpG coverage is represented
#
#   Input columns:
#     1 = chromosome
#     2 = position
#     3 = strand
#     4 = sequence context
#     5 = methylation level
#     6 = number of reads overlapping the site
#
#   CpG depth = column 6
#
# USAGE:
#   ./Calculate_CpG_Coverage_meth.sh output.tsv sample1.meth [sample2.meth ...]
#
# EXAMPLE:
#   ./Calculate_CpG_Coverage_meth.sh CpG_coverage_summary.tsv \
#       Sym.forSym.pparymem009.CpG_report.txt.meth
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail

# ------------------------------------------------------------------------------
# Check arguments
# ------------------------------------------------------------------------------

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 output.tsv sample1.meth [sample2.meth ...]" >&2
    exit 1
fi

OUT="$1"
shift

# ------------------------------------------------------------------------------
# Write output header
# ------------------------------------------------------------------------------

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "sample" \
    "reported_CpG_positions" \
    "mean_CpG_depth" \
    "CpGs_ge_5x" \
    "percent_CpGs_ge_5x" \
    "CpGs_ge_10x" \
    "percent_CpGs_ge_10x" \
    "CpGs_ge_20x" \
    "percent_CpGs_ge_20x" \
    > "$OUT"

# ------------------------------------------------------------------------------
# Process each meth file
# ------------------------------------------------------------------------------

for METH in "$@"; do

    if [[ ! -f "$METH" ]]; then
        echo "ERROR: File not found: $METH" >&2
        exit 1
    fi

    SAMPLE=$(basename "$METH")

    # Remove common file suffixes
    SAMPLE=${SAMPLE%.meth}

    echo "Processing: $SAMPLE" >&2

    awk -v sample="$SAMPLE" '
    BEGIN {
        OFS = "\t"
    }

    NF >= 6 {
        depth = $6

        n++
        total_depth += depth

        if (depth >= 5) {
            cov5++
        }

        if (depth >= 10) {
            cov10++
        }

        if (depth >= 20) {
            cov20++
        }
    }

    END {
        if (n == 0) {
            printf "%s\t0\tNA\t0\tNA\t0\tNA\t0\tNA\n", sample
            exit
        }

        printf "%s\t%d\t%.3f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n", \
            sample, \
            n, \
            total_depth / n, \
            cov5, 100 * cov5 / n, \
            cov10, 100 * cov10 / n, \
            cov20, 100 * cov20 / n
    }
    ' "$METH" >> "$OUT"

done

echo "Coverage summary written to: $OUT" >&2