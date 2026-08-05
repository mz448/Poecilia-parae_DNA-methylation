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
# VERSION: 04
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# V04: Only Nuclear chromosomes
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# DESCRIPTION:
#   Calculates CpG sequencing depth from DNMTools-style *.meth files.
#   - Do this in Symetric files! So the real CpG coverage is represented
#   - Zero-coverage CpG sites are INCLUDED in all calculations.
#   - Also reports Per sample SD and SEM
#
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

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "sample" \
    "reported_CpG_positions" \
    "mean_CpG_depth" \
    "SD_CpG_depth" \
    "SEM_CpG_depth" \
    "Q1_CpG_depth" \
    "Q2_CpG_depth" \
    "Q3_CpG_depth" \
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

    # Remove common file suffix
    SAMPLE=${SAMPLE%.meth}

    echo "Processing: $SAMPLE" >&2

    awk -v sample="$SAMPLE" '
        BEGIN {
            OFS = "\t"
        }
        
        $1 ~ /^Parae_(0[1-9]|1[0-9]|2[0-3])$/ && NF >= 6 {
            depth = $6
        
            n++
            total_depth += depth
            sumsq_depth += depth * depth
        
            depth_count[depth]++
        
            if (depth > max_depth) {
                max_depth = depth
            }
        
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
                printf "%s\t0\tNA\tNA\tNA\tNA\tNA\tNA\t0\tNA\t0\tNA\t0\tNA\n", sample
                exit
            }
        
            mean = total_depth / n
        
            if (n > 1) {
                variance = (sumsq_depth - (total_depth * total_depth / n)) / (n - 1)
        
                if (variance < 0) {
                    variance = 0
                }
        
                sd = sqrt(variance)
                sem = sd / sqrt(n)
            }
            else {
                sd = 0
                sem = 0
            }
        
            # Calculate Q1, Q2 (median), and Q3
            q1_target = n * 0.25
            q2_target = n * 0.50
            q3_target = n * 0.75
        
            cumulative = 0
            q1_found = 0
            q2_found = 0
            q3_found = 0
        
            for (d = 0; d <= max_depth; d++) {
        
                cumulative += depth_count[d]
        
                if (!q1_found && cumulative >= q1_target) {
                    q1 = d
                    q1_found = 1
                }
        
                if (!q2_found && cumulative >= q2_target) {
                    q2 = d
                    q2_found = 1
                }
        
                if (!q3_found && cumulative >= q3_target) {
                    q3 = d
                    q3_found = 1
                    break
                }
            }
        
            printf "%s\t%d\t%.3f\t%.3f\t%.6f\t%d\t%d\t%d\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n", \
                sample, \
                n, \
                mean, \
                sd, \
                sem, \
                q1, \
                q2, \
                q3, \
                cov5, 100 * cov5 / n, \
                cov10, 100 * cov10 / n, \
                cov20, 100 * cov20 / n
        }
        ' "$METH" >> "$OUT"
done

echo "Coverage summary written to: $OUT" >&2