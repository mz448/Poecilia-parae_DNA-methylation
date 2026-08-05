#!/usr/bin/env bash
#SBATCH --job-name=nuclear_depth
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=98G
#SBATCH --output=nuclear_depth_%j.out
#SBATCH --error=nuclear_depth_%j.err

# DATE: 20260717
# AUTHOR: MZF
# ==============================================================================
# SCRIPT: Calculate_Nuclear_Depth_parallel.sh
#
# DESCRIPTION:
#   Calculate nuclear genome depth and breadth of coverage from multiple
#   deduplicated Bismark BAM files in parallel.
#
# USAGE:
#   ./Calculate_Nuclear_Depth_parallel.sh output.tsv sample1.bam [sample2.bam ...]
#
# NUMBER OF PARALLEL JOBS:
#   1. Uses $JOBS when defined.
#   2. Otherwise uses $SLURM_CPUS_PER_TASK inside a SLURM job.
#   3. Otherwise defaults to 4.
#
# EXAMPLE:
#   JOBS=8 ./Calculate_Nuclear_Depth_parallel.sh \
#       Nuclear_depth_summary.tsv \
#       *.deduplicated.bam
# ==============================================================================

set -euo pipefail

# ------------------------------------------------------------------------------
# Check arguments and required programs
# ------------------------------------------------------------------------------

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 output.tsv sample1.bam [sample2.bam ...]" >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools was not found in PATH." >&2
    exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
    echo "ERROR: GNU parallel was not found in PATH." >&2
    exit 1
fi

OUT="$1"
shift

# for debugging
echo "Number of BAM files: $#" >&2
printf 'BAM: %s\n' "$@" >&2

# Use JOBS, SLURM CPUs, or 4 jobs by default
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-4}}"

if ! [[ "$JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: JOBS must be a positive integer." >&2
    exit 1
fi

# ------------------------------------------------------------------------------
# Write output header
# ------------------------------------------------------------------------------

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "sample" \
    "nuclear_reference_bases" \
    "mean_nuclear_depth" \
    "bases_ge_1x" \
    "percent_bases_ge_1x" \
    "bases_ge_5x" \
    "percent_bases_ge_5x" \
    "bases_ge_10x" \
    "percent_bases_ge_10x" \
    "bases_ge_20x" \
    "percent_bases_ge_20x" \
    > "$OUT"

# ------------------------------------------------------------------------------
# Function for processing one BAM
# ------------------------------------------------------------------------------

process_bam() {

    set -o pipefail

    local BAM="$1"
    local SAMPLE

    if [[ ! -f "$BAM" ]]; then
        echo "ERROR: File not found: $BAM" >&2
        return 1
    fi

    SAMPLE=$(basename "$BAM")

    # Remove common BAM suffixes
    SAMPLE=${SAMPLE%.bam}
    SAMPLE=${SAMPLE%.deduplicated}

    echo "Processing: $SAMPLE" >&2

    samtools depth -aa -s "$BAM" |
    awk -v sample="$SAMPLE" '
    BEGIN {
        OFS = "\t"
    }

    $1 ~ /^Parae_(0[1-9]|1[0-9]|2[0-3])$/ {
        depth = $3

        bases++
        total_depth += depth

        if (depth >= 1) {
            cov1++
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
        if (bases == 0) {
            printf "%s\t0\tNA\t0\tNA\t0\tNA\t0\tNA\t0\tNA\n", sample
            exit
        }

        printf "%s\t%d\t%.3f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\t%d\t%.2f\n", \
            sample, \
            bases, \
            total_depth / bases, \
            cov1, 100 * cov1 / bases, \
            cov5, 100 * cov5 / bases, \
            cov10, 100 * cov10 / bases, \
            cov20, 100 * cov20 / bases
    }'
}

export -f process_bam

# ------------------------------------------------------------------------------
# Process BAM files in parallel
# ------------------------------------------------------------------------------

echo "Running up to $JOBS BAM files simultaneously." >&2

# parallel \
#     --jobs "$JOBS" \
#     --keep-order \
#     --halt soon,fail=1 \
#     process_bam ::: "$@" \
#     >> "$OUT"

parallel \
    --jobs "$JOBS" \
    --keep-order \
    --tag \
    process_bam ::: "$@" \
    >> "$OUT"

echo "Nuclear depth summary written to: $OUT" >&2