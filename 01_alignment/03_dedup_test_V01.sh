#!/usr/bin/env bash


# DATE: 20260717
# AUTHOR: MZF
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# SCRIPT:   03_dedup_test_V01.sh
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# DESCRIPTION:
#   This script calculates depth from the 3 independent alingments of the same 
#   biological sample.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Conclussion:  I was manualy mismatiching the samples by name. The 2 runs are 
#               additive when aligned together using Bismark!
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# inside /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bismarkAlignment/bismark/bismark/deduplicated
# Run:  bash ../../../../scripts/03_dedup_test_V01.sh > pparimem012_allLibraries_dedupDepth.tsv
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

for BAM in \
    /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/01_nfcore_methylseq/bismark03/bismark/deduplicated/pparimem012.deduplicated.sorted.bam \
    /local/storage/Projects/ppar_emseq/analysis/010_emseq_pparae_muscle/data/01_alignment/bismark/bismark/deduplicated/pparimem012.deduplicated.sorted.bam \
    /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bismarkAlignment/bismark/bismark/deduplicated/pparfmem012.deduplicated.sorted.bam
do

    echo "$BAM"

    samtools depth -aa -s "$BAM" |
    awk '
    $1 ~ /^Parae_(0[1-9]|1[0-9]|2[0-3])$/ {
        n++
        sum += $3
    }
    END {
        print "nuclear_bases =", n
        print "total_depth   =", sum
        print "mean_depth    =", sum/n
    }'

done