#!/usr/bin/env bash
#SBATCH --job-name=01_differentialMethylation_methylasso
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=2
#SBATCH --mem=500G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:
#   Perform differential methylation detection using EM-seq samples processed
#   to be symmetric at CpG sites.
#
#   DMRs are calculated using:
#     1) Real sex contrasts
#     2) Real male-morph contrasts
#     3) Sex-preserved contrasts with shuffled male-morph composition
#     4) Male-morph-shuffled contrasts
#     5) Contrasts with sex and morph identities shuffled together
#
# Input:
#   MethyLasso-formatted methylation files:
#     methylasso.<sample>.txt
#
# Output:
#   One output directory per contrast under OUTBASE.
#
# Corrections:
#   - shuffled_sex sM_06 now uses ymem008 instead of duplicated ymem007.
#   - shuffled_morph sM_vs_sM_02 now uses ymem008 instead of duplicated ymem007.
#   - The corrected assignments ensure balanced use of all nine male samples
#     in the intended mixed-morph partitions.
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail


# =============================================================================
# Global parameters
# =============================================================================

SIF="/programs/methylasso-1/methylasso.sif"
SCRIPT="/MethyLasso.R"

OUTBASE="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/02_DMR/data/01_methylasso/methylasso_DMRs_V02"


# =============================================================================
# FUNCTION: run_methylasso
#
# Arguments:
#   1. Group 1 name
#   2. Group 1 comma-separated sample list
#   3. Group 2 name
#   4. Group 2 comma-separated sample list
#   5. Output directory
# =============================================================================

run_methylasso() {

  local group_name1="$1"
  local samples1="$2"
  local group_name2="$3"
  local samples2="$4"
  local outdir="$5"

  mkdir -p "${outdir}"

  /usr/bin/time -v \
    singularity exec \
      --bind "${PWD}" \
      --pwd "${PWD}" \
      "${SIF}" \
      Rscript "${SCRIPT}" \
        --n1 "${group_name1}" \
        --c1 "${samples1}" \
        --n2 "${group_name2}" \
        --c2 "${samples2}" \
        --cov 5 \
        --meth 4 \
        -c 5 \
        -d 0.1 \
        -p 0.0001 \
        --max_distance 500 \
        --min_width 100 \
        -o "${outdir}" \
        2>> "${outdir}/memory_usage.log" &
}


# =============================================================================
# real_sex
#
# Biological females are compared separately against each male morph.
# =============================================================================

run_methylasso \
  f \
  "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
  p \
  "methylasso.pmem004.txt,methylasso.pmem005.txt,methylasso.pmem006.txt" \
  "${OUTBASE}/real_sex/f_vs_p"

wait


run_methylasso \
  f \
  "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
  y \
  "methylasso.ymem007.txt,methylasso.ymem008.txt,methylasso.ymem009.txt" \
  "${OUTBASE}/real_sex/f_vs_y"

wait


run_methylasso \
  f \
  "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
  i \
  "methylasso.imem010.txt,methylasso.imem011.txt,methylasso.imem012.txt" \
  "${OUTBASE}/real_sex/f_vs_i"

wait


# =============================================================================
# real_morph
#
# Pairwise comparisons among the three biological male morphs.
# =============================================================================

run_methylasso \
  p \
  "methylasso.pmem004.txt,methylasso.pmem005.txt,methylasso.pmem006.txt" \
  i_01 \
  "methylasso.imem010.txt,methylasso.imem011.txt,methylasso.imem012.txt" \
  "${OUTBASE}/real_morph/p_vs_i_01"

wait


run_methylasso \
  y \
  "methylasso.ymem007.txt,methylasso.ymem008.txt,methylasso.ymem009.txt" \
  p_01 \
  "methylasso.pmem004.txt,methylasso.pmem005.txt,methylasso.pmem006.txt" \
  "${OUTBASE}/real_morph/y_vs_p_01"

wait


run_methylasso \
  i \
  "methylasso.imem010.txt,methylasso.imem011.txt,methylasso.imem012.txt" \
  y_01 \
  "methylasso.ymem007.txt,methylasso.ymem008.txt,methylasso.ymem009.txt" \
  "${OUTBASE}/real_morph/i_vs_y_01"

wait


# =============================================================================
# shuffled_sex
#
# Sex is preserved:
#   - Group 1 always contains the three biological females.
#   - Group 2 always contains three biological males.
#
# Male-morph composition is shuffled:
#   - Each artificial male group contains one immaculata, one parae and
#     one yellow individual.
#
# Two complete balanced partitions are used.
#
# Partition 1:
#   sM_01 = i10 + p05 + y09
#   sM_02 = i12 + p04 + y08
#   sM_03 = i11 + p06 + y07
#
# Partition 2:
#   sM_04 = i10 + p06 + y07
#   sM_05 = i11 + p04 + y09
#   sM_06 = i12 + p05 + y08
# =============================================================================

run_methylasso \
  rf \
  "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
  sM_01 \
  "methylasso.imem010.txt,methylasso.pmem005.txt,methylasso.ymem009.txt" \
  "${OUTBASE}/shuffled_sex/f_vs_sM_01"

wait


run_methylasso \
  rf \
  "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
  sM_02 \
  "methylasso.pmem004.txt,methylasso.ymem008.txt,methylasso.imem012.txt" \
  "${OUTBASE}/shuffled_sex/f_vs_sM_02"

wait


run_methylasso \
  rf \
  "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
  sM_03 \
  "methylasso.ymem007.txt,methylasso.imem011.txt,methylasso.pmem006.txt" \
  "${OUTBASE}/shuffled_sex/f_vs_sM_03"

wait


run_methylasso \
  rf \
  "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
  sM_04 \
  "methylasso.imem010.txt,methylasso.ymem007.txt,methylasso.pmem006.txt" \
  "${OUTBASE}/shuffled_sex/f_vs_sM_04"

wait


run_methylasso \
  rf \
  "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
  sM_05 \
  "methylasso.pmem004.txt,methylasso.imem011.txt,methylasso.ymem009.txt" \
  "${OUTBASE}/shuffled_sex/f_vs_sM_05"

wait


run_methylasso \
  rf \
  "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
  sM_06 \
  "methylasso.ymem008.txt,methylasso.pmem005.txt,methylasso.imem012.txt" \
  "${OUTBASE}/shuffled_sex/f_vs_sM_06"

wait


# =============================================================================
# shuffled_morph
#
# All samples are biological males.
#
# Three artificial groups are created, each containing:
#   - one immaculata
#   - one parae
#   - one yellow
#
# Artificial groups:
#   A = i10 + p05 + y09
#   B = i11 + p06 + y07
#   C = i12 + p04 + y08
#
# The three contrasts form the complete pairwise comparison:
#   A versus B
#   C versus A
#   B versus C
# =============================================================================

# A versus B
run_methylasso \
  sM_A \
  "methylasso.imem010.txt,methylasso.pmem005.txt,methylasso.ymem009.txt" \
  sM_B \
  "methylasso.ymem007.txt,methylasso.imem011.txt,methylasso.pmem006.txt" \
  "${OUTBASE}/shuffled_morph/sM_vs_sM_01"

wait


# C versus A
run_methylasso \
  sM_C \
  "methylasso.pmem004.txt,methylasso.ymem008.txt,methylasso.imem012.txt" \
  sM_A \
  "methylasso.imem010.txt,methylasso.pmem005.txt,methylasso.ymem009.txt" \
  "${OUTBASE}/shuffled_morph/sM_vs_sM_02"

wait


# B versus C
run_methylasso \
  sM_B \
  "methylasso.ymem007.txt,methylasso.imem011.txt,methylasso.pmem006.txt" \
  sM_C \
  "methylasso.pmem004.txt,methylasso.ymem008.txt,methylasso.imem012.txt" \
  "${OUTBASE}/shuffled_morph/sM_vs_sM_03"

wait


# =============================================================================
# shuffled_all
#
# Both sex and male-morph identities are mixed across artificial groups.
# Both sides contain one female and two males!
# Female pairings form a complete three-contrast cycle:
#   f001 vs f002
#   f002 vs f003
#   f003 vs f001
#
# =============================================================================

run_methylasso \
  sA \
  "methylasso.fmem001.txt,methylasso.imem011.txt,methylasso.ymem009.txt" \
  sA_01 \
  "methylasso.pmem004.txt,methylasso.fmem002.txt,methylasso.imem012.txt" \
  "${OUTBASE}/shuffled_all/sA_vs_sA_01"

wait


run_methylasso \
  sA \
  "methylasso.imem010.txt,methylasso.fmem002.txt,methylasso.pmem006.txt" \
  sA_02 \
  "methylasso.ymem007.txt,methylasso.pmem005.txt,methylasso.fmem003.txt" \
  "${OUTBASE}/shuffled_all/sA_vs_sA_02"

wait


run_methylasso \
  sA \
  "methylasso.pmem004.txt,methylasso.imem011.txt,methylasso.fmem003.txt" \
  sA_03 \
  "methylasso.fmem001.txt,methylasso.ymem007.txt,methylasso.pmem006.txt" \
  "${OUTBASE}/shuffled_all/sA_vs_sA_03"

wait


# Alernative Shuffled All (NOT IMPLEMENTED in this analysis, bc is morph imbalanced)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ shuffled_all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Both sides contain one female and two males.
# Female pairings form a complete three-contrast cycle:
#   f001 vs f002
#   f002 vs f003
#   f003 vs f001
#
# Both sides of each contrast have identical phenotype composition.
# Every sample appears at least once.
# One male from each morph appears twice, once on each contrast side.
# ---------------------------------------------------------------------------
# 
# # Female + immaculata + parae
# run_methylasso \
#   sA_01A \
#   "methylasso.fmem001.txt,methylasso.imem010.txt,methylasso.pmem004.txt" \
#   sA_01B \
#   "methylasso.fmem002.txt,methylasso.imem011.txt,methylasso.pmem005.txt" \
#   "$OUTBASE/shuffled_all/sA_vs_sA_01"
# 
# wait
# 
# 
# # Female + immaculata + yellow
# run_methylasso \
#   sA_02A \
#   "methylasso.fmem002.txt,methylasso.imem012.txt,methylasso.ymem007.txt" \
#   sA_02B \
#   "methylasso.fmem003.txt,methylasso.imem010.txt,methylasso.ymem008.txt" \
#   "$OUTBASE/shuffled_all/sA_vs_sA_02"
# 
# wait
# 
# 
# # Female + parae + yellow
# run_methylasso \
#   sA_03A \
#   "methylasso.fmem003.txt,methylasso.pmem006.txt,methylasso.ymem009.txt" \
#   sA_03B \
#   "methylasso.fmem001.txt,methylasso.pmem004.txt,methylasso.ymem007.txt" \
#   "$OUTBASE/shuffled_all/sA_vs_sA_03"

wait


echo "=== MethyLasso differential methylation analyses completed ==="
