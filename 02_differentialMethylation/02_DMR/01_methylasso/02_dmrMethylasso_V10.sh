#!/bin/bash
#SBATCH --job-name=02_dmrMethylasso_V10
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=4
#SBATCH --mem=300G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err


# DETAIL: Run DMR analysis using symmetric files, q=0.05, cov >= 3
# This script specifies source file format to correctly identify the columns 
OUTPREFIX="q0.05_cov3_sym_"

# Set variables
SIF="/programs/methylasso-1/methylasso.sif"
SCRIPT="/MethyLasso.R"
BASE_DIR="/local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/02_differential_methylation/03_methylasso_DMRs/data"
OUTBASE="${BASE_DIR}/dmr_output/${OUTPREFIX}run"

mkdir -p "$OUTBASE"

# Helper function for readability
run_methylasso () {
  name1="$1"
  group1="$2"
  name2="$3"
  group2="$4"
  logdir="$5"
  mkdir -p "$logdir"
  /usr/bin/time -v singularity exec --bind $PWD --pwd $PWD "$SIF" \
    Rscript "$SCRIPT" \
      --n1 "$name1" --c1 "$group1" \
      --n2 "$name2" --c2 "$group2" \
      --cov 5 --meth 4 \
      -q 0.05 -c 3 -t 2 \
      -o "$logdir" \
    2>> "$logdir/memory_usage.log" &
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ real_sex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso f "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               p "sym.formatted.pparpmem004.CpG_report.txt,sym.formatted.pparpmem005.CpG_report.txt,sym.formatted.pparpmem006.CpG_report.txt" \
               $OUTBASE/real_sex/f_vs_p
wait
run_methylasso f "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               y "sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparymem008.CpG_report.txt,sym.formatted.pparymem009.CpG_report.txt" \
               $OUTBASE/real_sex/f_vs_y
wait
run_methylasso f "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               i "sym.formatted.pparimem010.CpG_report.txt,sym.formatted.pparimem011.CpG_report.txt,sym.formatted.pparimem012.CpG_report.txt" \
               $OUTBASE/real_sex/f_vs_i
wai
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ real_morph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso p "sym.formatted.pparpmem004.CpG_report.txt,sym.formatted.pparpmem005.CpG_report.txt,sym.formatted.pparpmem006.CpG_report.txt" \
               i_01 "sym.formatted.pparimem010.CpG_report.txt,sym.formatted.pparimem011.CpG_report.txt,sym.formatted.pparimem012.CpG_report.txt" \
               $OUTBASE/real_morph/p_vs_i_01
wait
run_methylasso y "sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparymem008.CpG_report.txt,sym.formatted.pparymem009.CpG_report.txt" \
               p_01 "sym.formatted.pparpmem004.CpG_report.txt,sym.formatted.pparpmem005.CpG_report.txt,sym.formatted.pparpmem006.CpG_report.txt" \
               $OUTBASE/real_morph/y_vs_p_01
wait
run_methylasso i "sym.formatted.pparimem010.CpG_report.txt,sym.formatted.pparimem011.CpG_report.txt,sym.formatted.pparimem012.CpG_report.txt" \
               y_01 "sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparymem008.CpG_report.txt,sym.formatted.pparymem009.CpG_report.txt" \
               $OUTBASE/real_morph/i_vs_y_01
wait

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ shuffled_sex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso rf "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               sM_01 "sym.formatted.pparimem010.CpG_report.txt,sym.formatted.pparpmem005.CpG_report.txt,sym.formatted.pparymem009.CpG_report.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_01
wait
run_methylasso rf "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               sM_02 "sym.formatted.pparpmem004.CpG_report.txt,sym.formatted.pparymem008.CpG_report.txt,sym.formatted.pparimem012.CpG_report.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_02
wait
run_methylasso rf "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               sM_03 "sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparimem011.CpG_report.txt,sym.formatted.pparpmem006.CpG_report.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_03
wait
run_methylasso rf "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               sM_04 "sym.formatted.pparimem010.CpG_report.txt,sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparpmem006.CpG_report.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_04
wait
run_methylasso rf "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               sM_05 "sym.formatted.pparpmem004.CpG_report.txt,sym.formatted.pparimem011.CpG_report.txt,sym.formatted.pparymem009.CpG_report.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_05
wait
run_methylasso f "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               sM_06 "sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparpmem005.CpG_report.txt,sym.formatted.pparimem012.CpG_report.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_06
wait

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ shuffled_morph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso sM "sym.formatted.pparimem010.CpG_report.txt,sym.formatted.pparpmem005.CpG_report.txt,sym.formatted.pparymem009.CpG_report.txt" \
               sM_01 "sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparimem011.CpG_report.txt,sym.formatted.pparpmem006.CpG_report.txt" \
               $OUTBASE/shuffled_morph/sM_vs_sM_01
wait
run_methylasso sM "sym.formatted.pparpmem004.CpG_report.txt,sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparimem012.CpG_report.txt" \
               sM_02 "sym.formatted.pparimem010.CpG_report.txt,sym.formatted.pparpmem005.CpG_report.txt,sym.formatted.pparymem009.CpG_report.txt" \
               $OUTBASE/shuffled_morph/sM_vs_sM_02
wait
run_methylasso sM "sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparimem011.CpG_report.txt,sym.formatted.pparpmem006.CpG_report.txt" \
               sM_03 "sym.formatted.pparpmem004.CpG_report.txt,sym.formatted.pparymem008.CpG_report.txt,sym.formatted.pparimem012.CpG_report.txt" \
               $OUTBASE/shuffled_morph/sM_vs_sM_03
wait

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ shuffled_All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso sA "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparimem011.CpG_report.txt,sym.formatted.pparymem009.CpG_report.txt" \
               sA_01 "sym.formatted.pparpmem004.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparimem012.CpG_report.txt" \
               $OUTBASE/shuffled_all/sA_vs_sA_01
wait
run_methylasso sA "sym.formatted.pparimem010.CpG_report.txt,sym.formatted.pparfmem002.CpG_report.txt,sym.formatted.pparpmem006.CpG_report.txt" \
               sA_02 "sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparpmem005.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               $OUTBASE/shuffled_all/sA_vs_sA_02
wait
run_methylasso sA "sym.formatted.pparpmem004.CpG_report.txt,sym.formatted.pparimem011.CpG_report.txt,sym.formatted.pparfmem003.CpG_report.txt" \
               sA_03 "sym.formatted.pparfmem001.CpG_report.txt,sym.formatted.pparymem007.CpG_report.txt,sym.formatted.pparpmem006.CpG_report.txt" \
               $OUTBASE/shuffled_all/sA_vs_sA_03
wait

