#!/bin/bash
#SBATCH --job-name=01_differentialMethylation_methylasso.sh
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=2
#SBATCH --mem=500G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL: To perform differential methylation detection 
#       • Using emseq samples that have been processed to become symetrical at 
#         CpG sites (from BISMARk alignment and then DNMTools SYM func.)
#       • The detection is done in contrasts keeping the correct genotype, as well
#         as in contrasts with shuffled samples for sex or for morph 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

# Set variables
SIF="/programs/methylasso-1/methylasso.sif"
SCRIPT="/MethyLasso.R"
OUTBASE="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/02_DMR/data/01_methylasso/methylasso_DMRs"

# FUNCTION: run_methylasso ()
# → Runs differential methylation analysis for two groups and also outputs individual profiles
# → Takes:  group_name 1, comma-separated sample list 1, 
#           group_name 2, comma-separated sample list 2, output subfolder
run_methylasso () {
  group_name1="$1"
  samples1="$2"
  group_name2="$3"
  samples2="$4"
  outdir="$5"
  mkdir -p "$outdir"
  /usr/bin/time -v singularity exec --bind $PWD --pwd $PWD "$SIF" \
    Rscript "$SCRIPT" \
    --n1 "$group_name1" --c1 "$samples1" \
    --n2 "$group_name2" --c2 "$samples2"  \
    --cov 5 \
    --meth 4 \
    -c 5 \
    -d 0.1  \
    -p 0.0001  \
    --max_distance 500  \
    --min_width 100 \
    -o "$outdir" \
    2>> "$outdir/memory_usage.log" &
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ real_sex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso f "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
               p "methylasso.pmem004.txt,methylasso.pmem005.txt,methylasso.pmem006.txt" \
               $OUTBASE/real_sex/f_vs_p
wait
run_methylasso f "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
               y "methylasso.ymem007.txt,methylasso.ymem008.txt,methylasso.ymem009.txt" \
               $OUTBASE/real_sex/f_vs_y
wait
run_methylasso f "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
               i "methylasso.imem010.txt,methylasso.imem011.txt,methylasso.imem012.txt" \
               $OUTBASE/real_sex/f_vs_i
wait
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ real_morph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso p "methylasso.pmem004.txt,methylasso.pmem005.txt,methylasso.pmem006.txt" \
               i_01 "methylasso.imem010.txt,methylasso.imem011.txt,methylasso.imem012.txt" \
               $OUTBASE/real_morph/p_vs_i_01
wait
run_methylasso y "methylasso.ymem007.txt,methylasso.ymem008.txt,methylasso.ymem009.txt" \
               p_01 "methylasso.pmem004.txt,methylasso.pmem005.txt,methylasso.pmem006.txt" \
               $OUTBASE/real_morph/y_vs_p_01
wait
run_methylasso i "methylasso.imem010.txt,methylasso.imem011.txt,methylasso.imem012.txt" \
               y_01 "methylasso.ymem007.txt,methylasso.ymem008.txt,methylasso.ymem009.txt" \
               $OUTBASE/real_morph/i_vs_y_01
wait

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ shuffled_sex ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso rf "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
               sM_01 "methylasso.imem010.txt,methylasso.pmem005.txt,methylasso.ymem009.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_01
wait
run_methylasso rf "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
               sM_02 "methylasso.pmem004.txt,methylasso.ymem008.txt,methylasso.imem012.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_02
wait
run_methylasso rf "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
               sM_03 "methylasso.ymem007.txt,methylasso.imem011.txt,methylasso.pmem006.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_03
wait
run_methylasso rf "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
               sM_04 "methylasso.imem010.txt,methylasso.ymem007.txt,methylasso.pmem006.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_04
wait
run_methylasso rf "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
               sM_05 "methylasso.pmem004.txt,methylasso.imem011.txt,methylasso.ymem009.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_05
wait
run_methylasso rf "methylasso.fmem001.txt,methylasso.fmem002.txt,methylasso.fmem003.txt" \
               sM_06 "methylasso.ymem007.txt,methylasso.pmem005.txt,methylasso.imem012.txt" \
               $OUTBASE/shuffled_sex/f_vs_sM_06
wait

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ shuffled_morph ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso sM "methylasso.imem010.txt,methylasso.pmem005.txt,methylasso.ymem009.txt" \
               sM_01 "methylasso.ymem007.txt,methylasso.imem011.txt,methylasso.pmem006.txt" \
               $OUTBASE/shuffled_morph/sM_vs_sM_01
wait
run_methylasso sM "methylasso.pmem004.txt,methylasso.ymem007.txt,methylasso.imem012.txt" \
               sM_02 "methylasso.imem010.txt,methylasso.pmem005.txt,methylasso.ymem009.txt" \
               $OUTBASE/shuffled_morph/sM_vs_sM_02
wait
run_methylasso sM "methylasso.ymem007.txt,methylasso.imem011.txt,methylasso.pmem006.txt" \
               sM_03 "methylasso.pmem004.txt,methylasso.ymem008.txt,methylasso.imem012.txt" \
               $OUTBASE/shuffled_morph/sM_vs_sM_03
wait

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ shuffled_All ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
run_methylasso sA "methylasso.fmem001.txt,methylasso.imem011.txt,methylasso.ymem009.txt" \
               sA_01 "methylasso.pmem004.txt,methylasso.fmem002.txt,methylasso.imem012.txt" \
               $OUTBASE/shuffled_all/sA_vs_sA_01
wait
run_methylasso sA "methylasso.imem010.txt,methylasso.fmem002.txt,methylasso.pmem006.txt" \
               sA_02 "methylasso.ymem007.txt,methylasso.pmem005.txt,methylasso.fmem003.txt" \
               $OUTBASE/shuffled_all/sA_vs_sA_02
wait
run_methylasso sA "methylasso.pmem004.txt,methylasso.imem011.txt,methylasso.fmem003.txt" \
               sA_03 "methylasso.fmem001.txt,methylasso.ymem007.txt,methylasso.pmem006.txt" \
               $OUTBASE/shuffled_all/sA_vs_sA_03
wait

