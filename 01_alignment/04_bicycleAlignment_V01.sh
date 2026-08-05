#!/bin/bash
#SBATCH --job-name=bicycle
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=80
#SBATCH --mem=900G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err


# DATE: 20260804
# AUTHOR: MZF
# Script: bicycleAlignment.sh
# Version 01
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# GOAL:   Run the whole Bicycle pipeline 
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# PRODECURE: 
# 1. Prepare Files
# 2. Create project
# 3. Create the Watson and Crick in-silico bisulfited reference genomes
# 4. Create the reference index
# 5. Align to the reference genome
# 6. Perform methylation analysis and methylcytosine calling
# 7. Analyze differential methylation
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail


BICYCLEPATH="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment"

# 1) set-up folders and files ~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

mkdir ${BICYCLEPATH}
cd ${BICYCLEPATH}
# Create necessary folder structure and Populate with samples:
# Make the bicycle directories
mkdir ${BICYCLEPATH}/ref_genomes 
mkdir ${BICYCLEPATH}/reads 
# Copy the Genome 
cp /local/storage/Projects/ppar_emseq/data/011_emseq_pparae_muscle/genome/PparFemMitoLambdaPuc19Ver2024.fasta ./ref_genomes/PparFemMitoLambdaPuc19Ver2024.fa

# ------------------------------------------------------------------------------
# Reads source folder
READS_SOURCE="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bismarkAlignment/bismark/cat"
READS="${BICYCLEPATH}/reads"

# Make all the individual sample directories
mkdir ${READS}/pparfmem001 
mkdir ${READS}/pparfmem002 
mkdir ${READS}/pparfmem003 
mkdir ${READS}/pparpmem004
mkdir ${READS}/pparpmem005
mkdir ${READS}/pparpmem006
mkdir ${READS}/pparymem007
mkdir ${READS}/pparymem008
mkdir ${READS}/pparymem009
mkdir ${READS}/pparimem010
mkdir ${READS}/pparimem011
mkdir ${READS}/pparimem012

# copy samples into their own folder
echo "======= COPYING SAMPLES ======="
cp ${READS_SOURCE}/pparfmem001* ${READS}/pparfmem001/
cp ${READS_SOURCE}/pparfmem002* ${READS}/pparfmem002/
cp ${READS_SOURCE}/pparfmem003* ${READS}/pparfmem003/
cp ${READS_SOURCE}/pparfmem004* ${READS}/pparpmem004/
cp ${READS_SOURCE}/pparfmem005* ${READS}/pparpmem005/
cp ${READS_SOURCE}/pparfmem006* ${READS}/pparpmem006/
cp ${READS_SOURCE}/pparfmem007* ${READS}/pparymem007/
cp ${READS_SOURCE}/pparfmem008* ${READS}/pparymem008/
cp ${READS_SOURCE}/pparfmem009* ${READS}/pparymem009/
cp ${READS_SOURCE}/pparfmem010* ${READS}/pparimem010/
cp ${READS_SOURCE}/pparfmem011* ${READS}/pparimem011/
cp ${READS_SOURCE}/pparfmem012* ${READS}/pparimem012/
echo "======= [✓] SAMPLES COPIED  =======" &&

# Unzip samples
# echo "======= UNZIPING SAMPLES ======="
# date
# # gunzip ${READS}/*/pparfmem*.fastq.gz
# # gunzip ${READS}/*/pparpmem*.fastq.gz 
# # gunzip ${READS}/*/pparymem*.fastq.gz 
# # gunzip ${READS}/*/pparimem*.fastq.gz
# 
# 
# echo "======= SAMPLES READY ======="
# date

# Decompress FASTQ files in parallel
echo "======= DECOMPRESS SAMPLES ======="
date

# Number of FASTQ files decompressed simultaneously.
# Twelve is usually safer for disk I/O than launching all 24 at once.
GUNZIP_JOBS=12

find "$READS" \
    -mindepth 2 \
    -maxdepth 2 \
    -type f \
    -name '*.fastq.gz' \
    -print0 |
parallel \
    --null \
    --jobs "$GUNZIP_JOBS" \
    --halt soon,fail=1 \
    --joblog "${BICYCLEPATH}/gunzip.parallel.log" \
    'gunzip -- {}'

echo "======= [✓] SAMPLES DECOMPRESSED ======="
date

# rename paired sample names into their own folder
mv ${READS}/pparfmem001/pparfmem001_1_merged.fastq ${READS}/pparfmem001/pparfmem001_1.fastq  &&
mv ${READS}/pparfmem001/pparfmem001_2_merged.fastq ${READS}/pparfmem001/pparfmem001_2.fastq  &&
echo "[✓] pparfmem001: moved & renamed" &&
mv ${READS}/pparfmem002/pparfmem002_1_merged.fastq ${READS}/pparfmem002/pparfmem002_1.fastq  &&
mv ${READS}/pparfmem002/pparfmem002_2_merged.fastq ${READS}/pparfmem002/pparfmem002_2.fastq  &&
echo "[✓] pparfmem002: moved & renamed" &&
mv ${READS}/pparfmem003/pparfmem003_1_merged.fastq ${READS}/pparfmem003/pparfmem003_1.fastq  &&
mv ${READS}/pparfmem003/pparfmem003_2_merged.fastq ${READS}/pparfmem003/pparfmem003_2.fastq  &&
echo "[✓] pparfmem003: moved & renamed" &&
mv ${READS}/pparpmem004/pparfmem004_1_merged.fastq ${READS}/pparpmem004/pparpmem004_1.fastq  &&
mv ${READS}/pparpmem004/pparfmem004_2_merged.fastq ${READS}/pparpmem004/pparpmem004_2.fastq  &&
echo "[✓] pparpmem004: moved & renamed" &&
mv ${READS}/pparpmem005/pparfmem005_1_merged.fastq ${READS}/pparpmem005/pparpmem005_1.fastq  &&
mv ${READS}/pparpmem005/pparfmem005_2_merged.fastq ${READS}/pparpmem005/pparpmem005_2.fastq  &&
echo "[✓] pparpmem005: movedf& renamed" &&
mv ${READS}/pparpmem006/pparfmem006_1_merged.fastq ${READS}/pparpmem006/pparpmem006_1.fastq  &&
mv ${READS}/pparpmem006/pparfmem006_2_merged.fastq ${READS}/pparpmem006/pparpmem006_2.fastq  &&
echo "[✓] pparpmem006: moved & renamed" &&
mv ${READS}/pparymem007/pparfmem007_1_merged.fastq ${READS}/pparymem007/pparymem007_1.fastq  &&
mv ${READS}/pparymem007/pparfmem007_2_merged.fastq ${READS}/pparymem007/pparymem007_2.fastq  &&
echo "[✓] pparymem007: moved & renamed" &&
mv ${READS}/pparymem008/pparfmem008_1_merged.fastq ${READS}/pparymem008/pparymem008_1.fastq  &&
mv ${READS}/pparymem008/pparfmem008_2_merged.fastq ${READS}/pparymem008/pparymem008_2.fastq  &&
echo "[✓] pparymem008: moved & renamed" &&
mv ${READS}/pparymem009/pparfmem009_1_merged.fastq ${READS}/pparymem009/pparymem009_1.fastq  &&
mv ${READS}/pparymem009/pparfmem009_2_merged.fastq ${READS}/pparymem009/pparymem009_2.fastq  &&
echo "[✓] pparymem009: moved & renamed" &&
mv ${READS}/pparimem010/pparfmem010_1_merged.fastq ${READS}/pparimem010/pparimem010_1.fastq  &&
mv ${READS}/pparimem010/pparfmem010_2_merged.fastq ${READS}/pparimem010/pparimem010_2.fastq  &&
echo "[✓] pparimem010: moved & renamed" &&
mv ${READS}/pparimem011/pparfmem011_1_merged.fastq ${READS}/pparimem011/pparimem011_1.fastq  &&
mv ${READS}/pparimem011/pparfmem011_2_merged.fastq ${READS}/pparimem011/pparimem011_2.fastq  &&
echo "[✓] pparimem011: moved & renamed" &&
mv ${READS}/pparimem012/pparfmem012_1_merged.fastq ${READS}/pparimem012/pparimem012_1.fastq  &&
mv ${READS}/pparimem012/pparfmem012_2_merged.fastq ${READS}/pparimem012/pparimem012_2.fastq  &&
echo "[✓] pparimem012: moved & renamed" &&

echo  "=== [✓] SAMPLES COPIED & RENAMED ==="

# B) BYCYCLE  ~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# Step 1. Create project
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle create-project -p /bicycleAlignment/project -r /bicycleAlignment/ref_genomes -f /bicycleAlignment/reads --paired-mate1-regexp _1.fastq &&

# Step 2. Create the Watson and Crick in-silico bisulfited reference genomes
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle reference-bisulfitation -p /bicycleAlignment/project  &&

# Step 3. Create the reference index
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle reference-index -p /bicycleAlignment/project -t 80 &&

# Step 4. Align to the reference genome
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle align -p /bicycleAlignment/project -t 80 -q2 phred33 &&

# Step 5. Perform methylation analysis and methylcytosine calling
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-methylation -p /bicycleAlignment/project -n 80 -a

 ##---------------------------------------------
# Step 6. Analyze differential methylation
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /bicycleAlignment/project -c pparfmem001,pparfmem002,pparfmem003 -t pparpmem004,pparpmem005,pparpmem006 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /bicycleAlignment/project -c pparfmem001,pparfmem002,pparfmem003 -t pparymem007,pparymem008,pparymem009 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /bicycleAlignment/project -c pparfmem001,pparfmem002,pparfmem003 -t pparimem010,pparimem011,pparimem012 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /bicycleAlignment/project -c pparpmem004,pparpmem005,pparpmem006 -t pparimem010,pparimem011,pparimem012 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /bicycleAlignment/project -c pparpmem004,pparpmem005,pparpmem006 -t pparymem007,pparymem008,pparymem009 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bicycleAlignment:/bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /bicycleAlignment/project -c pparymem007,pparymem008,pparymem009 -t pparimem010,pparimem011,pparimem012 -x CG,CHG,CHH &&

current_date_time="`date "+%Y-%m-%d %H:%M:%S"`" && echo $current_date_time
