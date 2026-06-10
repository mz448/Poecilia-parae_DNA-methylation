AUTHOR: mz448
_________________________
GOAL: Run the bicycle alignment pipeline
_________________________
NOTES:

_________________________

```bash title:bicycleAlignment.sh
#!/bin/bash
#SBATCH --job-name=alignment_bicycle
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=80
#SBATCH --mem=900G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH -o /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/Logs/%x_%j.out
#SBATCH -e /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/Logs/%x_%j.err

cd /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/

# _____

# Step 1. Create project  
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle create-project -p /03_bicycleAlignment/project -r /03_bicycleAlignment/ref_genomes -f /03_bicycleAlignment/reads --paired-mate1-regexp _1.fastq &&

# Step 2. Create the Watson and Crick in-silico bisulfited reference genomes 
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle reference-bisulfitation -p /03_bicycleAlignment/project  &&

# Step 3. Create the reference index
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle reference-index -p /03_bicycleAlignment/project -t 80 &&

# Step 4. Align to the reference genome  
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle align -p /03_bicycleAlignment/project -t 80 -q2 phred33 &&

# Step 5. Perform methylation analysis and methylcytosine calling
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-methylation -p /03_bicycleAlignment/project -n 80 -a

 ##---------------------------------------------
# Step 6. Analyze differential methylation
singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /03_bicycleAlignment/project -c pparfmem001,pparfmem002,pparfmem003 -t pparpmem004,pparpmem005,pparpmem006 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /03_bicycleAlignment/project -c pparfmem001,pparfmem002,pparfmem003 -t pparymem007,pparymem008,pparymem009 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /03_bicycleAlignment/project -c pparfmem001,pparfmem002,pparfmem003 -t pparimem010,pparimem011,pparimem012 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /03_bicycleAlignment/project -c pparpmem004,pparpmem005,pparpmem006 -t pparimem010,pparimem011,pparimem012 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /03_bicycleAlignment/project -c pparpmem004,pparpmem005,pparpmem006 -t pparymem007,pparymem008,pparymem009 -x CG,CHG,CHH &&

singularity run --bind /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/03_bicycleAlignment:/03_bicycleAlignment /programs/bicycle-1.8.2/bicycle.sif  bicycle analyze-differential-methylation -p /03_bicycleAlignment/project -c pparymem007,pparymem008,pparymem009 -t pparimem010,pparimem011,pparimem012 -x CG,CHG,CHH && 

current_date_time="`date "+%Y-%m-%d %H:%M:%S"`" && echo $current_date_time

```