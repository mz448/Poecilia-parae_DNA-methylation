AUTHOR: mz448
_________________________
GOAL: Run **Nextflow Methylseq** pipeline using **Bismark alignment**
_________________________
NOTES:
1) This was done in a singularity container
_________________________

```bash title:bismark_singularityContainer.sh
#!/bin/bash
#SBATCH --job-name=bismark_test
#SBATCH --ntasks=90
#SBATCH --mem=900G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH -o ./Logs/%x_%j.out
#SBATCH -e ./Logs/%x_%j.err

module load nextflow

date

nextflow -v

nextflow run nf-core/methylseq \
-name alignment_bismark \
--input /local/storage/Projects/ppar_emseq/data/009_emseq_pparae_muscle/samples/pparemseqSamples_curated_all_V5.csv \
--bismark_index /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/01_nfcore_methylseq/reference_genome/BismarkIndex \
--fasta /local/storage/Projects/ppar_emseq/data/009_emseq_pparae_muscle/genome/PparFemMitoLambdaPuc19Ver2024.fasta \
--save_reference \
--outdir /local/storage/Projects/ppar_emseq/analysis/009_emseq_pparae_muscle/01_alignment/01_nfcore_methylseq/bismark03 \
--cytosine_report \
--save_align_intermeds \
--unmapped \
--save_trimmed \
--em_seq \
--email mz448@cornell.edu \
--aligner bismark \
--comprehensive \
--max_cpus 90 \
--max_memory 900.GB \
-profile singularity

date
```
