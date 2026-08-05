#!/bin/bash
#SBATCH --job-name=methylseq_bismark
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=90
#SBATCH --mem=900G
#SBATCH --partition=regular
#SBATCH --qos=regular
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mz448@cornell.edu
#SBATCH -o ./logs/%x_%j.out
#SBATCH -e ./logs/%x_%j.err

set -euo pipefail


# ==============================================================================
# Project paths
# ==============================================================================

PROJECT="/local/storage/Projects/ppar_emseq"
RUN_DIR="${PROJECT}/analysis/011_emseq_pparae_muscle/01_alignment/bismarkAlignment"

INPUT="${PROJECT}/data/011_emseq_pparae_muscle/allCuratedSamples/all_samples_V01.csv"

FASTA="${PROJECT}/data/011_emseq_pparae_muscle/genome/PparFemMitoLambdaPuc19Ver2024.fasta"

OUTDIR="${RUN_DIR}/bismark"
WORK_DIR="${RUN_DIR}/nextflow_work"

NXF_HOME_DIR="${RUN_DIR}/.nextflow"
NXF_CONTAINER_DIR="${RUN_DIR}/containers/nfcore_singularity"
SINGULARITY_CACHE="${RUN_DIR}/containers/singularity_cache"
APPTAINER_CACHE="${RUN_DIR}/containers/apptainer_cache"
CONTAINER_TMP="${RUN_DIR}/containers/tmp"

mkdir -p \
    "${OUTDIR}" \
    "${WORK_DIR}" \
    "${NXF_HOME_DIR}" \
    "${NXF_CONTAINER_DIR}" \
    "${SINGULARITY_CACHE}" \
    "${APPTAINER_CACHE}" \
    "${CONTAINER_TMP}"

cd "${RUN_DIR}"


# ==============================================================================
# Software available on the host
# ==============================================================================

module load nextflow

echo "Run started: $(date)"
echo "Host:        $(hostname)"
echo "Nextflow:    $(command -v nextflow)"
echo "Singularity: $(command -v singularity)"

nextflow -version
singularity --version

# ==============================================================================
# Singularity bind configuration
# ==============================================================================

# The FASTQ files, genome, output directory, work directory and container cache
# are all located below this project directory.
export APPTAINER_BINDPATH="${PROJECT}"
export SINGULARITY_BINDPATH="${PROJECT}"



# ==============================================================================
# Nextflow and container caches
# ==============================================================================

export NXF_HOME="${NXF_HOME_DIR}"
export NXF_SINGULARITY_CACHEDIR="${NXF_CONTAINER_DIR}"

export SINGULARITY_CACHEDIR="${SINGULARITY_CACHE}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHE}"

export SINGULARITY_TMPDIR="${CONTAINER_TMP}"
export APPTAINER_TMPDIR="${CONTAINER_TMP}"

# Prevent host Python packages from being imported inside the containers.
unset PYTHONPATH

echo "APPTAINER_BINDPATH=${APPTAINER_BINDPATH}"
echo "NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR}"
echo "WORK_DIR=${WORK_DIR}"
echo "OUTDIR=${OUTDIR}"

# ==============================================================================
# Validate required input files
# ==============================================================================

if [[ ! -r "${INPUT}" ]]; then
    echo "ERROR: Cannot read samplesheet: ${INPUT}" >&2
    exit 1
fi

if [[ ! -r "${FASTA}" ]]; then
    echo "ERROR: Cannot read reference genome: ${FASTA}" >&2
    exit 1
fi


# ==============================================================================
# Run nf-core/methylseq with Singularity
# ==============================================================================



nextflow run nf-core/methylseq \
    -r 4.2.0 \
    -name bismark_03 \
    -profile singularity \
    --input "${INPUT}" \
    --fasta "${FASTA}" \
    --outdir "${OUTDIR}" \
    --aligner bismark \
    --save_reference \
    --cytosine_report \
    --save_align_intermeds \
    --unmapped \
    --save_trimmed \
    --em_seq \
    --email mz448@cornell.edu \
    --comprehensive \
    --max_cpus 90 \
    --max_memory 900.GB

date