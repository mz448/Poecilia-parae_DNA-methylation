#!/bin/bash

# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# AUTHOR: MZF
# DATE: 20260718
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%
# NOTES: 
#   This version has hardcoded paths to preset the run via job manager
# ~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%~~~~~%

set -euo pipefail

SOURCE_DIR="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/01_alignment/bismarkAlignment/bismark/bismark/coverage2cytosine/reports"
OUTPUT_DIR="/local/storage/Projects/ppar_emseq/analysis/011_emseq_pparae_muscle/02_differentialMethylation/01_DMC/01_DNMTools/data"

# Create the output directory if it does not already exist
mkdir -p "$OUTPUT_DIR"

# Copy CpG report files
cp "$SOURCE_DIR"/*.gz "$OUTPUT_DIR"/
  
  # Decompress all copied files
  gunzip "$OUTPUT_DIR"/*.gz