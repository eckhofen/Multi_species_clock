#!/bin/bash
# Step 08: Aggregate QC across all stages with MultiQC v1.11.

set -euo pipefail
source "$(dirname "$0")/00_config.sh"

module load multiqc   # v1.11

multiqc "${FASTQC_RAW}"      -o "${FASTQC_RAW}"      -i Fastqc-Raw
multiqc "${FASTQC_TRIMMED}"  -o "${FASTQC_TRIMMED}"  -i Fastqc-Trimmed
multiqc "${ALIGNMENTS}"      -o "${MULTIQC_OUT}"     -i QC.alignments
multiqc "${METH_EXTRACT}"    -o "${MULTIQC_OUT}"     -i MethylationExtraction_QC
