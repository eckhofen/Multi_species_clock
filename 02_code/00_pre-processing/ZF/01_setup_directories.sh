#!/bin/bash
# Step 01: Create the project directory tree for the ZF RRBS pipeline.

set -euo pipefail
source "$(dirname "$0")/00_config.sh"

mkdir -p "${RAW_READS}" \
         "${FASTQC_RAW}" \
         "${TRIM_GALORE}" \
         "${FASTQC_TRIMMED}" \
         "${GENOME}" \
         "${ALIGNMENTS}" \
         "${METH_EXTRACT}" \
         "${MULTIQC_OUT}"

echo "Directory tree ready under ${PROJECT}"
