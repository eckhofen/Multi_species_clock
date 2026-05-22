#!/bin/bash
# Step 06: Align trimmed ZF reads to GRCz11 with Bismark.

set -euo pipefail
source "$(dirname "$0")/00_config.sh"

LOG="${ALIGNMENTS}/logs"
mkdir -p "${LOG}"

module load bismark/0.23.0

FILES=$(ls "${TRIM_GALORE}"/*.fq)

for file in ${FILES}; do
    NAME=$(basename "${file}")
    PREFIX=$(echo "${NAME}" | awk -F'[_]' '{print $1"_"$2"_Alignments"}')
    COMMAND="bismark \
        --genome ${GENOME} \
        ${file} \
        --o ${ALIGNMENTS}/${PREFIX} \
        --multicore 5 \
        --non_directional \
        --local"
    echo "${COMMAND}"
done | abatch -j "${ALIGNMENTS}/bis-align-local-ZF" --time 2-23:59:00 --mem 50G -c 1 | sbatch
