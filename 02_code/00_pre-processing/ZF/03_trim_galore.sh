#!/bin/bash
# Step 03: Trim Galore on ZF RRBS reads.

set -euo pipefail
source "$(dirname "$0")/00_config.sh"

module load trim_galore   # v0.4.1

FILES=$(basename -a "${RAW_READS}"/*.fastq)

for FILE in ${FILES}; do
    IN_FILE="${RAW_READS}/${FILE}"
    COMMAND="trim_galore \
        --rrbs \
        --quality 25 \
        --clip_R1 4 \
        --three_prime_clip_R1 2 \
        ${IN_FILE} \
        --output_dir ${TRIM_GALORE}"
    echo "${COMMAND}"
done | abatch -j "${TRIM_GALORE}/trim_galore-logs" --time 04:00:00 --mem 1G | sbatch
