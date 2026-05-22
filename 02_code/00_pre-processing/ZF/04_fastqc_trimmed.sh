#!/bin/bash
# Step 04: FastQC on trimmed ZF reads + MultiQC summary.

set -euo pipefail
source "$(dirname "$0")/00_config.sh"

module load FastQC
module load multiqc   # v1.11

FILES=$(ls "${TRIM_GALORE}"/*.fq)

for file in ${FILES}; do
    COMMAND="fastqc --nogroup -q -t 2 -o ${FASTQC_TRIMMED} ${file}"
    echo "${COMMAND}"
done | abatch -j "${FASTQC_TRIMMED}/fastqc-trimmed-logs" --time 02:00:00 --mem 1G | sbatch

# After all FastQC jobs finish, aggregate with MultiQC:
#   multiqc "${FASTQC_TRIMMED}" -o "${FASTQC_TRIMMED}" -i Fastqc-Trimmed
