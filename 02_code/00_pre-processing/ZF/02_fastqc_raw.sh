#!/bin/bash
# Step 02: FastQC on raw ZF RRBS reads + MultiQC summary.
# Methods: "Quality control was performed on both raw and trimmed reads using
# FastQC and assessed using MultiQC v1.11".

set -euo pipefail
source "$(dirname "$0")/00_config.sh"

module load FastQC
module load multiqc   # v1.11

FILES=$(find "${RAW_READS}" -maxdepth 1 -type f -name "*.fastq")

for file in ${FILES}; do
    COMMAND="fastqc --nogroup -q -t 2 -o ${FASTQC_RAW} ${file}"
    echo "${COMMAND}"
done | abatch -j "${FASTQC_RAW}/fastqc_raw" --time 01:00:00 --mem 20G -c 10 | sbatch

# After all FastQC jobs finish, aggregate with MultiQC:
#   multiqc "${FASTQC_RAW}" -o "${FASTQC_RAW}" -i Fastqc-Raw
