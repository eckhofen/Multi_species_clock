#!/bin/bash
# Step 07: Bismark methylation extraction.

set -euo pipefail
source "$(dirname "$0")/00_config.sh"

module load bismark/0.23.0

FILES=$(find "${ALIGNMENTS}/" -type f -name "*.bam")

for file in ${FILES}; do
    COMMAND="bismark_methylation_extractor \
        ${file} \
        -s \
        --merge_non_CpG \
        --cytosine_report \
        --scaffolds \
        --genome_folder ${GENOME} \
        -o ${METH_EXTRACT}"
    echo "${COMMAND}"
done | abatch -j "${METH_EXTRACT}/meth-extract" --time 23:59:00 --mem 10G | sbatch
