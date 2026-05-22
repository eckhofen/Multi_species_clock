#!/bin/bash
# Step 05: Prepare the GRCz11 reference for Bismark 

set -euo pipefail
source "$(dirname "$0")/00_config.sh"

module load bismark/0.23.0

# Place the GRCz11 FASTA in ${GENOME} before running (e.g.
# GCF_000002035.6_GRCz11_genomic.fna). Standard parameters per Methods.
bismark_genome_preparation "${GENOME}"
