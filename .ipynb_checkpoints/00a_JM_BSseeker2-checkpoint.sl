#!/bin/bash

#SBATCH --job-name=indexBS
#SBATCH --time=24:00:00
#SBATCH --mem=20G
#SBATCH --output=/workspace/cfngle/raw-data/JM/rgenome/output.log

PROJECT="/workspace/cfngle/raw-data/JM"
GENOME="/workspace/cfngle/raw-data/JM/rgenome"
OUT=$GENOME

module load conda
conda deactivate
conda activate cfngle_env-01

bs_seeker2-build.py \
    -f ${GENOME}/GCF_002234675.1_ASM223467v1_genomic.fasta \
    --aligner bowtie2 \
    -r \
    -d /workspace/cfngle/raw-data/JM/rgenome/BSseeker2-index