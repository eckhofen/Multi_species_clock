#!/bin/bash
######################################################
#SBATCH -J R_methyl_01
#SBATCH --time=01:59:59        # Walltime
#SBATCH --mem=70G     # memory/cpu 
#SBATCH -n 1                   # 1 core means serial
#SBATCH --output=/workspace/cfngle/raw-data/ZF/zzz_methyldata/job_output_bowtie%j.txt
#SBATCH --error=/workspace/cfngle/raw-data/ZF/zzz_methyldata/err_%j.txt
######################################################
###  Run the Parallel Program

module load R 

Rscript /workspace/cfngle/scripts/03c.1_ZF_methylation_extraction_SLURM_v-1.2.R