#!/bin/bash
######################################################
#SBATCH -J bowtie_index_human
#SBATCH --time=010:00:00        # Walltime
#SBATCH --mem=10G     # memory/cpu 
#SBATCH -n 1                   # 1 core means serial
#SBATCH --output=/workspace/cfngle/results-data/bowtie2/job_output_bowtie%j.txt
#SBATCH --error=/workspace/cfngle/results-data/bowtie2/job_output_bowtie%j.txt
######################################################
###  Run the Parallel Program

## BOWTIE2 ####

## modules
module load bowtie2

## setting up variables 
# paths
path_raw=/workspace/cfngle/raw-data/
path_results=/workspace/cfngle/results-data/bowtie2/

## indexing
bowtie2-build ${path_raw}zzz_human_rgenome/GRCh38.fa ${path_raw}zzz_human_rgenome/bowtie2-index/human_bowtie2
