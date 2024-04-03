#!/bin/bash
######################################################
#SBATCH -J cd_hit_01_75
#SBATCH --time=00:40:00        # Walltime
#SBATCH --mem-per-cpu=2G     # memory/cpu 
#SBATCH -n 20                   # 1 core means serial
#SBATCH --output=/workspace/cfngle/results-data/cd-hit/job_output_%j.txt
#SBATCH --error=/workspace/cfngle/results-data/cd-hit/job_output_%j.txt
######################################################
###  Run the Parallel Program

# Alignments for various species to reference genomes

# loading modules
module load bowtie2
module load minimap2 
module load cd-hit

# setting up variables 
# paths
path_raw=/workspace/cfngle/raw-data/
path_results=/workspace/cfngle/results-data/

# filenames
# 200 bp sequences
seq_200bp=("AC_CpG_200bp.fasta" "AS_CpG_200bp.fasta" "EH_CpG_200bp.fasta")

seq_100bp=("AC_CpG_100bp.fasta" "AS_CpG_100bp.fasta" "EH_CpG_100bp.fasta")

#### CD-HIT ####
# clustering sequences with cd-hit
cdhit_opt="-p 1 -c 0.8 -n 5 -T 0 -M 0"
cdhit_nmsfx="_8"

# AC
cd-hit-est-2d -i ${path_results}${seq_100bp[0]} -i2 ${path_results}${seq_100bp[1]} -o ${path_results}cd-hit/AC_AS_100_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
cd-hit-est-2d -i ${path_results}${seq_100bp[0]} -i2 ${path_results}${seq_100bp[2]} -o ${path_results}cd-hit/AC_EH_100_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
# AS
cd-hit-est-2d -i ${path_results}${seq_100bp[1]} -i2 ${path_results}${seq_100bp[2]} -o ${path_results}cd-hit/AS_EH_100_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
cd-hit-est-2d -i ${path_results}${seq_100bp[1]} -i2 ${path_results}${seq_100bp[0]} -o ${path_results}cd-hit/AS_AC_100_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
# EH
cd-hit-est-2d -i ${path_results}${seq_100bp[2]} -i2 ${path_results}${seq_100bp[0]} -o ${path_results}cd-hit/EH_AC_100_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
cd-hit-est-2d -i ${path_results}${seq_100bp[2]} -i2 ${path_results}${seq_100bp[1]} -o ${path_results}cd-hit/EH_AS_100_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt