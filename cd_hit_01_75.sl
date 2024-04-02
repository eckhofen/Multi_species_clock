#!/bin/bash
######################################################
#SBATCH -J cd_hit_01_75
#SBATCH --time=00:40:00        # Walltime
#SBATCH --mem-per-cpu=2G     # memory/cpu 
#SBATCH -n 20                   # 1 core means serial
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
cdhit_opt="-p 1 -c 0.75 -n 4 -T 0 -M 0"
cdhit_nmsfx="_75"

# logging
script ${path_results}cd-hit/log${cdhit_nmsfx}.log

# AC
cd-hit-est-2d -i ${path_results}${seq_200bp[0]} -i2 ${path_results}${seq_200bp[1]} -o ${path_results}cd-hit/AC_AS_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
cd-hit-est-2d -i ${path_results}${seq_200bp[0]} -i2 ${path_results}${seq_200bp[2]} -o ${path_results}cd-hit/AC_EH_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
# AS
cd-hit-est-2d -i ${path_results}${seq_200bp[1]} -i2 ${path_results}${seq_200bp[2]} -o ${path_results}cd-hit/AS_EH_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
cd-hit-est-2d -i ${path_results}${seq_200bp[1]} -i2 ${path_results}${seq_200bp[0]} -o ${path_results}cd-hit/AS_AC_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
# EH
cd-hit-est-2d -i ${path_results}${seq_200bp[2]} -i2 ${path_results}${seq_200bp[0]} -o ${path_results}cd-hit/EH_AC_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt
cd-hit-est-2d -i ${path_results}${seq_200bp[2]} -i2 ${path_results}${seq_200bp[1]} -o ${path_results}cd-hit/EH_AS_cd-hit${cdhit_nmsfx}.fasta $cdhit_opt

exit