#!/bin/bash
######################################################
#SBATCH -J MAFFT_03
#SBATCH --time=01:00:00        # Walltime
#SBATCH --mem-per-cpu=2G     # memory/cpu 
#SBATCH -n 40                   # 1 core means serial
#SBATCH --output=/workspace/cfngle/results-data/mafft/job_output_mafft%j.txt
#SBATCH --error=/workspace/cfngle/results-data/mafft/job_output_mafft%j.txt
######################################################
###  Run the Parallel Program

#### MAFFT ####

## modules
module load mafft

## setting up variables 

# paths
path_raw=/workspace/cfngle/raw-data/
path_results=/workspace/cfngle/results-data/
path_sequences=/workspace/cfngle/results-data/sequences/

# arguments
mafft_args="--auto --thread -1"

mafft $mafft_args ${path_results}02_conserved_seq/AC_AS_EH_1000_conserved.fasta > ${path_results}02_conserved_seq/AC_AS_EH_1000_conserved_MSA_mafft.fasta
#mafft $mafft_args ${path_sequences}AC_AS_EH_CpG_100bp.fasta > ${path_results}AC_AS_EH_CpG_100bp.fasta 
#mafft $mafft_args ${path_sequences}AC_AS_EH_CpG_200bp.fasta > ${path_results}AC_AS_EH_CpG_200bp.fasta 
#mafft $mafft_args ${path_sequences}AC_AS_EH_CpG_500bp.fasta > ${path_results}AC_AS_EH_CpG_500bp.fasta 
#mafft $mafft_args ${path_sequences}AC_AS_EH_CpG_1000bp.fasta > ${path_results}AC_AS_EH_CpG_1000bp.fasta 