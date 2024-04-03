#!/bin/bash
######################################################
#SBATCH -J bowtie_1000_all
#SBATCH --time=05:00:00        # Walltime
#SBATCH --mem=20G     # memory/cpu 
#SBATCH -n 4                   # 1 core means serial
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
path_sequences=/workspace/cfngle/results-data/sequences/

# filenames
seq_1000bp=("JM_243285_CpG_1000bp" "ZF_757883_CpG_1000bp")

# file ending
suffix=".fasta"

# arguments
bowtie2_args="--very-sensitive --local -p 8"


## rgenome Human
echo "human reference genome"
for (( i=0; i<${#seq_1000bp[@]}; i++ )); 
do 
   bowtie2 $bowtie2_args -x /workspace/cfngle/raw-data/zzz_human_rgenome/bowtie2-index/human_bowtie2 -f -U ${path_sequences}${seq_1000bp[$i]}$suffix -S ${path_results}human_${seq_1000bp[$i]}_bt2.sam -N 1
done
