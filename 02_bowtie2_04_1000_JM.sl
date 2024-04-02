#!/bin/bash
######################################################
#SBATCH -J bowtie_1000_JM
#SBATCH --time=02:00:00        # Walltime
#SBATCH --mem-per-cpu=2G     # memory/cpu 
#SBATCH -n 20                   # 1 core means serial
#SBATCH --output=/workspace/cfngle/results-data/bowtie2/job_output_bowtie%j.txt
#SBATCH --error=/workspace/cfngle/results-data/bowtie2/job_output_bowtie%j.txt
######################################################
###  Run the Parallel Program

## BOWTIE2 ####

## modules
module load bowtie

## setting up variables 
# paths
path_raw=/workspace/cfngle/raw-data/
path_results=/workspace/cfngle/results-data/bowtie2/
path_sequences=/workspace/cfngle/results-data/sequences/

# filenames
seq_100bp=("AC_CpG_100bp" "AS_CpG_100bp" "EH_CpG_100bp")
seq_200bp=("AC_CpG_200bp" "AS_CpG_200bp" "EH_CpG_200bp")
seq_500bp=("AC_CpG_500bp" "AS_CpG_500bp" "EH_CpG_500bp")
seq_1000bp=("AC_CpG_1000bp" "AS_CpG_1000bp" "EH_CpG_1000bp")

# file ending
suffix=".fasta"

# arguments
bowtie2_args="--very-sensitive --local -p 8"


bowtie2-build ${path_raw}JM/rgenome/GCF_002234675.1_ASM223467v1_genomic.fna ${path_raw}JM/rgenome/bowtie2-index/JM_bowtie2
echo "bowtie2 indexing for JM done"

## rgenome JM
echo "JM"
for (( i=0; i<${#seq_1000bp[@]}; i++ )); 
do 
    bowtie2 $bowtie2_args -x ${path_raw}JM/rgenome/bowtie2-index/JM_bowtie2 -f -U ${path_sequences}${seq_1000bp[$i]}$suffix -S ${path_results}JM_${seq_1000bp[$i]}_bt2_.sam -N 1
done