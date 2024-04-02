#!/bin/bash
######################################################
#SBATCH -J minimap2_02
#SBATCH --time=00:10:00        # Walltime
#SBATCH --mem-per-cpu=2G     # memory/cpu 
#SBATCH -n 20                   # 1 core means serial
#SBATCH --output=/workspace/cfngle/results-data/minimap2/job_output_minimap_500.txt
#SBATCH --error=/workspace/cfngle/results-data/minimap2/job_output_minimap_500.txt
######################################################
###  Run the Parallel Program

# Alignments for various species to reference genomes

# loading modules
module load minimap2 


# setting up variables 
# paths
path_raw=/workspace/cfngle/raw-data/
path_results=/workspace/cfngle/results-data/
path_sequences=/workspace/cfngle/results-data/sequences/

# filenames
seq_1000bp=("JM_CpG_1000bp.fasta")

#### MINIMAP2 ####

#ZF
echo "ZF"
#JM
minimap2 -ax map-ont ${path_raw}ZF/rgenome/GCF_000002035.6_GRCz11_genomic.fna ${path_sequences}${seq_1000bp[0]} > ${path_results}minimap2/ZF_JM_1000_minimap.sam

path_rgenome="AC/GCF_902167405.1_gadMor3.0_genomic.fasta"
#AC
echo "AC"
#JM
minimap2 -ax map-ont ${path_raw}$path_rgenome ${path_sequences}${seq_1000bp[0]} > ${path_results}minimap2/AC_JM_1000_minimap.sam

path_rgenome="AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta"
#AS
echo "AS"
#JM
minimap2 -ax map-ont ${path_raw}$path_rgenome ${path_sequences}${seq_1000bp[0]} > ${path_results}minimap2/AS_JM_1000_minimap.sam

path_rgenome="EH/fMerMel2.1_cnag1.scaffolds.fa"
#EH
echo "EH"
#JM
minimap2 -ax map-ont ${path_raw}$path_rgenome ${path_sequences}${seq_1000bp[0]} > ${path_results}minimap2/EH_JM_1000_minimap.sam

path_rgenome="JM/rgenome/GCF_002234675.1_ASM223467v1_genomic.fna"
#JM
echo "JM"
#JM
minimap2 -ax map-ont ${path_raw}$path_rgenome ${path_sequences}${seq_1000bp[0]} > ${path_results}minimap2/JM_JM_1000_minimap.sam
