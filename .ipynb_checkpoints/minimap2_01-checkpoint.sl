#!/bin/bash
######################################################
#SBATCH -J minimap2
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

#### MINIMAP2 ####

#AC
minimap2 -ax map-ont ${path_raw}ZF/rgenome/GCF_000002035.6_GRCz11_genomic.fna ${path_results}${seq_200bp[0]} > ${path_results}minimap2/ZF_AC_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}ZF/rgenome/GCF_000002035.6_GRCz11_genomic.fna ${path_results}${seq_100bp[0]} > ${path_results}minimap2/ZF_AC_100_minimap.sam -t 20
#AS
minimap2 -ax map-ont ${path_raw}ZF/rgenome/GCF_000002035.6_GRCz11_genomic.fna ${path_results}${seq_200bp[1]} > ${path_results}minimap2/ZF_AS_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}ZF/rgenome/GCF_000002035.6_GRCz11_genomic.fna ${path_results}${seq_100bp[1]} > ${path_results}minimap2/ZF_AS_100_minimap.sam -t 20
#EH
minimap2 -ax map-ont ${path_raw}ZF/rgenome/GCF_000002035.6_GRCz11_genomic.fna ${path_results}${seq_200bp[2]} > ${path_results}minimap2/ZF_EH_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}ZF/rgenome/GCF_000002035.6_GRCz11_genomic.fna ${path_results}${seq_100bp[2]} > ${path_results}minimap2/ZF_EH_100_minimap.sam -t 20