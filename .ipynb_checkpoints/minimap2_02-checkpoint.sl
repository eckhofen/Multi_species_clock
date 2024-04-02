#!/bin/bash
######################################################
#SBATCH -J minimap2_02
#SBATCH --time=00:40:00        # Walltime
#SBATCH --mem-per-cpu=2G     # memory/cpu 
#SBATCH -n 20                   # 1 core means serial
#SBATCH --output=/workspace/cfngle/results-data/minimap2/job_output_%j.txt
#SBATCH --error=/workspace/cfngle/results-data/minimap2/job_output_%j.txt
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

path_rgenome="AC/GCF_902167405.1_gadMor3.0_genomic.fasta"
#AC
#AC
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}${seq_200bp[0]} > ${path_results}minimap2/AC_AC_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}${seq_100bp[0]} > ${path_results}minimap2/AC_AC_100_minimap.sam -t 20
#AS
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}${seq_200bp[1]} > ${path_results}minimap2/AC_AS_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}${seq_100bp[1]} > ${path_results}minimap2/AC_AS_100_minimap.sam -t 20
#EH
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}${seq_200bp[2]} > ${path_results}minimap2/AC_EH_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}${seq_100bp[2]} > ${path_results}minimap2/AC_EH_100_minimap.sam -t 20 

path_rgenome="AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta"
#AS
#AC
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}${seq_200bp[0]} > ${path_results}minimap2/AS_AC_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}${seq_100bp[0]} > ${path_results}minimap2/AS_AC_100_minimap.sam -t 20
#AS
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}${seq_200bp[1]} > ${path_results}minimap2/AS_AS_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}${seq_100bp[1]} > ${path_results}minimap2/AS_AS_100_minimap.sam -t 20
#EH
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}${seq_200bp[2]} > ${path_results}minimap2/AS_EH_200_minimap.sam -t 20 
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}${seq_100bp[2]} > ${path_results}minimap2/AS_EH_100_minimap.sam -t 20 

path_rgenome="EH/fMerMel2.1_cnag1.scaffolds.fa"
#EH\
#AC
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}${seq_200bp[0]} > ${path_results}minimap2/EH_AC_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}${seq_100bp[0]} > ${path_results}minimap2/EH_AC_100_minimap.sam -t 20
#AS
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}${seq_200bp[1]} > ${path_results}minimap2/EH_AS_200_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}${seq_100bp[1]} > ${path_results}minimap2/EH_AS_100_minimap.sam -t 20
#EH
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}${seq_200bp[2]} > ${path_results}minimap2/EH_EH_200_minimap.sam -t 20 
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}${seq_100bp[2]} > ${path_results}minimap2/EH_EH_100_minimap.sam -t 20 
