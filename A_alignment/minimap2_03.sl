#!/bin/bash
######################################################
#SBATCH -J minimap2_02
#SBATCH --time=00:40:00        # Walltime
#SBATCH --mem-per-cpu=2G     # memory/cpu 
#SBATCH -n 20                   # 1 core means serial
#SBATCH --output=/workspace/cfngle/results-data/minimap2/job_output_minimap_500.txt
#SBATCH --error=/workspace/cfngle/results-data/minimap2/job_output_minimap_500.txt
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
path_sequences=/workspace/cfngle/results-data/sequences/

# filenames
seq_100bp=("AC_CpG_100bp.fasta" "AS_CpG_100bp.fasta" "EH_CpG_100bp.fasta")

seq_200bp=("AC_CpG_200bp.fasta" "AS_CpG_200bp.fasta" "EH_CpG_200bp.fasta")

seq_500bp=("AC_CpG_500bp.fasta" "AS_CpG_500bp.fasta" "EH_CpG_500bp.fasta")

seq_1000bp=("AC_CpG_1000bp.fasta" "AS_CpG_1000bp.fasta" "EH_CpG_1000bp.fasta")

#### MINIMAP2 ####

path_rgenome="AC/GCF_902167405.1_gadMor3.0_genomic.fasta"
## AC rgenome
#AC
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}sequences/${seq_1000bp[0]} > ${path_results}minimap2/AC_AC_1000_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}sequences/${seq_500bp[0]} > ${path_results}minimap2/AC_AC_500_minimap.sam -t 20
#AS
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}sequences/${seq_1000bp[1]} > ${path_results}minimap2/AC_AS_1000_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}sequences/${seq_500bp[1]} > ${path_results}minimap2/AC_AS_500_minimap.sam -t 20
#EH
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}sequences/${seq_1000bp[2]} > ${path_results}minimap2/AC_EH_1000_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AC/GCF_902167405.1_gadMor3.0_genomic.fasta ${path_results}sequences/${seq_500bp[2]} > ${path_results}minimap2/AC_EH_500_minimap.sam -t 20 

path_rgenome="AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta"
## AS rgenome
#AC
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}sequences/${seq_1000bp[0]} > ${path_results}minimap2/AS_AC_1000_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}sequences/${seq_500bp[0]} > ${path_results}minimap2/AS_AC_500_minimap.sam -t 20
#AS
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}sequences/${seq_1000bp[1]} > ${path_results}minimap2/AS_AS_1000_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}sequences/${seq_500bp[1]} > ${path_results}minimap2/AS_AS_500_minimap.sam -t 20
#EH
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}sequences/${seq_1000bp[2]} > ${path_results}minimap2/AS_EH_1000_minimap.sam -t 20 
minimap2 -ax map-ont ${path_raw}AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta ${path_results}sequences/${seq_500bp[2]} > ${path_results}minimap2/AS_EH_500_minimap.sam -t 20 

path_rgenome="EH/fMerMel2.1_cnag1.scaffolds.fa"
## EH rgenome
#AC
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_1000bp[0]} > ${path_results}minimap2/EH_AC_1000_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_500bp[0]} > ${path_results}minimap2/EH_AC_500_minimap.sam -t 20
#AS
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_1000bp[1]} > ${path_results}minimap2/EH_AS_1000_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_500bp[1]} > ${path_results}minimap2/EH_AS_500_minimap.sam -t 20
#EH
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_1000bp[2]} > ${path_results}minimap2/EH_EH_1000_minimap.sam -t 20 
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_500bp[2]} > ${path_results}minimap2/EH_EH_500_minimap.sam -t 20 

## ZF rgenome
path_rgenome=ZF/rgenome/GCF_000002035.6_GRCz11_genomic.fna

#AC
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_1000bp[0]} > ${path_results}minimap2/ZF_AC_1000_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_500bp[0]} > ${path_results}minimap2/ZF_AC_500_minimap.sam -t 20
#AS
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_1000bp[1]} > ${path_results}minimap2/ZF_AS_1000_minimap.sam -t 20
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_500bp[1]} > ${path_results}minimap2/ZF_AS_500_minimap.sam -t 20
#EH
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_1000bp[2]} > ${path_results}minimap2/ZF_EH_1000_minimap.sam -t 20 
minimap2 -ax map-ont ${path_raw}EH/fMerMel2.1_cnag1.scaffolds.fa ${path_results}sequences/${seq_500bp[2]} > ${path_results}minimap2/ZF_EH_500_minimap.sam -t 20 