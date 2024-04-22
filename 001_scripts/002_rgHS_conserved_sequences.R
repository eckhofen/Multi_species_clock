#### Overview ####
# Conserved sequences will be extracted from the aligned sequences

#### settings ####
data_folder <- "/powerplant/workspace/cfngle/script_GH/Multi_species_clock/000_data/"
save_folder <- paste0(data_folder, "002_conserved_seq/") # folder where extracted sequences will be saved
suffix <- ".fasta"

#### Preparation ####
# loading libraries
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(GenomicAlignments)
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(dplyr)
library(tidyr)
library(Rsamtools)
library(ggplot2)

#### Loading data ####
# working directory
setwd("/powerplant/workspace/cfngle")

# loading bam files
HS_AC_1000_bt2 <- readGAlignments("results-data/bowtie2/human_AC_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_AS_1000_bt2 <- readGAlignments("results-data/bowtie2/human_AS_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_EH_1000_bt2 <- readGAlignments("results-data/bowtie2/human_EH_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_JM_1000_bt2 <- readGAlignments("results-data/bowtie2/human_JM_243285_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_ZF_1000_bt2 <- readGAlignments("results-data/bowtie2/human_ZF_757883_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))

# adding metadata
AC_metadata <- read.csv("script_GH/Multi_species_clock/000_data/001_sequences/AC_metadata_1000bp.csv")
AC_metadata_matched <- AC_metadata[match(names(HS_AC_1000_bt2), AC_metadata$seq),]
mcols(HS_AC_1000_bt2) <- data.frame(mcols(HS_AC_1000_bt2), AC_metadata_matched$methyl_pos,AC_metadata_matched$methyl_n)
colnames(mcols(HS_AC_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

AS_metadata <- read.csv("script_GH/Multi_species_clock/000_data/001_sequences/AS_metadata_1000bp.csv")
AS_metadata_matched <- AS_metadata[match(names(HS_AS_1000_bt2), AS_metadata$seq),]
mcols(HS_AS_1000_bt2) <- data.frame(mcols(HS_AS_1000_bt2), AS_metadata_matched$methyl_pos,AS_metadata_matched$methyl_n)
colnames(mcols(HS_AS_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

EH_metadata <- read.csv("script_GH/Multi_species_clock/000_data/001_sequences/EH_metadata_1000bp.csv")
EH_metadata_matched <- EH_metadata[match(names(HS_EH_1000_bt2), EH_metadata$seq),]
mcols(HS_EH_1000_bt2) <- data.frame(mcols(HS_EH_1000_bt2), EH_metadata_matched$methyl_pos,EH_metadata_matched$methyl_n)
colnames(mcols(HS_EH_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

JM_metadata <- read.csv("script_GH/Multi_species_clock/000_data/001_sequences/JM_metadata_243285_1000bp.csv")
JM_metadata_matched <- JM_metadata[match(names(HS_JM_1000_bt2), JM_metadata$seq),]
mcols(HS_JM_1000_bt2) <- data.frame(mcols(HS_JM_1000_bt2), JM_metadata_matched$methyl_pos,JM_metadata_matched$methyl_n)
colnames(mcols(HS_JM_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

ZF_metadata <- read.csv("script_GH/Multi_species_clock/000_data/001_sequences/ZF_metadata_757883_1000bp.csv")
ZF_metadata_matched <- ZF_metadata[match(names(HS_ZF_1000_bt2), ZF_metadata$seq),]
mcols(HS_ZF_1000_bt2) <- data.frame(mcols(HS_ZF_1000_bt2), ZF_metadata_matched$methyl_pos,ZF_metadata_matched$methyl_n)
colnames(mcols(HS_ZF_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")


#### Finding overlapping sequences ####
## function to find overlaps from multiple GRanges objects
find.Overlap <- function(...) {
  seq_list <- list(...)
  seqs <- seq_list[[1]]
  for(i in 1:length(seq_list)) {
    overlaps <- findOverlaps(seqs,seq_list[[i]])
    seqs <- seqs[unique(queryHits(overlaps))]
  }
  
  # creating a list of all shared sequences 
  seq_list_overlap <- list()
  
  for(i in 1:length(seq_list)) {
    overlaps <- findOverlaps(seqs,seq_list[[i]])
    seq_list_overlap[[i]] <- seq_list[[i]][unique(subjectHits(overlaps))]
  }
  
  return(seq_list_overlap)
}

#### Run sequences through function ####
# list for aligned sequences
HS_seqs <- list(HS_AC_1000_bt2, HS_AS_1000_bt2, HS_EH_1000_bt2, HS_JM_1000_bt2, HS_ZF_1000_bt2)

# getting overlaps
HS_overlap_seqs <- find.Overlap(HS_seqs[[1]],HS_seqs[[2]],HS_seqs[[3]],HS_seqs[[5]]) # without JM

# ## optional
# HS_overlap_seqs <- find.Overlap(HS_seqs[[1]],HS_seqs[[3]],HS_seqs[[5]]) # testing
# 
# # checking how many seqs in overlaps
# unlist(lapply(HS_overlap_seqs, function(x) length(x)))
# 
# # check how many SMRs
# HS_gr_overlap_seqs <- lapply(HS_overlap_seqs, function(x) granges(x))
# HS_group_gr_overlap <- do.call(c, HS_gr_overlap_seqs)
# length(GenomicRanges::reduce(HS_group_gr_overlap))

# saving overlaps 
save_path <- paste0(save_folder, "HS_AC_AS_EH_ZF_overlaps.Rdata")
save(HS_overlap_seqs, file = save_path)
