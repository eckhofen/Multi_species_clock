#### Overview ####
# Conserved sequences will be extracted from the aligned sequences

#### Settings ####
data_folder <- "/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/Multi_species_clock/000_data/"
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
setwd("/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/")

# loading bam files
EH_AC_1000_bt2 <- readGAlignments("Temp/EH_AC_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
EH_AS_1000_bt2 <- readGAlignments("Temp/EH_AS_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
EH_EH_1000_bt2 <- readGAlignments("Temp/EH_EH_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
EH_JM_1000_bt2 <- readGAlignments("Temp/EH_JM_243285_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
EH_ZF_1000_bt2 <- readGAlignments("Temp/EH_ZF_757883_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))

# adding metadata
setwd("/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/Multi_species_clock/")

AC_metadata <- read.csv("000_data/001_sequences/AC_metadata_1000bp.csv")
AC_metadata_matched <- AC_metadata[match(names(EH_AC_1000_bt2), AC_metadata$seq),]
mcols(EH_AC_1000_bt2) <- data.frame(mcols(EH_AC_1000_bt2), AC_metadata_matched$methyl_pos,AC_metadata_matched$methyl_n)
colnames(mcols(EH_AC_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

AS_metadata <- read.csv("000_data/001_sequences/AS_metadata_1000bp.csv")
AS_metadata_matched <- AS_metadata[match(names(EH_AS_1000_bt2), AS_metadata$seq),]
mcols(EH_AS_1000_bt2) <- data.frame(mcols(EH_AS_1000_bt2), AS_metadata_matched$methyl_pos,AS_metadata_matched$methyl_n)
colnames(mcols(EH_AS_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

EH_metadata <- read.csv("000_data/001_sequences/EH_metadata_1000bp.csv")
EH_metadata_matched <- EH_metadata[match(names(EH_EH_1000_bt2), EH_metadata$seq),]
mcols(EH_EH_1000_bt2) <- data.frame(mcols(EH_EH_1000_bt2), EH_metadata_matched$methyl_pos,EH_metadata_matched$methyl_n)
colnames(mcols(EH_EH_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

JM_metadata <- read.csv("000_data/001_sequences/JM_metadata_243285_1000bp.csv")
JM_metadata_matched <- JM_metadata[match(names(EH_JM_1000_bt2), JM_metadata$seq),]
mcols(EH_JM_1000_bt2) <- data.frame(mcols(EH_JM_1000_bt2), JM_metadata_matched$methyl_pos,JM_metadata_matched$methyl_n)
colnames(mcols(EH_JM_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

ZF_metadata <- read.csv("000_data/001_sequences/ZF_metadata_757883_1000bp.csv")
ZF_metadata_matched <- ZF_metadata[match(names(EH_ZF_1000_bt2), ZF_metadata$seq),]
mcols(EH_ZF_1000_bt2) <- data.frame(mcols(EH_ZF_1000_bt2), ZF_metadata_matched$methyl_pos,ZF_metadata_matched$methyl_n)
colnames(mcols(EH_ZF_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")


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
EH_seqs <- list(EH_AC_1000_bt2, EH_AS_1000_bt2, EH_EH_1000_bt2, EH_JM_1000_bt2, EH_ZF_1000_bt2)

# getting overlaps
EH_overlap_seqs <- find.Overlap(EH_seqs[[1]],EH_seqs[[2]],EH_seqs[[3]],EH_seqs[[4]],EH_seqs[[5]])
EH_overlap_seqs <- find.Overlap(EH_seqs[[1]],EH_seqs[[2]],EH_seqs[[3]],EH_seqs[[5]]) # without JM
EH_overlap_seqs

## optional
EH_overlap_seqs <- find.Overlap(EH_seqs[[1]],EH_seqs[[2]],EH_seqs[[3]]) # testing

# checking how many seqs in overlaps
unlist(lapply(EH_overlap_seqs, function(x) length(x)))

# check how many SMRs
EH_gr_overlap_seqs <- lapply(EH_overlap_seqs, function(x) granges(x))
EH_group_gr_overlap <- do.call(c, EH_gr_overlap_seqs)
length(GenomicRanges::reduce(EH_group_gr_overlap))

# saving overlaps 
save_path <- paste0(save_folder, "HS_AC_AS_EH_ZF_overlaps.Rdata")
save(HS_overlap_seqs, file = save_path)
