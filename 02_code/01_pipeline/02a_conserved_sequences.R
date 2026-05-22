# Metadata ----------------------------------------------------------------
# Project: Shared Methylation Regions
# Description: Extracting conserved sequences from aligned sequences
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026-05

#### Settings ####
save_folder <- "01_data/02_external_server_outputs/02_conserved_seq/"
dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)
suffix <- ".fasta"

#### Preparation ####
# loading libraries
library(GenomicRanges)
library(GenomicAlignments)
library(Biostrings)
library(tidyverse)
library(ggVennDiagram)
library(Rsamtools)

#### Loading data ####
# loading color palettes
load("01_data/04_metadata/color_palettes.RData")

# loading bam files
HS_AC_1000_bt2 <- readGAlignments("01_data/02_external_server_outputs/02_conserved_seq/HS_AC_CpG_1000bp.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_AS_1000_bt2 <- readGAlignments("01_data/02_external_server_outputs/02_conserved_seq/HS_AS_CpG_1000bp.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_EH_1000_bt2 <- readGAlignments("01_data/02_external_server_outputs/02_conserved_seq/HS_EH_CpG_1000bp.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_ZF_1000_bt2 <- readGAlignments("01_data/02_external_server_outputs/02_conserved_seq/HS_ZF_CpG_1000bp.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))

# adding metadata
AC_metadata <- read.csv("01_data/02_external_server_outputs/01_sequences/AC_metadata_1000bp.csv")
AC_metadata_matched <- AC_metadata[match(names(HS_AC_1000_bt2), AC_metadata$seq),]
mcols(HS_AC_1000_bt2) <- data.frame(mcols(HS_AC_1000_bt2), AC_metadata_matched$methyl_pos,AC_metadata_matched$methyl_n)
colnames(mcols(HS_AC_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

AS_metadata <- read.csv("01_data/02_external_server_outputs/01_sequences/AS_metadata_1000bp.csv")
AS_metadata_matched <- AS_metadata[match(names(HS_AS_1000_bt2), AS_metadata$seq),]
mcols(HS_AS_1000_bt2) <- data.frame(mcols(HS_AS_1000_bt2), AS_metadata_matched$methyl_pos,AS_metadata_matched$methyl_n)
colnames(mcols(HS_AS_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

EH_metadata <- read.csv("01_data/02_external_server_outputs/01_sequences/EH_metadata_1000bp.csv")
EH_metadata_matched <- EH_metadata[match(names(HS_EH_1000_bt2), EH_metadata$seq),]
mcols(HS_EH_1000_bt2) <- data.frame(mcols(HS_EH_1000_bt2), EH_metadata_matched$methyl_pos,EH_metadata_matched$methyl_n)
colnames(mcols(HS_EH_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

ZF_metadata <- read.csv("01_data/02_external_server_outputs/01_sequences/ZF_metadata_757883_1000bp.csv")
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
HS_seqs <- list(HS_AC_1000_bt2, HS_AS_1000_bt2, HS_EH_1000_bt2, HS_ZF_1000_bt2)

# getting overlaps
HS_overlap_seqs <- find.Overlap(HS_seqs[[1]],HS_seqs[[2]],HS_seqs[[3]],HS_seqs[[4]])

## optional
# check how many SMRs
HS_gr_overlap_seqs <- lapply(HS_overlap_seqs, function(x) granges(x))
HS_group_gr_overlap <- do.call(c, HS_gr_overlap_seqs)
length(GenomicRanges::reduce(HS_group_gr_overlap))

# saving overlaps 
save_path <- paste0(save_folder, "HS_AC_AS_EH_ZF_overlaps_bt2.Rdata")
save(HS_overlap_seqs, file = save_path)