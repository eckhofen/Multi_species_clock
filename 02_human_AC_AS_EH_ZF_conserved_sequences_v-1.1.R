#### Overview ####
# Conserved sequences will be extracted from the aligned sequences

#### Preparation ####
# loading libraries
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(GenomicAlignments)
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
# library(ggbio) # https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
library(dplyr)
library(tidyr)
library(Rsamtools)
library(ggplot2)
#require(BiocManager)


#### loading data ####
setwd("/powerplant/workspace/cfngle")

# defining objects 
save_path <- "/workspace/cfngle/results-data/02_conserved_seq/"
suffix <- ".fasta"

# data bowtie2 
HS_AC_1000_bt2 <- readGAlignments("results-data/bowtie2/human_AC_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_AS_1000_bt2 <- readGAlignments("results-data/bowtie2/human_AS_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_EH_1000_bt2 <- readGAlignments("results-data/bowtie2/human_EH_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_JM_1000_bt2 <- readGAlignments("results-data/bowtie2/human_JM_243285_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
HS_ZF_1000_bt2 <- readGAlignments("results-data/bowtie2/human_ZF_757883_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))


## BOWTIE2 ##JM
AC_metadata <- read.csv("results-data/sequences/AC_metadata_1000bp.csv")
AC_metadata_matched <- AC_metadata[match(names(HS_AC_1000_bt2), AC_metadata$seq),]
mcols(HS_AC_1000_bt2) <- data.frame(mcols(HS_AC_1000_bt2), AC_metadata_matched$methyl_pos,AC_metadata_matched$methyl_n)
colnames(mcols(HS_AC_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

AS_metadata <- read.csv("results-data/sequences/AS_metadata_1000bp.csv")
AS_metadata_matched <- AS_metadata[match(names(HS_AS_1000_bt2), AS_metadata$seq),]
mcols(HS_AS_1000_bt2) <- data.frame(mcols(HS_AS_1000_bt2), AS_metadata_matched$methyl_pos,AS_metadata_matched$methyl_n)
colnames(mcols(HS_AS_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

EH_metadata <- read.csv("results-data/sequences/EH_metadata_1000bp.csv")
EH_metadata_matched <- EH_metadata[match(names(HS_EH_1000_bt2), EH_metadata$seq),]
mcols(HS_EH_1000_bt2) <- data.frame(mcols(HS_EH_1000_bt2), EH_metadata_matched$methyl_pos,EH_metadata_matched$methyl_n)
colnames(mcols(HS_EH_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

JM_metadata <- read.csv("results-data/sequences/JM_metadata_243285_1000bp.csv")
JM_metadata_matched <- JM_metadata[match(names(HS_JM_1000_bt2), JM_metadata$seq),]
mcols(HS_JM_1000_bt2) <- data.frame(mcols(HS_JM_1000_bt2), JM_metadata_matched$methyl_pos,JM_metadata_matched$methyl_n)
colnames(mcols(HS_JM_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

ZF_metadata <- read.csv("results-data/sequences/ZF_metadata_7578831000bp.csv")
ZF_metadata_matched <- ZF_metadata[match(names(HS_ZF_1000_bt2), ZF_metadata$seq),]
mcols(HS_ZF_1000_bt2) <- data.frame(mcols(HS_ZF_1000_bt2), ZF_metadata_matched$methyl_pos,ZF_metadata_matched$methyl_n)
colnames(mcols(HS_ZF_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

# This is just checking which chromosomes/scaffolds/contigs are shared between all aligned seqs
# shared_AC_1000_mini <- AC_AC_1000_mini[seqnames(AC_AC_1000_mini) %in% seqnames(AC_AS_1000_mini) & seqnames(AC_AC_1000_mini) %in% seqnames(AC_EH_1000_mini)]

#### Finding overlapping sequences ####
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

## Bowtie2
HS_overlap_seqs_bt2 <- find.Overlap(HS_AC_1000_bt2,
                                 HS_AS_1000_bt2,
                                 HS_EH_1000_bt2, 
                                 HS_JM_1000_bt2,
                                 HS_ZF_1000_bt2
                                 )
save(HS_overlap_seqs_bt2, file = "results-data/02_conserved_seq/HS_AC_EH_ZF_overlaps_bt2.Rdata")

