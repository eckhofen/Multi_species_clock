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
AC_AC_1000_bt2 <- readGAlignments("results-data/bowtie2/AC_AC_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AC_AS_1000_bt2 <- readGAlignments("results-data/bowtie2/AC_AS_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AC_EH_1000_bt2 <- readGAlignments("results-data/bowtie2/AC_EH_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AC_JM_1000_bt2 <- readGAlignments("results-data/bowtie2/AC_JM_243285_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AC_ZF_1000_bt2 <- readGAlignments("results-data/bowtie2/AC_ZF_757883_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))

AS_AC_1000_bt2 <- readGAlignments("results-data/bowtie2/AS_AC_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AS_AS_1000_bt2 <- readGAlignments("results-data/bowtie2/AS_AS_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AS_EH_1000_bt2 <- readGAlignments("results-data/bowtie2/AS_EH_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AS_JM_1000_bt2 <- readGAlignments("results-data/bowtie2/AS_JM_243285_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AS_ZF_1000_bt2 <- readGAlignments("results-data/bowtie2/AS_ZF_757883_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))

EH_AC_1000_bt2 <- readGAlignments("results-data/bowtie2/EH_AC_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
EH_AS_1000_bt2 <- readGAlignments("results-data/bowtie2/EH_AS_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
EH_EH_1000_bt2 <- readGAlignments("results-data/bowtie2/EH_EH_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
EH_JM_1000_bt2 <- readGAlignments("results-data/bowtie2/EH_JM_243285_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
EH_ZF_1000_bt2 <- readGAlignments("results-data/bowtie2/EH_ZF_757883_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))

JM_AC_1000_bt2 <- readGAlignments("results-data/bowtie2/JM_AC_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
JM_AS_1000_bt2 <- readGAlignments("results-data/bowtie2/JM_AS_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
JM_EH_1000_bt2 <- readGAlignments("results-data/bowtie2/JM_EH_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
JM_JM_1000_bt2 <- readGAlignments("results-data/bowtie2/JM_JM_243285_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
JM_ZF_1000_bt2 <- readGAlignments("results-data/bowtie2/JM_ZF_757883_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))

ZF_AC_1000_bt2 <- readGAlignments("results-data/bowtie2/ZF_AC_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
ZF_AS_1000_bt2 <- readGAlignments("results-data/bowtie2/ZF_AS_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
ZF_EH_1000_bt2 <- readGAlignments("results-data/bowtie2/ZF_EH_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
ZF_JM_1000_bt2 <- readGAlignments("results-data/bowtie2/ZF_JM_243285_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
ZF_ZF_1000_bt2 <- readGAlignments("results-data/bowtie2/ZF_ZF_757883_CpG_1000bp_bt2.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))

## BOWTIE2
AC_metadata <- read.csv("results-data/sequences/AC_metadata_1000bp.csv")
AC_metadata_matched <- AC_metadata[match(names(AS_AC_1000_bt2), AC_metadata$seq),]
mcols(AS_AC_1000_bt2) <- data.frame(mcols(AS_AC_1000_bt2), AC_metadata_matched$methyl_pos,AC_metadata_matched$methyl_n)
colnames(mcols(AS_AC_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

AS_metadata <- read.csv("results-data/sequences/AS_metadata_1000bp.csv")
AS_metadata_matched <- AS_metadata[match(names(AS_AS_1000_bt2), AS_metadata$seq),]
mcols(AS_AS_1000_bt2) <- data.frame(mcols(AS_AS_1000_bt2), AS_metadata_matched$methyl_pos,AS_metadata_matched$methyl_n)
colnames(mcols(AS_AS_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

EH_metadata <- read.csv("results-data/sequences/EH_metadata_1000bp.csv")
EH_metadata_matched <- EH_metadata[match(names(AS_EH_1000_bt2), EH_metadata$seq),]
mcols(AS_EH_1000_bt2) <- data.frame(mcols(AS_EH_1000_bt2), EH_metadata_matched$methyl_pos,EH_metadata_matched$methyl_n)
colnames(mcols(AS_EH_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

JM_metadata <- read.csv("results-data/sequences/JM_metadata_243285_1000bp.csv")
JM_metadata_matched <- JM_metadata[match(names(AS_JM_1000_bt2), JM_metadata$seq),]
mcols(AS_JM_1000_bt2) <- data.frame(mcols(AS_JM_1000_bt2), JM_metadata_matched$methyl_pos,JM_metadata_matched$methyl_n)
colnames(mcols(AS_JM_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

ZF_metadata <- read.csv("results-data/sequences/ZF_metadata_7578831000bp.csv")
ZF_metadata_matched <- ZF_metadata[match(names(AS_ZF_1000_bt2), ZF_metadata$seq),]
mcols(AS_ZF_1000_bt2) <- data.frame(mcols(AS_ZF_1000_bt2), ZF_metadata_matched$methyl_pos,ZF_metadata_matched$methyl_n)
colnames(mcols(AS_ZF_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

## BOWTIE2 ##JM
AC_metadata <- read.csv("results-data/sequences/AC_metadata_1000bp.csv")
AC_metadata_matched <- AC_metadata[match(names(JM_AC_1000_bt2), AC_metadata$seq),]
mcols(JM_AC_1000_bt2) <- data.frame(mcols(JM_AC_1000_bt2), AC_metadata_matched$methyl_pos,AC_metadata_matched$methyl_n)
colnames(mcols(JM_AC_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

AS_metadata <- read.csv("results-data/sequences/AS_metadata_1000bp.csv")
AS_metadata_matched <- AS_metadata[match(names(JM_AS_1000_bt2), AS_metadata$seq),]
mcols(JM_AS_1000_bt2) <- data.frame(mcols(JM_AS_1000_bt2), AS_metadata_matched$methyl_pos,AS_metadata_matched$methyl_n)
colnames(mcols(JM_AS_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

EH_metadata <- read.csv("results-data/sequences/EH_metadata_1000bp.csv")
EH_metadata_matched <- EH_metadata[match(names(JM_EH_1000_bt2), EH_metadata$seq),]
mcols(JM_EH_1000_bt2) <- data.frame(mcols(JM_EH_1000_bt2), EH_metadata_matched$methyl_pos,EH_metadata_matched$methyl_n)
colnames(mcols(JM_EH_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

JM_metadata <- read.csv("results-data/sequences/JM_metadata_243285_1000bp.csv")
JM_metadata_matched <- JM_metadata[match(names(JM_JM_1000_bt2), JM_metadata$seq),]
mcols(JM_JM_1000_bt2) <- data.frame(mcols(JM_JM_1000_bt2), JM_metadata_matched$methyl_pos,JM_metadata_matched$methyl_n)
colnames(mcols(JM_JM_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

ZF_metadata <- read.csv("results-data/sequences/ZF_metadata_7578831000bp.csv")
ZF_metadata_matched <- ZF_metadata[match(names(JM_ZF_1000_bt2), ZF_metadata$seq),]
mcols(JM_ZF_1000_bt2) <- data.frame(mcols(JM_ZF_1000_bt2), ZF_metadata_matched$methyl_pos,ZF_metadata_matched$methyl_n)
colnames(mcols(JM_ZF_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

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


## Bowtie2

##AC
overlap_seqs_bt2_AC <- find.Overlap(AC_AC_1000_bt2,
                                    AC_AS_1000_bt2,
                                    AC_EH_1000_bt2, 
                                    AC_JM_1000_bt2,
                                    AC_ZF_1000_bt2
                                    )

##AS
overlap_seqs_bt2_AS <- find.Overlap(AS_AC_1000_bt2,
                                    AS_AS_1000_bt2,
                                    AS_EH_1000_bt2, 
                                    AS_JM_1000_bt2,
                                    AS_ZF_1000_bt2
                                    )

##EH
overlap_seqs_bt2_EH <- find.Overlap(EH_AC_1000_bt2,
                                    EH_AS_1000_bt2,
                                    EH_EH_1000_bt2, 
                                    EH_JM_1000_bt2,
                                    EH_ZF_1000_bt2
                                    )

##JM
overlap_seqs_bt2_JM <- find.Overlap(JM_AC_1000_bt2,
                                    JM_AS_1000_bt2,
                                    JM_EH_1000_bt2, 
                                    JM_JM_1000_bt2,
                                    JM_ZF_1000_bt2
                                    )

##ZF
overlap_seqs_bt2_ZF <- find.Overlap(ZF_AC_1000_bt2,
                                    ZF_AS_1000_bt2,
                                    ZF_EH_1000_bt2, 
                                    ZF_JM_1000_bt2,
                                    ZF_ZF_1000_bt2
                                    )

# unlist(lapply(overlap_seqs_bt2_AC, function(x) length(x)))
# unlist(lapply(overlap_seqs_bt2_AS, function(x) length(x)))
# unlist(lapply(overlap_seqs_bt2_EH, function(x) length(x)))
# unlist(lapply(overlap_seqs_bt2_JM, function(x) length(x)))
# unlist(lapply(overlap_seqs_bt2_ZF, function(x) length(x)))
# 
# 
# table(unlist(lapply(overlap_seqs_bt2_AC, function(x) width(x))) > 250)
# table(unlist(lapply(overlap_seqs_bt2_AS, function(x) width(x))) > 250)
# table(unlist(lapply(overlap_seqs_bt2_EH, function(x) width(x))) > 250)
# table(unlist(lapply(overlap_seqs_bt2_JM, function(x) width(x))) > 250)
# table(unlist(lapply(overlap_seqs_bt2_ZF, function(x) width(x))) > 250)
# 
# 
# table(unlist(lapply(overlap_seqs_bt2_JM, function(x) width(x))) > 100)
# table(unlist(lapply(overlap_seqs_bt2_JM, function(x) width(x))) > 250)
# table(unlist(lapply(overlap_seqs_bt2_JM, function(x) width(x))) > 250)
# table(unlist(lapply(overlap_seqs_bt2_JM, function(x) width(x))) > 250)
# overlap_seqs_bt2_filtered <- filter.GAlignments(overlap_seqs_bt2, 20, 200, input = "list")

save(overlap_seqs_bt2_AC, file = "results-data/02_conserved_seq/AC_overlaps_bt2.RData")
save(overlap_seqs_bt2_AS, file = "results-data/02_conserved_seq/AS_overlaps_bt2.RData")
save(overlap_seqs_bt2_EH, file = "results-data/02_conserved_seq/EH_overlaps_bt2.RData")
save(overlap_seqs_bt2_JM, file = "results-data/02_conserved_seq/JM_overlaps_bt2.RData")
save(overlap_seqs_bt2_ZF, file = "results-data/02_conserved_seq/ZF_overlaps_bt2.RData")
      