#### Overview ####
# Extraction of sequences of methylation data for the AC, EH and AS

#### Settings ####
bp_ext <- 200 # this will be the length of the extracted sequences around the CpG 
save_folder <- "/powerplant/workspace/cfngle/results-data/" # folder where extracted sequences will be saved
file_ext <- ".fasta" # which file extension will be used for the sequences

#### Preparation ####
# loading libraries
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(ggbio) # https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
library(dplyr)

# loading data for European hake (EH) and Atlantic cod (AC) and setting wd
setwd("/powerplant/workspace/cfngle")

## EH
EH_raw <- read.table("/powerplant/workspace/cfngle/raw-data/EH/BisRAD-CpGs-Hake.txt", sep = "\t", header = TRUE)
EH_rgenome <- readDNAStringSet("raw-data/EH/fMerMel2.1_cnag1.scaffolds.fa")

# renaming the rgenome
EH_rgenome_nvec <- names(EH_rgenome) %>%
  gsub("fMerMel2.1_cnag1_", "", .)
names(EH_rgenome) <- EH_rgenome_nvec

## AC
AC_raw <- read.table("/powerplant/workspace/cfngle/raw-data/AC/BisRAD-CpGs-Cod.txt", sep = "\t", header = TRUE)
AC_rgenome <- readDNAStringSet("raw-data/AC/GCF_902167405.1_gadMor3.0_genomic.fasta")
# renaming the rgenome
AC_rgenome_nvec <- names(AC_rgenome) %>% 
  gsub(", gadMor3.0, whole genome shotgun sequence", "", .) %>% 
  gsub(" Gadus morhua unplaced genomic scaffold", "", .)
AC_rgenome_nvec[1:23] <- gsub("^.{36}", "", AC_rgenome_nvec[1:23])
AC_rgenome_nvec <- gsub(" Gadus morhua mitochondrion, complete genome", "", AC_rgenome_nvec)
names(AC_rgenome) <- AC_rgenome_nvec

## AS
AS_raw <- read.table("/powerplant/workspace/cfngle/raw-data/AS/BisRAD-CpGs-Snapper.txt", sep = "\t", header = TRUE)
AS_rgenome <- readDNAStringSet("raw-data/AS/Chrysophrys_auratus.v.1.0.all.male.map.fasta")
# renaming the rgenome 
AS_rgenome_nvec <- names(AS_rgenome) %>% 
  gsub(" size.*$", "", .)
names(AS_rgenome) <- AS_rgenome_nvec

#### Manipulation Genomic Ranges ####

## EH
# load data as GRanges class
EH <- GRanges(
  seqnames = Rle(EH_raw$chr),
  ranges = IRanges(c(start = EH_raw$start), end = c(EH_raw$end), names = 1:length(EH_raw$chr)),
  strand = Rle(EH_raw$strand))

# to include the neighboring bp, "resize" can be used. to get the flanking bp, use "flank"
EH <- resize(EH,bp_ext) 

# "reduce" can be used to get merge overlapping sequences. "coverage" can be used to identify how much bp are overlapping
coverage(EH)
EH <- reduce(EH)

## AC
AC <- GRanges(
  seqnames = Rle(AC_raw$chr),
  ranges = IRanges(c(start = AC_raw$start), end = c(AC_raw$end), names = 1:length(AC_raw$chr)),
  strand = Rle(AC_raw$strand))

# to include the neighboring bp, "resize" can be used. to get the flanking bp, use "flank"
AC <- resize(AC,bp_ext) 

# "reduce" can be used to get merge overlapping sequences "coverage" can be used to identify how much bp are overlapping
coverage(AC)
AC <- reduce(AC)

## AS
AS <- GRanges(
  seqnames = Rle(AS_raw$chr),
  ranges = IRanges(c(start = AS_raw$start), end = c(AS_raw$end), names = 1:length(AS_raw$chr)),
  strand = Rle(AS_raw$strand))

# to include the neighboring bp, "resize" can be used. to get the flanking bp, use "flank"
AS <- resize(AS,bp_ext) 

# in this case this led to some sequences being positioned in nonexistent indexes (negative).
# This happens when they are close to the start and then being extended 

# function to fix this problem (see explanation at appendix)
fix.seq() <- function(seq, rgenome, seq_width) {
  seq[seq@ranges@start <= 0] <- shift(seq[seq@ranges@start <= 0], -1*c(seq@ranges@start[seq@ranges@start <= 0])+1)
  matching_rgenome <- rgenome[unique(rgenome@ranges@NAMES) %in% unique(seq@seqnames@values)]
  rg_max_length <- width(matching_rgenome[seq.character(seqnames(seq))])
  end(seq) <- pmin(end(seq), rg_max_length)
  start(seq[width(seq) < seq_width]) <- start(seq[width(seq) < seq_width]) - c(seq_width-width(seq[width(seq) < seq_width]))
}

# "reduce" can be used to get merge overlapping sequences "coverage" can be used to identify how much bp are overlapping
coverage(AS)
AS <- reduce(AS)

#### Check ####
# this just makes sure all the sequenced locations are also represented in the rgenome
unique(AS_raw$chr) %in% names(AS_rgenome)
unique(AC_raw$chr) %in% names(AC_rgenome)
unique(EH_raw$chr) %in% names(EH_rgenome)

#### Extraction ####

## EH
# extracting sequences using chromosome location and rgenome with the "Biostrings" package
EH_seq <- getSeq(EH_rgenome, EH)

#adding names to the sequences
names(EH_seq) <- paste0("EH_", as.character(seqnames(EH)), "_", start(EH),":", end(EH))

# saving file
EH_filename <- paste0(save_folder, "EH_CpG_", bp_ext, "bp", file_ext)
writeXStringSet(EH_seq, file = EH_filename)

## AC
# extracting sequences using chromosome location and rgenome with the "Biostrings" package
AC_seq <- getSeq(AC_rgenome, AC)

#adding names to the sequences
names(AC_seq) <- paste0("AC_", as.character(seqnames(AC)), "_", start(AC),":", end(AC))

# saving file
AC_filename <- paste0(save_folder, "AC_CpG_", bp_ext, "bp", file_ext)
writeXStringSet(AC_seq, file = AC_filename)

## AS
# extracting sequences using chromosome location and rgenome with the "Biostrings" package
AS_seq <- getSeq(AS_rgenome, AS)

#adding names to the sequences
names(AS_seq) <- paste0("AS_", as.character(seqnames(AS)), "_", start(AS),":", end(AS))

# saving file
AS_filename <- paste0(save_folder, "AS_CpG_", bp_ext, "bp", file_ext)
writeXStringSet(AS_seq, file = AS_filename)


#### Visualization and testing ####
#print(autoplot(subset(EH, seqnames(EH) == "Chr1")))
print(autoplot(EH[EH@strand == "+"&EH@seqnames == "Chr1"]))

EH_rev <- reverseComplement(EH_seq)
writeXStringSet(EH_seq, file = "/workspace/cfngle/results-data/EH_CpG_100bp_rev.fasta")
      
#### Appendix ####

## function to fix out-of-boundary sequences
"""
# FIX start: this shifts the sequences which are out of range into the sequence range again so that it starts with 1
AS[AS@ranges@start <= 0] <- shift(AS[AS@ranges@start <= 0], -1*c(AS@ranges@start[AS@ranges@start <= 0])+1)

# FIX end: the same goes for the ends of the sequences. Some have reached beyond the reference genome
# filtering only the regions of the rgenome which are also in the methylation data
AS_matching <- AS_rgenome[unique(AS_rgenome@ranges@NAMES) %in% unique(AS@seqnames@values)]

# extracting the maximum possible sequence lengths for the methylation data 
AS_max_length <- width(AS_matching[as.character(seqnames(AS))])

# updating the ends of AS
end(AS) <- pmin(end(AS), AS_max_length)

# adjusting length to bp_ext again 
start(AS[width(AS) < bp_ext]) <- start(AS[width(AS) < bp_ext]) - c(bp_ext-width(AS[width(AS) < bp_ext]))
"""