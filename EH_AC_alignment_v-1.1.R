#### Overview ####
# Comparing methylated sites for two species 

#### Preparation ####
# loading libraries
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(ggbio) # https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
library(dplyr)

# loading data for European hake (EH) and Atlantic cod (AC) and setting wd
setwd("/powerplant/workspace/cfngle")

# EH
EH_raw <- read.table("/powerplant/workspace/cfngle/raw-data/EH/BisRAD-CpGs-Hake.txt", sep = "\t", header = TRUE)
EH_rgenome <- readDNAStringSet("raw-data/EH/fMerMel2.1_cnag1.scaffolds.fa")

# renaming the rgenome
EH_rgenome_nvec <- names(EH_rgenome) %>%
  gsub("fMerMel2.1_cnag1_", "", .)
names(EH_rgenome) <- EH_rgenome_nvec

#AC
AC_raw <- read.table("/powerplant/workspace/cfngle/raw-data/AC/BisRAD-CpGs-Cod.txt", sep = "\t", header = TRUE)
AC_rgenome <- readDNAStringSet("raw-data/AC/GCF_902167405.1_gadMor3.0_genomic.fasta")
# renaming the rgenome
AC_rgenome_nvec <- names(AC_rgenome) %>% 
  gsub(", gadMor3.0, whole genome shotgun sequence", "", .) %>% 
  gsub(" Gadus morhua unplaced genomic scaffold", "", .)
AC_rgenome_nvec[1:23] <- gsub("^.{36}", "", AC_rgenome_nvec[1:23])
AC_rgenome_nvec <- gsub(" Gadus morhua mitochondrion, complete genome", "", AC_rgenome_nvec)
names(AC_rgenome) <- AC_rgenome_nvec

#### Manipulation Genomic Ranges ####

# EH
# load data as GRanges class
EH <- GRanges(
  seqnames = Rle(EH_raw$chr),
  ranges = IRanges(c(start = EH_raw$start), end = c(EH_raw$end), names = 1:length(EH_raw$chr)),
  strand = Rle(EH_raw$strand))

# to include the neighboring bp, "resize" can be used. to get the flanking bp, use "flank"
EH <- resize(EH,200) 

# "reduce" can be used to get merge overlapping sequences. "coverage" can be used to identify how much bp are overlapping
coverage(EH)
EH <- reduce(EH)

#AC
AC <- GRanges(
  seqnames = Rle(AC_raw$chr),
  ranges = IRanges(c(start = AC_raw$start), end = c(AC_raw$end), names = 1:length(AC_raw$chr)),
  strand = Rle(AC_raw$strand))

# to include the neighboring bp, "resize" can be used. to get the flanking bp, use "flank"
AC <- resize(AC,200) 

# "reduce" can be used to get merge overlapping sequences "coverage" can be used to identify how much bp are overlapping
coverage(AC)
AC <- reduce(AC)

#### Extraction ####

# extracting sequences using chromosome location and rgenome with the "Biostrings" package
EH_seq <- getSeq(EH_rgenome, EH)

#adding names to the sequences
names(EH_seq) <- paste0("EH_", as.character(seqnames(EH)), "_", start(EH),":", end(EH))

# saving file
writeXStringSet(EH_seq, file = "results-data/EH_CpG_200bp.fasta")

# extracting sequences using chromosome location and rgenome with the "Biostrings" package
AC_seq <- getSeq(AC_rgenome, AC)

#adding names to the sequences
names(AC_seq) <- paste0("AC_", as.character(seqnames(AC)), "_", start(AC),":", end(AC))

# saving file
writeXStringSet(AC_seq, file = "results-data/AC_CpG_200bp.fasta")


#### Visualization ####
#print(autoplot(subset(EH, seqnames(EH) == "Chr1")))
print(autoplot(EH[EH@strand == "+"&EH@seqnames == "Chr1"]))

      