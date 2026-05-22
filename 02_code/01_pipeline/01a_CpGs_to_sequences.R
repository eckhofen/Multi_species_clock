# Metadata ----------------------------------------------------------------
# Project: Shared Methylation Regions
# Description: Extracting 1000bp windows from CpG sites 
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026-05

#### Settings ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(ggbio) # https://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf
library(dplyr)
library(ggplot2)
library(tidyr)

# Genomic window
bp_ext <- 1000 
data_folder <- "01_data/01_raw_private/"
save_folder <- "01_data/02_external_server_outputs/01_sequences" 

#### Loading data ####
## EH
EH_raw <- read.table("01_data/01_raw_private/BisRAD-CpGs-Hake.txt", sep = "\t", header = TRUE)
EH_rgenome <- readDNAStringSet("01_data/01_raw_private/genomes/fMerMel2.1_cnag1.scaffolds.fa")
EH_rgenome_nvec <- names(EH_rgenome) %>%
  gsub("fMerMel2.1_cnag1_", "", .)
names(EH_rgenome) <- EH_rgenome_nvec

## AC
AC_raw <- read.table("01_data/01_raw_private//BisRAD-CpGs-Cod.txt", sep = "\t", header = TRUE)
AC_rgenome <- readDNAStringSet("01_data/01_raw_private/genomes/GCF_902167405.1_gadMor3.0_genomic.fasta")
AC_rgenome_nvec <- names(AC_rgenome) %>% 
  gsub(", gadMor3.0, whole genome shotgun sequence", "", .) %>% 
  gsub(" Gadus morhua unplaced genomic scaffold", "", .)
AC_rgenome_nvec[1:23] <- gsub("^.{36}", "", AC_rgenome_nvec[1:23])
AC_rgenome_nvec <- gsub(" Gadus morhua mitochondrion, complete genome", "", AC_rgenome_nvec)
names(AC_rgenome) <- AC_rgenome_nvec

## AS
AS_raw <- read.table("01_data/01_raw_private/BisRAD-CpGs-Snapper.txt", sep = "\t", header = TRUE)
AS_rgenome <- readDNAStringSet("01_data/01_raw_private/genomes/Chrysophrys_auratus.v.1.0.all.male.map.fasta")
AS_rgenome_nvec <- names(AS_rgenome) %>% 
  gsub(" size.*$", "", .)
names(AS_rgenome) <- AS_rgenome_nvec

## ZF
load("01_data/01_raw_private/ZF_methylpos.RData")
ZF_raw <- ZF_methyl_pos
ZF_rgenome <- readDNAStringSet("01_data/01_raw_private/genomes/GCF_000002035.6_GRCz11_genomic.fasta")
ZF_rgenome_nvec <- names(ZF_rgenome) %>% 
  gsub(" Danio.*$", "", .)
names(ZF_rgenome) <- ZF_rgenome_nvec

#### Extending genomic windows ####
# helper functions
df.to.granges <- function(df) {
  granges_obj <- GRanges(
    seqnames = Rle(df$chr),
    ranges = IRanges(
      c(start = df$start), 
      end = c(df$end), 
      names = 1:length(df$chr)),
    strand = Rle(df$strand))
  return(granges_obj)
}

# function to fix overextending sequences
fix.seq <- function(seq, rgenome, seq_width) {
  seq[seq@ranges@start <= 0] <- shift(seq[seq@ranges@start <= 0], -1*c(seq@ranges@start[seq@ranges@start <= 0])+1)
  matching_rgenome <- rgenome[unique(rgenome@ranges@NAMES) %in% unique(seq@seqnames@values)]
  rg_max_length <- width(matching_rgenome[as.character(seqnames(seq))])
  end(seq) <- pmin(end(seq), rg_max_length)
  return(seq)
}

extend.windows <- function(granges_obj, rgenome, bp_ext) {
  granges_obj <- fix.seq(granges_obj, rgenome, bp_ext)
  granges_obj <- GenomicRanges::reduce(granges_obj)
  return(granges_obj)
}

# Run for all species
# AC
AC_methyl <- df.to.granges(AC_raw)
AC <- extend.windows(AC_methyl, AC_rgenome, bp_ext)

# AS
AS_methyl <- df.to.granges(AS_raw)
AS <- extend.windows(AS_methyl, AS_rgenome, bp_ext)

# EH
EH_methyl <- df.to.granges(EH_raw)
EH <- extend.windows(EH_methyl, EH_rgenome, bp_ext)

# ZF
ZF_methyl <- df.to.granges(ZF_raw)
ZF <- extend.windows(ZF_methyl, ZF_rgenome, bp_ext)

# Sanity check
unique(AS_raw$chr) %in% names(AS_rgenome)
unique(AC_raw$chr) %in% names(AC_rgenome)
unique(EH_raw$chr) %in% names(EH_rgenome)
unique(ZF_raw$chr) %in% names(ZF_rgenome)

#### Extraction ####

# helper functions
extract.sequences <- function(granges_obj, rgenome, species) {
  seq <- getSeq(rgenome, granges_obj)
  names(seq) <- paste0(species, "_", as.character(seqnames(granges_obj)), "_", start(granges_obj),":", end(granges_obj))
  return(seq)
}

# creating dataframe for methylation positions
create.MethylPos <- function(seqs, seqs_GR, methylsites, name = "CpGs") {
  overlaps <- findOverlaps(seqs_GR, methylsites)
  seq_names <- names(seqs)
  methyl_sites <- start(methylsites[subjectHits(overlaps)])
  concat_methyl_sites <- tapply(methyl_sites, INDEX = queryHits(overlaps), FUN = function(x) paste0(x, collapse = ","))
  df <- data.frame(
    seq = seq_names,
    methyl_pos = concat_methyl_sites,
    methyl_n = as.vector(table(queryHits(overlaps)))
  )
  return(df)
}

# EH
EH_seq <- extract.sequences(EH, EH_rgenome, "EH")
EH_filename <- paste0(save_folder, "EH_CpG_", bp_ext, "bp.fasta")
writeXStringSet(EH_seq, file = EH_filename)

# metadata
EH_metadata <- create.MethylPos(EH_seq, EH, EH_methyl)
EH_metadata_filename <- paste0(save_folder, "EH_metadata_", bp_ext, "bp.csv")
write.csv(EH_metadata, EH_metadata_filename)

# AC
AC_seq <- extract.sequences(AC, AC_rgenome, "AC")
AC_filename <- paste0(save_folder, "AC_CpG_", bp_ext, "bp.fasta")
writeXStringSet(AC_seq, file = AC_filename)

# metadata
AC_metadata <- create.MethylPos(AC_seq, AC, AC_methyl)
AC_metadata_filename <- paste0(save_folder, "AC_metadata_", bp_ext, "bp.csv")
write.csv(AC_metadata, AC_metadata_filename)

# AS
AS_seq <- extract.sequences(AS, AS_rgenome, "AS")
AS_filename <- paste0(save_folder, "AS_CpG_", bp_ext, "bp.fasta")
writeXStringSet(AS_seq, file = AS_filename)

# metadata
AS_metadata <- create.MethylPos(AS_seq, AS, AS_methyl)
AS_metadata_filename <- paste0(save_folder, "AS_metadata_", bp_ext, "bp.csv")
write.csv(AS_metadata, AS_metadata_filename)

# ZF
ZF_seq <- getSeq(ZF_rgenome, ZF)
names(ZF_seq) <- paste0("ZF_", as.character(seqnames(ZF)), "_", start(ZF),":", end(ZF))
ZF_filename <- paste0(save_folder, "ZF_88_CpG_", bp_ext, "bp.fasta")
writeXStringSet(ZF_seq, file = ZF_filename)

# metadata
ZF_metadata <- create.MethylPos(ZF_seq, ZF, ZF_methyl)
ZF_metadata_filename <- paste0(save_folder, "ZF_metadata_", bp_ext, "bp.csv")
write.csv(ZF_metadata, ZF_metadata_filename)