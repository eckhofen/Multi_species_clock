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

# data minimap2
AC_EH_1000_mini <- readGAlignments("results-data/minimap2/AC_EH_1000_minimap.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq"))) #creates another type of object which is similar to granges object
AC_AC_1000_mini <- readGAlignments("results-data/minimap2/AC_AC_1000_minimap.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AC_AS_1000_mini <- readGAlignments("results-data/minimap2/AC_AS_1000_minimap.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AC_JM_1000_mini <- readGAlignments("results-data/minimap2/AC_JM_1000_minimap.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))

# data bowtie2 
AC_AC_1000_bt2 <- readGAlignments("results-data/bowtie2/AC_AC_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AC_AS_1000_bt2 <- readGAlignments("results-data/bowtie2/AC_AS_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AC_EH_1000_bt2 <- readGAlignments("results-data/bowtie2/AC_EH_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))
AC_JM_1000_bt2 <- readGAlignments("results-data/bowtie2/AC_JM_CpG_1000bp_bt2_.bam", use.names = TRUE, param = ScanBamParam(what = c("mapq")))


#### add metadata do GA object ####
## MINIMAP2 
AC_metadata <- read.csv("results-data/sequences/AC_metadata_1000bp.csv")
AC_metadata_matched <- AC_metadata[match(names(AC_AC_1000_mini), AC_metadata$seq),]
mcols(AC_AC_1000_mini) <- data.frame(mcols(AC_AC_1000_mini), AC_metadata_matched$methyl_pos,AC_metadata_matched$methyl_n)
colnames(mcols(AC_AC_1000_mini)) <- c("mapq", "methyl_pos", "methyl_n")

AS_metadata <- read.csv("results-data/sequences/AS_metadata_1000bp.csv")
AS_metadata_matched <- AS_metadata[match(names(AC_AS_1000_mini), AS_metadata$seq),]
mcols(AC_AS_1000_mini) <- data.frame(mcols(AC_AS_1000_mini), AS_metadata_matched$methyl_pos,AS_metadata_matched$methyl_n)
colnames(mcols(AC_AS_1000_mini)) <- c("mapq", "methyl_pos", "methyl_n")

EH_metadata <- read.csv("results-data/sequences/EH_metadata_1000bp.csv")
EH_metadata_matched <- EH_metadata[match(names(AC_EH_1000_mini), EH_metadata$seq),]
mcols(AC_EH_1000_mini) <- data.frame(mcols(AC_EH_1000_mini), EH_metadata_matched$methyl_pos,EH_metadata_matched$methyl_n)
colnames(mcols(AC_EH_1000_mini)) <- c("mapq", "methyl_pos", "methyl_n")

JM_metadata <- read.csv("results-data/sequences/JM_metadata_1000bp.csv")
JM_metadata_matched <- JM_metadata[match(names(AC_JM_1000_mini), JM_metadata$seq),]
mcols(AC_JM_1000_mini) <- data.frame(mcols(AC_JM_1000_mini), JM_metadata_matched$methyl_pos,JM_metadata_matched$methyl_n)
colnames(mcols(AC_JM_1000_mini)) <- c("mapq", "methyl_pos", "methyl_n")


## BOWTIE2
AC_metadata <- read.csv("results-data/sequences/AC_metadata_1000bp.csv")
AC_metadata_matched <- AC_metadata[match(names(AC_AC_1000_bt2), AC_metadata$seq),]
mcols(AC_AC_1000_bt2) <- data.frame(mcols(AC_AC_1000_bt2), AC_metadata_matched$methyl_pos,AC_metadata_matched$methyl_n)
colnames(mcols(AC_AC_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

AS_metadata <- read.csv("results-data/sequences/AS_metadata_1000bp.csv")
AS_metadata_matched <- AS_metadata[match(names(AC_AS_1000_bt2), AS_metadata$seq),]
mcols(AC_AS_1000_bt2) <- data.frame(mcols(AC_AS_1000_bt2), AS_metadata_matched$methyl_pos,AS_metadata_matched$methyl_n)
colnames(mcols(AC_AS_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

EH_metadata <- read.csv("results-data/sequences/EH_metadata_1000bp.csv")
EH_metadata_matched <- EH_metadata[match(names(AC_EH_1000_bt2), EH_metadata$seq),]
mcols(AC_EH_1000_bt2) <- data.frame(mcols(AC_EH_1000_bt2), EH_metadata_matched$methyl_pos,EH_metadata_matched$methyl_n)
colnames(mcols(AC_EH_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")

JM_metadata <- read.csv("results-data/sequences/JM_metadata_1000bp.csv")
JM_metadata_matched <- JM_metadata[match(names(AC_JM_1000_bt2), JM_metadata$seq),]
mcols(AC_JM_1000_bt2) <- data.frame(mcols(AC_JM_1000_bt2), JM_metadata_matched$methyl_pos,JM_metadata_matched$methyl_n)
colnames(mcols(AC_JM_1000_bt2)) <- c("mapq", "methyl_pos", "methyl_n")


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

#### Filtering ####

## using mapq score from bam files
# AC_AS_1000_bt2[mcols(AC_AS_1000_bt2)$mapq > 20]

# filtering by alignment width
# AC_AS_1000_bt2[width(AC_AS_1000_bt2) > 200 & mcols(AC_AS_1000_bt2)$mapq > 20] 

## function combining filter parameters

filter.GAlignments <- function(seq, min_mapq, min_alignwidth, input = "single") {
  if(input == "single") {
    seq <- seq[mcols(seq)$mapq >= min_mapq]
    seq <- seq[width(seq) >= min_alignwidth]
    return(seq)
  }
  else if(input == "list") {
    filtered_list <- list()
    for (i in 1:length(seq)) {
      x <- seq[[i]]
      # Apply the existing filtering criteria
      x_filter <- x[mcols(x)$mapq >= min_mapq]
      x_filter <- x_filter[width(x_filter) >= min_alignwidth]
      # Store the filtered GAlignments in the new list
      filtered_list[[i]] <- x_filter
    }
    return(filtered_list)
  }
  else {
    stop("please check input argument (\"single\" for a single sequence, \"list\" when sequences are stored in a list. Must be GAlignments object(s)")
  }
  
}


#### Run sequences through function ####

## Bowtie2
overlap_seqs_bt2 <- find.Overlap(AC_AC_1000_bt2,AC_AS_1000_bt2,AC_EH_1000_bt2, AC_JM_1000_bt2)
overlap_seqs_bt2_filtered <- filter.GAlignments(overlap_seqs_bt2, 20, 200, input = "list")

save(overlap_seqs_bt2, file = "results-data/02_conserved_seq/AC_AS_EH_JM_overlaps_bt2.R")

## Minimap2
overlap_seqs_mini <- find.Overlap(AC_AC_1000_mini,AC_AS_1000_mini,AC_EH_1000_mini,AC_JM_1000_mini)
overlap_seqs_mini_filtered <- filter.GAlignments(overlap_seqs_mini, 20, 200, input = "list")
save(overlap_seqs_mini, file = "results-data/02_conserved_seq/AC_AS_EH_JM_overlaps_mini.R")

## Combined
overlap_seqs_both <- list(overlap_seqs_bt2_filtered[[1]][names(overlap_seqs_bt2_filtered[[1]]) %in% names(overlap_seqs_mini_filtered[[1]])],
                          overlap_seqs_bt2_filtered[[2]][names(overlap_seqs_bt2_filtered[[2]]) %in% names(overlap_seqs_mini_filtered[[2]])],
                          overlap_seqs_bt2_filtered[[3]][names(overlap_seqs_bt2_filtered[[3]]) %in% names(overlap_seqs_mini_filtered[[3]])])

#### Creating csv methyldata overview ####
df_methyl_data_bt2 <- data.frame(
  c("bt2"),
  c("AC","AS","EH","JM"),
  c(sum(mcols(overlap_seqs_bt2[[1]])$methyl_n),
    sum(mcols(overlap_seqs_bt2[[2]])$methyl_n),
    sum(mcols(overlap_seqs_bt2[[3]])$methyl_n),
    sum(mcols(overlap_seqs_bt2[[4]])$methyl_n)),
  c(sum(mcols(overlap_seqs_bt2_filtered[[1]])$methyl_n),
    sum(mcols(overlap_seqs_bt2_filtered[[2]])$methyl_n),
    sum(mcols(overlap_seqs_bt2_filtered[[3]])$methyl_n),
    sum(mcols(overlap_seqs_bt2_filtered[[4]])$methyl_n)),
  c(length(overlap_seqs_bt2[[1]]),
    length(overlap_seqs_bt2[[2]]),
    length(overlap_seqs_bt2[[3]]),
    length(overlap_seqs_bt2[[4]])),
  c(length(overlap_seqs_bt2_filtered[[1]]),
    length(overlap_seqs_bt2_filtered[[2]]),
    length(overlap_seqs_bt2_filtered[[3]]),
    length(overlap_seqs_bt2_filtered[[4]])),
  c(sum(mcols(AC_AC_1000_bt2)$methyl_n),
    sum(mcols(AC_AS_1000_bt2)$methyl_n),
    sum(mcols(AC_EH_1000_bt2)$methyl_n),
    sum(mcols(AC_JM_1000_bt2)$methyl_n)),
  c(length(AC_AC_1000_bt2),
    length(AC_AS_1000_bt2),
    length(AC_EH_1000_bt2),
    length(AC_JM_1000_bt2))
)
colnames(df_methyl_data_bt2) <- c("method","species", "CpGs","CpGs_filtered","seqs_overlap","seqs_filtered", "CpGs_aligned", "seqs_aligned")

df_methyl_data_mini <- data.frame(
  c("mini"),
  c("AC","AS","EH"),
  c(sum(mcols(overlap_seqs_mini[[1]])$methyl_n),
    sum(mcols(overlap_seqs_mini[[2]])$methyl_n),
    sum(mcols(overlap_seqs_mini[[3]])$methyl_n)),
  c(sum(mcols(overlap_seqs_mini_filtered[[1]])$methyl_n),
    sum(mcols(overlap_seqs_mini_filtered[[2]])$methyl_n),
    sum(mcols(overlap_seqs_mini_filtered[[3]])$methyl_n)),
  c(length(overlap_seqs_mini[[1]]),
    length(overlap_seqs_mini[[2]]),
    length(overlap_seqs_mini[[3]])),
  c(length(overlap_seqs_mini_filtered[[1]]),
    length(overlap_seqs_mini_filtered[[2]]),
    length(overlap_seqs_mini_filtered[[3]])),
  c(sum(mcols(AC_AC_1000_mini)$methyl_n),
    sum(mcols(AC_AS_1000_mini)$methyl_n),
    sum(mcols(AC_EH_1000_mini)$methyl_n)),
  c(length(AC_AC_1000_mini),
    length(AC_AS_1000_mini),
    length(AC_EH_1000_mini))
)
colnames(df_methyl_data_mini) <- c("method","species", "CpGs","CpGs_filtered","seqs_overlap","seqs_filtered",  "CpGs_aligned", "seqs_aligned")

df_methyl_data_both <- data.frame(
  c("both"),
  c("AC","AS","EH"),
  c(NA,
    NA,
    NA),
  c(sum(mcols(overlap_seqs_both[[1]])$methyl_n),
    sum(mcols(overlap_seqs_both[[2]])$methyl_n),
    sum(mcols(overlap_seqs_both[[3]])$methyl_n)),
  c(NA,
    NA,
    NA),
  c(length(overlap_seqs_both[[1]]),
    length(overlap_seqs_both[[2]]),
    length(overlap_seqs_both[[3]])),
  c(NA,
    NA,
    NA),
  c(NA,
    NA,
    NA)
)
colnames(df_methyl_data_both) <- c("method","species", "CpGs","CpGs_filtered","seqs_overlap","seqs_filtered", "CpGs_aligned", "seqs_aligned")
# combine dfs
df_methyl_data <- rbind(df_methyl_data_bt2,df_methyl_data_mini, df_methyl_data_both)

df_save_name <- paste0(save_path, "methyl_overview.csv")
write.csv(df_methyl_data, df_save_name)

#### visualising data ####
df_long_CpG <- pivot_longer(df_methyl_data[1:3], cols = -species, names_to = "type", values_to = "count")
df_long_seq <- pivot_longer(df_methyl_data[c(-2:-3)], cols = -species, names_to = "type", values_to = "count")

ggplot(df_long_CpG, aes(x = species, y = count, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + 
  theme_classic() +
  labs(title = "Comparison of Original and Filtered Counts",
       x = "Species",
       y = "Count",
       fill = "Measurement")

ggplot(df_long_seq, aes(x = species, y = count, fill = type)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Set2") + 
  theme_classic() +
  labs(title = "Comparison of Original and Filtered Counts",
       x = "Species",
       y = "Count",
       fill = "Measurement")

#### other things ####

# test_overlaps <- findOverlaps(ZF_EH_1000_mini, ZF_AS_1000_mini) # returns the indexes of overlapping sequences 
# test_overlaps_count <- countOverlaps(ZF_EH_1000_mini, ZF_AS_1000_mini) # count how many ranges were overlapping how many times
# filtered_overlaps <- test_overlaps[test_overlaps_count>0]

ZF_EH_1000_AC_ol_index <- findOverlaps(ZF_EH_1000_mini, ZF_AC_1000_mini)
ZF_EH_1000_AS_ol_index <- findOverlaps(ZF_EH_1000_mini, ZF_AS_1000_mini)

ZF_AC_1000_AS_ol_index <- findOverlaps(ZF_AC_1000_mini, ZF_AS_1000_mini)
ZF_AC_1000_EH_ol_index <- findOverlaps(ZF_AC_1000_mini, ZF_EH_1000_mini)

ZF_AS_1000_AC_ol_index <- findOverlaps(ZF_AS_1000_mini, ZF_AC_1000_mini)
ZF_AS_1000_EH_ol_index <- findOverlaps(ZF_AS_1000_mini, ZF_EH_1000_mini)

# checking which sequences are overlapping in EH_AC and EH_AS
unique_ZF_EH_AC_AS_index_EH <- unique(queryHits(ZF_EH_1000_AC_ol_index)[queryHits(ZF_EH_1000_AC_ol_index) %in% queryHits(ZF_EH_1000_AS_ol_index)])
unique_ZF_EH_AC_AS_index_AC <- unique(queryHits(ZF_AC_1000_AS_ol_index)[queryHits(ZF_AC_1000_AS_ol_index) %in% queryHits(ZF_AC_1000_EH_ol_index)])
unique_ZF_EH_AC_AS_index_AS <- unique(queryHits(ZF_AS_1000_AC_ol_index)[queryHits(ZF_AS_1000_AC_ol_index) %in% queryHits(ZF_AS_1000_EH_ol_index)])

# making sure that the queryHits have the same number and objects in both AC and AS 
unique(queryHits(findOverlaps(ZF_EH_1000_mini[unique_ZF_EH_AC_AS_index_EH], ZF_AC_1000_mini)))
unique(queryHits(findOverlaps(ZF_EH_1000_mini[unique_ZF_EH_AC_AS_index_EH], ZF_AS_1000_mini)))

# subjectHits do have different values and amounts because the sequences are overlapping various times
unique(subjectHits(findOverlaps(ZF_EH_1000_mini[unique_ZF_EH_AC_AS_index_EH], ZF_AC_1000_mini)))
unique(subjectHits(findOverlaps(ZF_EH_1000_mini[unique_ZF_EH_AC_AS_index_EH], ZF_AS_1000_mini)))

# getting sequence names which were overlapping
seqnames_AC_1000 <-names(ZF_AC_1000_mini[unique_ZF_EH_AC_AS_index_AC])
seqnames_AS_1000 <-names(ZF_AS_1000_mini[unique_ZF_EH_AC_AS_index_AS])
seqnames_EH_1000 <- names(ZF_EH_1000_mini[unique_ZF_EH_AC_AS_index_EH])

#### getting sequences which overlap ####

#loading sequences
AC_1000_seq <- readDNAStringSet("results-data/sequences/AC_CpG_1000bp.fasta")
AS_1000_seq <- readDNAStringSet("results-data/sequences/AS_CpG_1000bp.fasta")
EH_1000_seq <- readDNAStringSet("results-data/sequences/EH_CpG_1000bp.fasta")

conserved_AC_1000 <- AC_1000_seq[AC_1000_seq@ranges@NAMES %in% unique(seqnames_AC_1000)]
conserved_AS_1000 <- AS_1000_seq[AS_1000_seq@ranges@NAMES %in% unique(seqnames_AS_1000)]
conserved_EH_1000 <- EH_1000_seq[EH_1000_seq@ranges@NAMES %in% unique(seqnames_EH_1000)]

# write sequences as fasta
AC_1000_filename <- paste0(save_path, "AC_1000_conserved", suffix)
writeXStringSet(conserved_AC_1000, AC_1000_filename)

AS_1000_filename <- paste0(save_path, "AS_1000_conserved", suffix)
writeXStringSet(conserved_AS_1000, AS_1000_filename)

EH_1000_filename <- paste0(save_path, "EH_1000_conserved", suffix)
writeXStringSet(conserved_EH_1000, EH_1000_filename)

#combine fasta files

AC_AS_EH_1000_filename <- paste0(save_path, "AC_AS_EH_1000_conserved", suffix)
writeXStringSet(c(conserved_AC_1000,conserved_AS_1000,conserved_EH_1000), AC_AS_EH_1000_filename)


# xx <- as.data.frame(mergeByOverlaps(ZF_EH_1000_mini, ZF_AC_1000_mini)) #merges overlapping sequences (non overlapping are retained as well) 

#### overlap settings ####
# this is identifying what the amount of bps are which are overlapping
# Extract the overlapping ranges
overlappingRanges <- pintersect(granges(ZF_AC_1000_mini)[queryHits(ZF_AC_1000_AS_ol_index)], granges(ZF_AS_1000_mini)[subjectHits(ZF_AC_1000_AS_ol_index)])

# Calculate widths of the overlaps
overlapWidths <- width(overlappingRanges)
length(overlapWidths[overlapWidths > 600])

#### visualizing overlaps ####
df <- data.frame(overlapWidths[overlapWidths > 600])
ggplot(df, aes(x = df$overlapWidths.overlapWidths...600.)) +
  geom_histogram(binwidth = 10, fill = "darkblue") + # Adjust binwidth as needed
  theme_classic() +
  labs(title = paste0("Distribution of Overlap Lengths n= ", length(df$overlapWidths)),
       x = "Overlap Length",
       y = "Frequency")


overlap_ZF_EH_AC <- ZF_EH_1000_mini[unique(queryHits(ZF_EH_1000_AC_ol_index))]
overlap_ZF_EH_AC_AS_index <- findOverlaps(overlap_ZF_EH_AC, ZF_AS_1000_mini)
overlap_ZF_EH_AC_AS <- overlap_ZF_EH_AC[unique(queryHits(overlap_ZF_EH_AC_AS_index))]


# ZF_EH_1000_mini[queryHits(ZF_EH_1000_AS_ol_index)]


# queryHits() # This returns the indexes of the query (the first sequence) which was overlapping with something
# subjectHits()

ZF_EH_1000_mini[queryHits(test_overlaps)] 
ZF_AC_1000_mini[subjectHits(test_overlaps)]


# ZF_EH_5000 <- scanBam("results-data/minimap2/ZF_EH_5000_minimap.bam")



#### function ####
#### seq.len must be changed! The width refers to the aligned length and not sequence length

extract.seq <- function(seq.read, seq.len, name = "test", save.path = save_path, suffix = ".fasta", write = FALSE) {
  extracted_seq <- seq.read[[1]][["seq"]][!is.na(seq.read[[1]][["pos"]]) & width(seq.read[[1]][["seq"]]) >= seq.len]
  names(extracted_seq) <- seq.read[[1]][["qname"]][!is.na(seq.read[[1]][["pos"]]) & width(seq.read[[1]][["seq"]]) >= seq.len]
  if(write == TRUE){
    file_name <- paste0(save.path, name,"_conserved", suffix)
    writeXStringSet(extracted_seq, file = file_name)
  }
  return(extracted_seq)
}

#### Extraction ####

extract.seq(ZF_EH_100_mini, 100, name = "ZF_EH_100_mini", write = TRUE)
extract.seq(ZF_EH_100_bt2, 100, name = "ZF_EH_100_bt2", write = TRUE)

ZF_EH_1000_mini <- extract.seq(ZF_EH_1000_mini, 1000)
table(width(ZF_EH_1000_mini)>999)
sort(ZF_EH_1000_mini[[1]][["pos"]][!is.na(ZF_EH_1000_mini[[1]][["pos"]]) & ZF_EH_1000_mini[[1]][["rname"]]=="NC_007112.7"])

# extracted_seq <- ZF_EH_100_mini[[1]][["seq"]][(ZF_EH_100_mini[[1]][["flag"]]) != 4]
# 
# extracted_seq <- ZF_EH_5000[[1]][["seq"]][!is.na(ZF_EH_5000[[1]][["pos"]]) & ZF_EH_5000[[1]][["qwidth"]] > 0]
# 
# names(extracted_seq) <- ZF_EH_5000[[1]][["qname"]][!is.na(ZF_EH_5000[[1]][["pos"]])]
# writeXStringSet(extracted_seq, file = paste0(save.path, name,"_conserved", suffix))




# extract.seq <- function(seq.read, name, save.path = save_path,  suffix = suffix) {
#   extracted_seq <- seq.read[[1]][["seq"]][!is.na(seq.read[[1]][["pos"]])]
#   names(extracted_seq) <- seq.read[[1]][["qname"]][!is.na(seq.read[[1]][["pos"]])]
#   writeXStringSet(extracted_seq, file = paste0(save.path, name,"_conserved", suffix))
# }
# 
# extract.seq(ZF_EH_5000, ZF_EH_5000, save_path, suffix)

writeXStringSet(extracted_seq, file = "/workspace/cfngle/results-data/02_conserved_seq/ZF_EH_100_conserved.fasta")
                      