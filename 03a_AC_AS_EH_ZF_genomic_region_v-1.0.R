library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)

setwd("/powerplant/workspace/cfngle")

# defining objects 
save_path <- "/workspace/cfngle/results-data/02_conserved_seq/"
suffix <- ".fasta"
AC_annotations <- import("raw-data/AC/rgenome/annotations/ncbi_dataset/data/GCF_902167405.1/genomic.gff")

AC_all_mini <- c(overlap_seqs_mini[[1]],overlap_seqs_mini[[3]],overlap_seqs_mini[[3]])
sort(table(seqnames(AC_all_mini)))

AC_all_mini_ano_overlaps_index <- findOverlaps(AC_all_mini, AC_annotations)
AC_all_mini_ano_overlaps <- AC_annotations[subjectHits(AC_all_mini_ano_overlaps_index)]

table(mcols(AC_all_mini_ano_overlaps)$type)

AC_all_mini_genes <- AC_all_mini_ano_overlaps[mcols(AC_all_mini_ano_overlaps)$type == "gene"]
AC_all_mini_genes_names <- unique(mcols(AC_all_mini_genes)$Name)

##bowtie2
AC_all_bt2 <- c(overlap_seqs_bt2[[1]],overlap_seqs_bt2[[3]],overlap_seqs_bt2[[3]])
sort(table(seqnames(AC_all_bt2)))

AC_all_bt2_ano_overlaps_index <- findOverlaps(AC_all_bt2, AC_annotations)
AC_all_bt2_ano_overlaps <- AC_annotations[subjectHits(AC_all_bt2_ano_overlaps_index)]

table(mcols(AC_all_bt2_ano_overlaps)$type)

AC_all_bt2_genes <- AC_all_bt2_ano_overlaps[mcols(AC_all_bt2_ano_overlaps)$type == "gene"]
AC_all_bt2_genes_names <- unique(mcols(AC_all_bt2_genes)$Name)


unique(mcols(AC_annotations)$type)


test_granges <- granges(overlap_seqs_mini_filtered[[2]])
