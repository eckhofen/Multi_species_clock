#### Overview ####
# Extraction of methylation data based on overlapping sequences

#### Settings ####
setwd("/powerplant/workspace/cfngle")

#### Preparation ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(tidyverse)

## loading data 

load("results-data/02_conserved_seq/AC_AS_EH_overlaps_bt2.R")
load("results-data/02_conserved_seq/AC_AS_EH_overlaps_mini.R")

overlap_seqs_bt2
overlap_seqs_mini

overlap_AC_AS_bt2 <- overlap_seqs_bt2[[2]]
overlap_AC_EH_bt2 <- overlap_seqs_bt2[[3]]



# getting start of alignment of CIGAR
cigar_sep <- cigar(overlap_AC_AS_bt2) %>% 
  sapply(function(x) regmatches(x, gregexpr("\\d+[A-Z]", x)))

start_AS_AC <- lapply(cigar_sep, function(x) if(grepl("S$", x[1])) {
  x[1] <- as.integer(strsplit(x[1], split = "S")[1])
}) %>% 
  gsub('NULL', '0',.) %>% as.integer()

# getting end of alignment of CIGAR
end_AS_AC <- lapply(cigar_sep, function(x)
  if(grepl("S$", x[length(x)])) {
    x[length(x)] <- as.integer(strsplit(x[length(x)], split = "S")[1])
  }) %>% 
  gsub('NULL', '0',.) %>% as.integer()

seq_end <- qwidth(overlap_AC_AS_bt2) - end_AS_AC 

df_start_end <- data.frame(row.names = names(overlap_AC_AS_bt2),start_AS_AC, seq_end, seq_end-start_AS_AC, width(overlap_AC_AS_bt2),width(overlap_AC_AS_bt2)-(seq_end-start_AS_AC))
names(df_start_end) <- c("start","end","width","width_align","diff")



## getting methylation position on aligned sequences 
# getting the start of sequence which was aligned
seq_start_pos <- names(overlap_AC_AS_bt2) %>% 
  gsub(":.*$", "",.) %>% 
  gsub("^.+_", "",.) %>% as.integer()

seq_chr_name <- names(overlap_AC_AS_bt2) %>% 
  gsub("AS_", "",.) %>% 
  gsub("_[^_]*$", "",., perl = TRUE)

# getting methylation site an normalizing to the aligned sequence part
methyl_sites <- mcols(overlap_AC_AS_bt2)$methyl_pos %>% 
  strsplit(",") %>% lapply(function(x) as.integer(x))
normalized_methyl_sites <- mapply(function(x, start_pos) x - start_pos, methyl_sites, seq_start_pos, SIMPLIFY = FALSE)
mapped_methyl_sites <- mapply(function(x, start_pos) x - start_pos, normalized_methyl_sites, df_start_end$start, SIMPLIFY = FALSE)

methyl_on_aligned <- mapply(function(x, w) x <= w & x > 0, mapped_methyl_sites, width(overlap_AC_AS_bt2), SIMPLIFY = FALSE)

# giving the elements in the list the corresponding chromosome names 
Map(function(l, v) list(l, list(v)), mapped_methyl_sites, seq_chr_name)

final_methyl_sites <- lapply(seq_along(mapped_methyl_sites), function(i) {
  df <- data.frame(pos_rgenome = methyl_sites[[i]],pos_aglin = mapped_methyl_sites[[i]], width_align = width(overlap_AC_AS_bt2)[i], chr_align = seqnames(overlap_AC_AS_bt2)[i], aligned = methyl_on_aligned[[i]])
  df$Chr <- seq_chr_name[i]
  df <- cbind(df[length(df)], df[-length(df)])
  return(df)
})


testt <- cigar(overlap_AC_AS_bt2[names(overlap_AC_AS_bt2)=="AS_LG3_4867166:4868268"])
regmatches(testt, gregexpr("\\d+[A-Z]", testt))
           
filtered_mapped_methyl_sites <- mapped_methyl_sites


