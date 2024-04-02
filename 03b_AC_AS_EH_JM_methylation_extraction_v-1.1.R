#### Overview ####
# Extraction of methylation data based on overlapping sequences

#### Settings ####
setwd("/powerplant/workspace/cfngle")

#### Preparation ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(tidyverse)

## loading data 

load("results-data/02_conserved_seq/AC_AS_EH_JM_overlaps_bt2.R")
load("results-data/02_conserved_seq/AC_AS_EH_JM_overlaps_mini.R")

overlap_seqs_bt2
overlap_seqs_mini

overlap_AC_AC_bt2 <- overlap_seqs_bt2[[1]]
overlap_AC_AS_bt2 <- overlap_seqs_bt2[[2]]
overlap_AC_EH_bt2 <- overlap_seqs_bt2[[3]]
overlap_AC_JM_bt2 <- overlap_seqs_bt2[[4]]

overlap_AC_AC_mini <- overlap_seqs_mini[[1]]
overlap_AC_AS_mini <- overlap_seqs_mini[[2]]
overlap_AC_EH_mini <- overlap_seqs_mini[[3]]
overlap_AC_JM_mini <- overlap_seqs_mini[[4]]

#### Extracting alignment data from CIGAR code ####
# (see appendix)
## function
cigar.to.width <- function(seqs_aligned) {
  cigar_sep <- cigar(seqs_aligned) %>% 
    sapply(function(x) regmatches(x, gregexpr("\\d+[A-Z]", x)))
  
  start_temp <- lapply(cigar_sep, function(x) if(grepl("S$", x[1])) {
    x[1] <- as.integer(strsplit(x[1], split = "S")[1])
  }) %>% 
    gsub('NULL', '0',.) %>% as.integer()
  
  # getting end of alignment of CIGAR
  end_AS_AC <- lapply(cigar_sep, function(x)
    if(grepl("S$", x[length(x)])) {
      x[length(x)] <- as.integer(strsplit(x[length(x)], split = "S")[1])
    }) %>% 
    gsub('NULL', '0',.) %>% as.integer()
  
  seq_end <- qwidth(seqs_aligned) - end_AS_AC 
  
  cigar_width_df <- data.frame(seq_names = names(seqs_aligned),start_temp, seq_end, seq_end-start_temp, width(seqs_aligned),width(seqs_aligned)-(seq_end-start_temp))
  names(cigar_width_df) <- c("seq_names", "start","end","width","width_align","diff")
  return(cigar_width_df)
}

#### getting methylation position on aligned sequences ####
# see appendix for explanation
### function
map.methyl.to.align <- function(seqs_aligned, cigar_width_df, species = "undefined") {
  # getting the start of sequence which was aligned
  seq_start_pos <- names(seqs_aligned) %>% 
    gsub(":.*$", "",.) %>% 
    gsub("^.+_", "",.) %>% as.integer()
  
  seq_chr_name <- names(seqs_aligned) %>% 
    gsub("AS_", "",.) %>% 
    gsub("_[^_]*$", "",., perl = TRUE)
  
  # getting methylation site an normalizing to the aligned sequence part
  methyl_sites <- mcols(seqs_aligned)$methyl_pos %>% 
    strsplit(",") %>% lapply(function(x) as.integer(x))
  normalized_methyl_sites <- mapply(function(x, start_pos) x - start_pos, methyl_sites, seq_start_pos, SIMPLIFY = FALSE)
  mapped_methyl_sites <- mapply(function(x, start_pos) x - start_pos, normalized_methyl_sites, cigar_width_df$start, SIMPLIFY = FALSE)
  aligned_methyl_sites <- mapply(function(x, start_pos, align_pos) x - start_pos + align_pos, mapped_methyl_sites, cigar_width_df$start, start(seqs_aligned), SIMPLIFY = FALSE)
  
  methyl_on_aligned <- mapply(function(x, w) x <= w & x > 0, mapped_methyl_sites, width(seqs_aligned), SIMPLIFY = FALSE)
  
  final_methyl_sites_AS <- lapply(seq_along(mapped_methyl_sites), function(i) {
    df <- data.frame(seq_names = cigar_width_df$seq_names[i], 
                     pos_rgenome = methyl_sites[[i]],
                     pos_seq = mapped_methyl_sites[[i]],
                     pos_align = aligned_methyl_sites[[i]], 
                     width_align = width(seqs_aligned)[i], 
                     chr_align = seqnames(seqs_aligned)[i], 
                     aligned = methyl_on_aligned[[i]], 
                     seq_number = sprintf("seq_%03d", 1:length(aligned_methyl_sites))[i], 
                     species = species)
    df$Chr <- seq_chr_name[i]
    df <- cbind(df[length(df)], df[-length(df)])
    return(df)
  })
  
  df_final_methyl_sites_AS <- bind_rows(final_methyl_sites_AS)
  return(df_final_methyl_sites_AS)
  
}

### master function

get.methyl.sites <- function(seqs_aligned, species = "undefined") {
  ### A) extractinformation from CIGAR code
  cigar_sep <- cigar(seqs_aligned) %>% 
    sapply(function(x) regmatches(x, gregexpr("\\d+[A-Z]", x)))
  
  start_temp <- lapply(cigar_sep, function(x) if(grepl("S$", x[1])) {
    x[1] <- as.integer(strsplit(x[1], split = "S")[1])
  }) %>% 
    gsub('NULL', '0',.) %>% as.integer()
  
  # getting end of alignment of CIGAR
  end_AS_AC <- lapply(cigar_sep, function(x)
    if(grepl("S$", x[length(x)])) {
      x[length(x)] <- as.integer(strsplit(x[length(x)], split = "S")[1])
    }) %>% 
    gsub('NULL', '0',.) %>% as.integer()
  
  seq_end <- qwidth(seqs_aligned) - end_AS_AC 
  
  cigar_width_df <- data.frame(seq_names = names(seqs_aligned),start_temp, seq_end, seq_end-start_temp, width(seqs_aligned),width(seqs_aligned)-(seq_end-start_temp))
  names(cigar_width_df) <- c("seq_names", "start","end","width","width_align","diff")
  
  ### B) getting the start of sequence which was aligned
  seq_start_pos <- names(seqs_aligned) %>% 
    gsub(":.*$", "",.) %>% 
    gsub("^.+_", "",.) %>% as.integer()
  
  seq_chr_name <- names(seqs_aligned) %>% 
    gsub("AS_", "",.) %>% 
    gsub("_[^_]*$", "",., perl = TRUE)
  
  # getting methylation site an normalizing to the aligned sequence part
  methyl_sites <- mcols(seqs_aligned)$methyl_pos %>% 
    strsplit(",") %>% lapply(function(x) as.integer(x))
  normalized_methyl_sites <- mapply(function(x, start_pos) x - start_pos, methyl_sites, seq_start_pos, SIMPLIFY = FALSE)
  mapped_methyl_sites <- mapply(function(x, start_pos) x - start_pos, normalized_methyl_sites, cigar_width_df$start, SIMPLIFY = FALSE)
  aligned_methyl_sites <- mapply(function(x, start_pos, align_pos) x - start_pos + align_pos, mapped_methyl_sites, cigar_width_df$start, start(seqs_aligned), SIMPLIFY = FALSE)
  
  methyl_on_aligned <- mapply(function(x, w) x <= w & x > 0, mapped_methyl_sites, width(seqs_aligned), SIMPLIFY = FALSE)
  
  final_methyl_sites_AS <- lapply(seq_along(mapped_methyl_sites), function(i) {
    df <- data.frame(seq_names = cigar_width_df$seq_names[i], 
                     pos_rgenome = methyl_sites[[i]],
                     pos_seq = mapped_methyl_sites[[i]],
                     pos_align = aligned_methyl_sites[[i]], 
                     width_align = width(seqs_aligned)[i], 
                     chr_align = seqnames(seqs_aligned)[i], 
                     aligned = methyl_on_aligned[[i]], 
                     seq_number = sprintf("seq_%03d", 1:length(aligned_methyl_sites))[i], 
                     species = species)
    df$Chr <- seq_chr_name[i]
    df <- cbind(df[length(df)], df[-length(df)])
    return(df)
  })
  df_final_methyl_sites_AS <- bind_rows(final_methyl_sites_AS)
  return(df_final_methyl_sites_AS)
}
#### getting all species transformed ####

### BT2
#AC
df_cigar_AC_bt2 <- cigar.to.width(overlap_AC_AC_bt2)
AC_methyl_df_bt2 <- map.methyl.to.align(overlap_AC_AC_bt2, df_cigar_AC_bt2, "AC")

#AS
df_cigar_AS_bt2 <- cigar.to.width(overlap_AC_AS_bt2)
AS_methyl_df_bt2 <- map.methyl.to.align(overlap_AC_AS_bt2, df_cigar_AS_bt2, "AS")

#EH
df_cigar_EH_bt2 <- cigar.to.width(overlap_AC_EH_bt2)
EH_methyl_df_bt2 <- map.methyl.to.align(overlap_AC_EH_bt2, df_cigar_EH_bt2, "EH")

#JM
df_cigar_JM_bt2 <- cigar.to.width(overlap_AC_JM_bt2)
JM_methyl_df_bt2 <- map.methyl.to.align(overlap_AC_JM_bt2, df_cigar_JM_bt2, "JM")

#ALL Bowtie2
combined_df_bt2 <- bind_rows(AC_methyl_df_bt2, AS_methyl_df_bt2, EH_methyl_df_bt2,JM_methyl_df_bt2)

AC_chr_names <- sort(unique(combined_df_bt2$chr_align))
AC_chr_names_simple <-  paste0("Chr ", 1:23)

combined_df_bt2 <- cbind(combined_df_bt2, data.frame(Chr_align_simple = AC_chr_names_simple[match(combined_df_bt2$chr_align, AC_chr_names)], align_method = "Bowtie2")) 

### Mini
#AC
df_cigar_AC_mini <- cigar.to.width(overlap_AC_AC_mini)
AC_methyl_df_mini <- map.methyl.to.align(overlap_AC_AC_mini, df_cigar_AC_mini, "AC")

#AS
df_cigar_AS_mini <- cigar.to.width(overlap_AC_AS_mini)
AS_methyl_df_mini <- map.methyl.to.align(overlap_AC_AS_mini, df_cigar_AS_mini, "AS")

#EH
df_cigar_EH_mini <- cigar.to.width(overlap_AC_EH_mini)
EH_methyl_df_mini <- map.methyl.to.align(overlap_AC_EH_mini, df_cigar_EH_mini, "EH")

#JM
df_cigar_JM_mini <- cigar.to.width(overlap_AC_JM_mini)
JM_methyl_df_mini <- map.methyl.to.align(overlap_AC_JM_mini, df_cigar_JM_mini, "JM")

#ALL Mini
combined_df_mini <- bind_rows(AC_methyl_df_mini, AS_methyl_df_mini, EH_methyl_df_mini,JM_methyl_df_mini)

AC_chr_names <- sort(unique(combined_df_mini$chr_align))
AC_chr_names_simple <-  paste0("Chr ", 1:23)

combined_df_mini <- cbind(combined_df_mini, data.frame(Chr_align_simple = AC_chr_names_simple[match(combined_df_mini$chr_align, AC_chr_names)], align_method ="Minimap")) 

### MINI & BOWTIE2
combined_df_bt2_mini <- rbind(combined_df_bt2, combined_df_mini)
combined_df_bt2_mini$Chr_align_simple <- factor(combined_df_bt2_mini$Chr_align_simple, levels = AC_chr_names_simple)
combined_df_bt2$Chr_align_simple <- factor(combined_df_bt2_mini$Chr_align_simple, levels = AC_chr_names_simple)
combined_df_mini$Chr_align_simple <- factor(combined_df_bt2_mini$Chr_align_simple, levels = AC_chr_names_simple)


#### plotting data ####
##bt2
library(ggplot2)
# library(ggpattern)

ggplot(combined_df_bt2, aes(x = pos_align, fill = species)) +
  geom_histogram(alpha = 0.5, position = "identity") +  
  facet_wrap(~ Chr_align_simple, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "All CpGs") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

ggplot(subset(combined_df_bt2, chr_align == "NC_044048.1"), aes(x = pos_align, fill = species)) +
  geom_histogram(alpha = 0.5, position = "identity") +  
  labs(x = "Position", y = "CpGs",title = "CpGs on Chromosome 1") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

ggplot(subset(combined_df_bt2, chr_align == "NC_044048.1" & aligned == TRUE), aes(x = pos_align, fill = species)) +
  geom_histogram(alpha = 0.5, position = "identity") +  
  labs(x = "Position", y = "CpGs",title = "CpGs on Chromosome 1 on SMRs") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

ggplot(subset(combined_df_bt2, aligned == TRUE), aes(x = pos_align, fill = species)) +
  geom_histogram(alpha = 0.5, position = "identity") +  
  facet_wrap(~ chr_align, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "All CpGs on SMRs") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

##combineds
ggplot(combined_df_bt2_mini, aes(x = pos_align, fill = align_method)) +
  geom_histogram(alpha = 0.5, position = "identity") +  
  facet_wrap(~ Chr_align_simple, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "All shared methylation regions for bowtie2 and minimap") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

ggplot(combined_df_bt2_mini, aes(x = pos_align, fill = species, shape = align_method)) +
  geom_histogram(alpha = 1) +  
  facet_wrap(~ Chr_align_simple, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "All shared methylation regions for bowtie2 and minimap") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

ggplot(subset(combined_df_bt2_mini, Chr_align_simple == "Chr 5"), aes(x = pos_align, fill = species)) +
  geom_histogram(alpha = 1, position = "dodge") +  
  labs(x = "Position", y = "CpGs",title = "CpGs on one Chromosome") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  



## mini only
ggplot(subset(combined_df_bt2_mini, align_method == "Minimap"), aes(x = pos_align, fill = species)) +
  geom_histogram(alpha = 0.5, position = "identity") +  
  facet_wrap(~ Chr_align_simple, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "All CpGs") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

# stacked
ggplot(subset(combined_df_bt2_mini, align_method == "Minimap"), aes(x = pos_align, fill = species)) +
  geom_histogram(alpha = 1) +  
  facet_wrap(~ Chr_align_simple, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "All CpGs") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

ggplot(subset(combined_df_bt2_mini, Chr_align_simple == "Chr 20"), aes(x = pos_align, fill = align_method)) +
  geom_histogram(alpha = 1, position = "dodge") +  
  labs(x = "Position", y = "CpGs",title = "CpGs on Chromosome 1") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  


## Density plots 
ggplot(combined_df_bt2_mini, aes(x = pos_align, fill = align_method, alpha = 0.5)) +
  geom_density() +  
  facet_wrap(~ Chr_align_simple, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "All CpGs") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))

ggplot(subset(combined_df_bt2_mini, Chr_align_simple == "Chr 20"), aes(x = pos_align, fill = align_method, alpha = 0.5)) +
  geom_density() +  
  labs(x = "Position", y = "CpGs",title = "CpGs on Chromosome 1") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

#### FIND OUT >>> why the width for some ranges vary from the calculated width based on the CIGAR code ANSWER: it was the inserts (add to width) and deletion (take from width)


#### Appendix ####

# # getting start of alignment of CIGAR
# cigar_sep <- cigar(overlap_AC_AS_bt2) %>% 
#   sapply(function(x) regmatches(x, gregexpr("\\d+[A-Z]", x)))
# 
# start_AS_AC <- lapply(cigar_sep, function(x) if(grepl("S$", x[1])) {
#   x[1] <- as.integer(strsplit(x[1], split = "S")[1])
# }) %>% 
#   gsub('NULL', '0',.) %>% as.integer()
# 
# # getting end of alignment of CIGAR
# end_AS_AC <- lapply(cigar_sep, function(x)
#   if(grepl("S$", x[length(x)])) {
#     x[length(x)] <- as.integer(strsplit(x[length(x)], split = "S")[1])
#   }) %>% 
#   gsub('NULL', '0',.) %>% as.integer()
# 
# seq_end <- qwidth(overlap_AC_AS_bt2) - end_AS_AC 
# 
# cigar_width_df <- data.frame(row.names = names(overlap_AC_AS_bt2),start_AS_AC, seq_end, seq_end-start_AS_AC, width(overlap_AC_AS_bt2),width(overlap_AC_AS_bt2)-(seq_end-start_AS_AC))
# names(cigar_width_df) <- c("start","end","width","width_align","diff")
# 
# # getting the start of sequence which was aligned
# seq_start_pos <- names(overlap_AC_AS_bt2) %>% 
#   gsub(":.*$", "",.) %>% 
#   gsub("^.+_", "",.) %>% as.integer()
# 
# seq_chr_name <- names(overlap_AC_AS_bt2) %>% 
#   gsub("AS_", "",.) %>% 
#   gsub("_[^_]*$", "",., perl = TRUE)
# 
# # getting methylation site and normalizing to the aligned sequence part
# methyl_sites <- mcols(overlap_AC_AS_bt2)$methyl_pos %>% 
#   strsplit(",") %>% lapply(function(x) as.integer(x))
# normalized_methyl_sites <- mapply(function(x, start_pos) x - start_pos, methyl_sites, seq_start_pos, SIMPLIFY = FALSE)
# mapped_methyl_sites <- mapply(function(x, start_pos) x - start_pos, normalized_methyl_sites, cigar_width_df$start, SIMPLIFY = FALSE)
# aligned_methyl_sites <- mapply(function(x, start_pos, align_pos) x - start_pos + align_pos, mapped_methyl_sites, cigar_width_df$start, start(overlap_AC_AS_bt2), SIMPLIFY = FALSE)
# 
# methyl_on_aligned <- mapply(function(x, w) x <= w & x > 0, mapped_methyl_sites, width(overlap_AC_AS_bt2), SIMPLIFY = FALSE)
# # giving the elements in the list the corresponding chromosome names 
# # Map(function(l, v) list(l, list(v)), mapped_methyl_sites, seq_chr_name)
# 
# # creating a list of dataframes to be able to identify the CpGs and where they were coming from
# final_methyl_sites_AS <- lapply(seq_along(mapped_methyl_sites), function(i) {
#   df <- data.frame(pos_rgenome = methyl_sites[[i]], 
#                    pos_seq = mapped_methyl_sites[[i]], 
#                    pos_align = aligned_methyl_sites[[i]], 
#                    width_align = width(overlap_AC_AS_bt2)[i], 
#                    chr_align = seqnames(overlap_AC_AS_bt2)[i], 
#                    aligned = methyl_on_aligned[[i]], 
#                    seq_number = sprintf("seq_%03d", 1:length(aligned_methyl_sites))[i])
#   df$Chr <- seq_chr_name[i]
#   df <- cbind(df[length(df)], df[-length(df)])
#   return(df)
# })
# 
# df_final_methyl_sites_AS <- bind_rows(final_methyl_sites_AS)
# df_final <- df_final_methyl_sites_AS
