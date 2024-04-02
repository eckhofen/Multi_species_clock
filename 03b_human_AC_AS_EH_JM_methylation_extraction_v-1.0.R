#### Overview ####
# Extraction of methylation data based on overlapping sequences

#### Settings ####
setwd("/powerplant/workspace/cfngle")

#### Preparation ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(tidyverse)

## loading data 

load("results-data/02_conserved_seq/HS_AC_AS_EH_JM_overlaps_bt2.Rdata")


HS_overlap_seqs_bt2


overlap_HS_AC_bt2 <- HS_overlap_seqs_bt2[[1]]
overlap_HS_AS_bt2 <- HS_overlap_seqs_bt2[[2]]
overlap_HS_EH_bt2 <- HS_overlap_seqs_bt2[[3]]
overlap_HS_JM_bt2 <- HS_overlap_seqs_bt2[[4]]

## getting SMR 
HS_gr_overlap_seqs_bt2 <- lapply(HS_overlap_seqs_bt2, function(x) granges(x))

HS_group_gr_overlap_bt2 <-  c(HS_gr_overlap_seqs_bt2[[1]], HS_gr_overlap_seqs_bt2[[2]], HS_gr_overlap_seqs_bt2[[3]], HS_gr_overlap_seqs_bt2[[4]])
HS_SMR_b_bt2 <- GenomicRanges::reduce(HS_group_gr_overlap_bt2)
names(HS_SMR_b_bt2) <- sprintf("HS_SMR_b_bt2_%03d", 1:length(HS_SMR_b_bt2))

HS_AC_AS_EH_JM_SMR_b_bt2 <- HS_SMR_b_bt2

save(HS_AC_AS_EH_JM_SMR_b_bt2, file = "/workspace/cfngle/results-data/04_SMRs/HS_AC_AS_EH_JM_SMR_b_bt2.Rdata")

### master function

get.methyl.sites <- function(seqs_aligned, species = "undefined", SMRs = "undefined") {
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
  
  ### C) finding overlaps between sequences and SMRs
  SMR_index <- subjectHits(findOverlaps(seqs_aligned, SMRs))
  
  
  final_methyl_sites_AS <- lapply(seq_along(mapped_methyl_sites), function(i) {
    df <- data.frame(seq_names = cigar_width_df$seq_names[i], 
                     pos_rgenome = methyl_sites[[i]],
                     pos_seq = mapped_methyl_sites[[i]],
                     pos_align = aligned_methyl_sites[[i]], 
                     width_align = width(seqs_aligned)[i], 
                     chr_align = seqnames(seqs_aligned)[i], 
                     aligned = methyl_on_aligned[[i]], 
                     seq_number = sprintf("seq_%03d", 1:length(aligned_methyl_sites))[i],
                     SMR = sprintf("SMR_%03d", SMR_index[i]),
                     species = species)
    df$Chr <- seq_chr_name[i]
    df <- cbind(df[length(df)], df[-length(df)])
    return(df)
  })
  df_final_methyl_sites_AS <- bind_rows(final_methyl_sites_AS)
  return(df_final_methyl_sites_AS)
}
#### getting all species transformed ####
#ALL Bowtie2

AC_methyl_df_bt2 <- get.methyl.sites(overlap_HS_AC_bt2, species = "AC", SMRs = HS_SMR_b_bt2)
AS_methyl_df_bt2 <- get.methyl.sites(overlap_HS_AS_bt2, species = "AS", SMRs = HS_SMR_b_bt2)
EH_methyl_df_bt2 <- get.methyl.sites(overlap_HS_EH_bt2, species = "EH", SMRs = HS_SMR_b_bt2)
JM_methyl_df_bt2 <- get.methyl.sites(overlap_HS_JM_bt2, species = "JM", SMRs = HS_SMR_b_bt2)

combined_df_bt2 <- bind_rows(AC_methyl_df_bt2, AS_methyl_df_bt2, EH_methyl_df_bt2,JM_methyl_df_bt2)

HS_chr_names <- sort(unique(combined_df_bt2$chr_align))
# AC_chr_names_simple <-  paste0("Chr ", 1:23)

combined_df_bt2 <- cbind(combined_df_bt2, data.frame(align_method = "Bowtie2")) 

#### plotting data ####
##bt2
library(ggplot2)
# library(ggpattern)

ggplot(combined_df_bt2, aes(x = pos_align, fill = species)) +
  geom_histogram(alpha = 0.5, position = "identity") +  
  facet_wrap(~ chr_align, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "All CpGs") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0))  

#### get methyl sitesand value ####
get.methyl.sites <- function(seqs_aligned, species = "undefined", SMRs = "undefined") {
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
  
  ### C) finding overlaps between sequences and SMRs
  SMR_index <- subjectHits(findOverlaps(seqs_aligned, SMRs))
  
  
  final_methyl_sites_AS <- lapply(seq_along(mapped_methyl_sites), function(i) {
    df <- data.frame(seq_names = cigar_width_df$seq_names[i], 
                     pos_rgenome = methyl_sites[[i]],
                     pos_seq = mapped_methyl_sites[[i]],
                     pos_align = aligned_methyl_sites[[i]], 
                     width_align = width(seqs_aligned)[i], 
                     chr_align = seqnames(seqs_aligned)[i], 
                     aligned = methyl_on_aligned[[i]], 
                     seq_number = sprintf("seq_%03d", 1:length(aligned_methyl_sites))[i],
                     SMR = sprintf("SMR_%03d", SMR_index[i]),
                     species = species)
    df$Chr <- seq_chr_name[i]
    df <- cbind(df[length(df)], df[-length(df)])
    return(df)
  })
  df_final_methyl_sites_AS <- bind_rows(final_methyl_sites_AS)
  return(df_final_methyl_sites_AS)
}

# getting methyl sites for all species

AC_methyl_sites <- get.methyl.sites(overlap_HS_AC_bt2, species = "AC", SMRs = HS_SMR_b_bt2)
AS_methyl_sites <- get.methyl.sites(overlap_HS_AS_bt2, species = "AS", SMRs = HS_SMR_b_bt2)
EH_methyl_sites <- get.methyl.sites(overlap_HS_EH_bt2, species = "EH", SMRs = HS_SMR_b_bt2)
JM_methyl_sites <- get.methyl.sites(overlap_HS_JM_bt2, species = "JM", SMRs = HS_SMR_b_bt2)

### load methylation data for samples as well as age vector for each species

##AC
# xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/meth-corrected-batchcorrected-cod.Rdata")
xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/Meth-complete-nobatchcorrection-cod.RData")
assign("AC_meth_data", get(xx))
AC_meth_data <- as.data.frame(AC_meth_data)
# tail(colnames(AC_meth_data))
AC_age <- AC_meth_data$age

##AS
xx <- load("/workspace/cfngle/raw-data/AS/zzz_methyl_data/Meth-complete-snapper.RData")
assign("AS_meth_data", get(xx))
# tail(colnames(AS_meth_data))
AS_age <- AS_meth_data$age

##EH
xx <- load("/workspace/cfngle/raw-data/EH/zzz-methyl_data/Meth-complete-hake.RData")
assign("EH_meth_data", get(xx))
# tail(colnames(EH_meth_data))
EH_age <- EH_meth_data$age

##JM
# JM_meth_data <- read.csv("/workspace/cfngle/raw-data/JM/zzz-methyldata/00_JM_methyldata_66079_CpGs.csv")
JM_meth_data <- load("/workspace/cfngle/raw-data/JM/zzz-methyldata/00_JM_methyldata_179818_CpGs.Rdata")
JM_meth_data <- JM_24_methyl_data

JM_age <- JM_meth_data$age

#### extract methylation data for all samples ####
# Not all the datasets have the same naming structure, hence the steps are different and are done one by one

##AC
meth_sites_names_tmp <- paste0(AC_methyl_sites$Chr, ".", AC_methyl_sites$pos_align) %>% 
  gsub("AC_", "Chr", .)
AC_meth_data_test <- gsub("X", "Chr", colnames(AC_meth_data))
table(meth_sites_names_tmp %in% AC_meth_data_test)

meth_columns_tmp <- sapply(meth_sites_names_tmp[meth_sites_names_tmp %in% AC_meth_data_test], function(x) grep(x, AC_meth_data_test)) %>% 
  as.vector()

AC_meth_values <- AC_meth_data[,meth_columns_tmp]

##AS
meth_sites_names_tmp <- paste0(AS_methyl_sites$Chr, "-", AS_methyl_sites$pos_rgenome)
meth_sites_names_tmp %in% colnames(AS_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(AS_meth_data))) %>% 
  as.vector()

AS_meth_values <- AS_meth_data[,meth_columns_tmp] 

# EH_meth_values_t <- t(EH_meth_data[,meth_columns_tmp]) %>% 
#   cbind(., data.frame(SMR = EH_methyl_sites$SMR))

##EH
meth_sites_names_tmp <- paste0(gsub("EH_", "",EH_methyl_sites$Chr), ".", EH_methyl_sites$pos_rgenome)
meth_sites_names_tmp %in% colnames(EH_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(EH_meth_data))) %>% 
  as.vector()

EH_meth_values <- EH_meth_data[,meth_columns_tmp] 

# EH_meth_values_t <- t(EH_meth_data[,meth_columns_tmp]) %>% 
#   cbind(., data.frame(SMR = EH_methyl_sites$SMR))

##JM
meth_sites_names_tmp <- paste0(gsub("JM_", "",JM_methyl_sites$Chr), ":", JM_methyl_sites$pos_rgenome)
meth_sites_names_tmp %in% colnames(JM_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(JM_meth_data))) %>% 
  as.vector()

JM_meth_values <- JM_meth_data[,meth_columns_tmp] 

# JM_meth_values_t <- t(JM_meth_data[,meth_columns_tmp]) %>% 
#   cbind(., data.frame(SMR = JM_methyl_sites$SMR))
