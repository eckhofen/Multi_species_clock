# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: Methylation extraction
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026-03

#### Overview ####
# Extraction of methylation data based on overlapping sequences

#### Settings ####
data_folder <- "01_data/"
save_folder <- "01_data/03_intermediate/03_SMR/"
results_folder <- "03_results/01_figures/"

raw_private_folder <- paste0(data_folder, "01_raw_private/")
annotation_folder <- paste0(data_folder, "03_intermediate/08_annotation/")
comparison_folder <- paste0(data_folder, "03_intermediate/07_data_comparison/")

dir.create(annotation_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(comparison_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

load("01_data/04_metadata/color_palettes.RData")

#### Preparation ####
library(GenomicRanges)
library(Biostrings)
library(GenomicAlignments)
library(tidyverse)
library(ggforce)
library(ggfortify)

## loading data 
load("01_data/02_external_server_outputs/02_conserved_seq/HS_AC_AS_EH_ZF_overlaps_bt2.Rdata")

# assign overlapping sequences for each species
overlap_HS_AC <- HS_overlap_seqs_bt2[[1]]
overlap_HS_AS <- HS_overlap_seqs_bt2[[2]]
overlap_HS_EH <- HS_overlap_seqs_bt2[[3]]
overlap_HS_ZF <- HS_overlap_seqs_bt2[[4]]

# transforming aligned reads into GRanges object
HS_gr_overlap_seqs <- lapply(HS_overlap_seqs_bt2, function(x) granges(x))

# using the overlap of the sequences to get the shared methylation regions (SMRs)
HS_group_gr_overlap <- do.call(c, HS_gr_overlap_seqs)
HS_SMR_b <- GenomicRanges::reduce(HS_group_gr_overlap)
names(HS_SMR_b) <- sprintf("HS_SMR_b_bt2_%03d", seq_along(HS_SMR_b))

# Saving file
save(HS_SMR_b, file = paste0(save_folder, "HS_SMR.RData"))

#### Methylation extraction ####

## function for extraction
get.methyl.sites <- function(seqs_aligned, species = "undefined", SMRs = "undefined") {
  ### A) extract information from CIGAR code
  cigar_sep <- cigar(seqs_aligned) %>% 
    sapply(function(x) regmatches(x, gregexpr("\\d+[A-Z]", x)))
  
  start_temp <- lapply(cigar_sep, function(x) if(grepl("S$", x[1])) {
    x[1] <- as.integer(strsplit(x[1], split = "S")[1])
  }) %>% 
    gsub('NULL', '0',.) %>% as.integer()
  
  # getting end of alignment from CIGAR
  end <- lapply(cigar_sep, function(x)
    if(grepl("S$", x[length(x)])) {
      x[length(x)] <- as.integer(strsplit(x[length(x)], split = "S")[1])
    }) %>% 
    gsub('NULL', '0',.) %>% as.integer()
  
  seq_end <- qwidth(seqs_aligned) - end 
  
  cigar_width_df <- data.frame(seq_names = names(seqs_aligned),
                               start_temp, 
                               seq_end, 
                               seq_end-start_temp, 
                               width(seqs_aligned),
                               width(seqs_aligned)-(seq_end-start_temp))
  names(cigar_width_df) <- c("seq_names", "start","end","width","width_align","diff")
  
  ### B) getting the start of sequence which was aligned
  seq_start_pos <- names(seqs_aligned) %>% 
    gsub(":.*$", "",.) %>% 
    gsub("^.+_", "",.) %>% as.integer()
  
  seq_chr_name <- names(seqs_aligned) %>% 
    gsub("AS_", "",.) %>% 
    gsub("_[^_]*$", "",., perl = TRUE)
  
  # getting methylation site and normalizing to the aligned sequence
  methyl_sites <- mcols(seqs_aligned)$methyl_pos %>% 
    strsplit(",") %>% lapply(function(x) as.integer(x))
  normalized_methyl_sites <- mapply(function(x, start_pos) x - start_pos, methyl_sites, seq_start_pos, SIMPLIFY = FALSE)
  mapped_methyl_sites <- mapply(function(x, start_pos) x - start_pos, normalized_methyl_sites, cigar_width_df$start, SIMPLIFY = FALSE)
  aligned_methyl_sites <- mapply(function(x, start_pos, align_pos) x - start_pos + align_pos, mapped_methyl_sites, cigar_width_df$start, start(seqs_aligned), SIMPLIFY = FALSE)
  
  methyl_on_aligned <- mapply(function(x, w) x <= w & x > 0, mapped_methyl_sites, width(seqs_aligned), SIMPLIFY = FALSE)
  
  ### C) finding overlaps between sequences and SMRs
  SMR_index <- subjectHits(findOverlaps(seqs_aligned, SMRs))
  
  final_methyl_sites <- lapply(seq_along(mapped_methyl_sites), function(i) {
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
  df_final_methyl_sites <- bind_rows(final_methyl_sites)
  return(df_final_methyl_sites)
}
#### Getting all species transformed ####

AC_methyl_sites <- get.methyl.sites(overlap_HS_AC, species = "AC", SMRs = HS_SMR_b)
AS_methyl_sites <- get.methyl.sites(overlap_HS_AS, species = "AS", SMRs = HS_SMR_b)
EH_methyl_sites <- get.methyl.sites(overlap_HS_EH, species = "EH", SMRs = HS_SMR_b)
ZF_methyl_sites <- get.methyl.sites(overlap_HS_ZF, species = "ZF", SMRs = HS_SMR_b)

methyl_sites_combined <- bind_rows(AC_methyl_sites, AS_methyl_sites, EH_methyl_sites, ZF_methyl_sites)
methyl_sites_combined$species <- as.factor(methyl_sites_combined$species)
HS_chr_names <- sort(unique(methyl_sites_combined$chr_align))

if ("SMR" %in% colnames(methyl_sites_combined) && !("SMR" %in% colnames(methyl_sites_combined))) {
  methyl_sites_combined <- dplyr::rename(methyl_sites_combined, SMR = SMR)
}

save(methyl_sites_combined, file = paste0(annotation_folder, "methylsites_all.RData"))

#### Methylation metadata ####
message("methylation metadata...")
## AC
xx <- load(paste0(raw_private_folder, "AC_methyldata.RData"))
assign("AC_meth_data", get(xx))
AC_meth_data <- as.data.frame(AC_meth_data)
AC_age <- AC_meth_data$age
AC_age[AC_age == 0] <- 0.01

## AS
xx <- load(paste0(raw_private_folder, "AS_methyldata.RData"))
assign("AS_meth_data", get(xx))
AS_age <- AS_meth_data$age

## EH
xx <- load(paste0(raw_private_folder, "EH_methyldata.RData"))
assign("EH_meth_data", get(xx))
EH_metadata_samples <- read.csv(paste0(raw_private_folder, "EH_meta_data.txt"), sep = "\t")
EH_age <- EH_metadata_samples$age

## ZF
ZF_meth_data <- load(paste0(raw_private_folder, "ZF_methyldata.RData"))
ZF_meth_data <- ZF_methyl_data

ZF_age <- ZF_meth_data$age/52

message("methylation metadata done")

#### Extract methylation data for all samples ####
message("extracting methylation data...") 
# Not all the datasets have the same naming structure, hence the steps are different and are done one by one

## AC
meth_sites_names_tmp_AC <- paste0(AC_methyl_sites$Chr, ".", AC_methyl_sites$pos_rgenome) %>% 
  gsub("AC_", "Chr", .)
AC_meth_data_test <- gsub("X", "Chr", colnames(AC_meth_data))
table(meth_sites_names_tmp_AC %in% AC_meth_data_test)
meth_columns_tmp <- sapply(meth_sites_names_tmp_AC[meth_sites_names_tmp_AC %in% AC_meth_data_test], function(x) grep(x, AC_meth_data_test)) %>% 
  as.vector()

AC_meth_values <- AC_meth_data[ ,meth_columns_tmp]
AC_methyl_sites <- AC_methyl_sites[meth_sites_names_tmp_AC %in% AC_meth_data_test,]

## AS
meth_sites_names_tmp <- paste0(AS_methyl_sites$Chr, "-", AS_methyl_sites$pos_rgenome)
meth_sites_names_tmp %in% colnames(AS_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(AS_meth_data))) %>% 
  as.vector()

AS_meth_values <- AS_meth_data[,meth_columns_tmp] 

## EH
meth_sites_names_tmp <- paste0(gsub("EH_", "",EH_methyl_sites$Chr), ".", EH_methyl_sites$pos_rgenome)
meth_sites_names_tmp %in% colnames(EH_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(EH_meth_data))) %>% 
  as.vector()

# EH_meth_values <- EH_meth_data[,meth_columns_tmp]
EH_meth_values <- EH_meth_data[meth_sites_names_tmp]

## ZF
meth_sites_names_tmp <- paste0(gsub("ZF_", "",ZF_methyl_sites$Chr), ":", ZF_methyl_sites$pos_rgenome)
table(meth_sites_names_tmp %in% colnames(ZF_meth_data))
meth_columns_tmp <- ZF_meth_data[meth_sites_names_tmp] %>% as.vector()
meth_columns_tmp <- unlist(meth_columns_tmp)
# ZF_meth_values_JM <- ZF_meth_data[,meth_columns_tmp] 
ZF_meth_values <- ZF_meth_data[meth_sites_names_tmp]

## saving data
save_folder <- paste0(data_folder, "03_intermediate/04_methyl_values/")
dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)

message("extracting methylation data done")

### saving methylation VALUES
message("saving methylation data...")
## AC
write.csv(AC_meth_values, file = paste0(save_folder, "HS_AC_meth_values.csv") )
save(AC_meth_values, file = paste0(save_folder, "HS_AC_meth_values.Rdata"))

## AS
write.csv(AS_meth_values, file = paste0(save_folder, "HS_AS_meth_values.csv") )
save(AS_meth_values, file = paste0(save_folder, "HS_AS_meth_values.Rdata"))

## EH
write.csv(EH_meth_values, file = paste0(save_folder, "HS_EH_meth_values.csv") )
save(EH_meth_values, file = paste0(save_folder, "HS_EH_meth_values.Rdata"))

## ZF
# imputing ZF missing values
library(zoo)
set.seed(123)
ZF_meth_values <- na.aggregate(ZF_meth_values)

write.csv(ZF_meth_values, file = paste0(save_folder, "HS_ZF_meth_values.csv"))
save(ZF_meth_values, file = paste0(save_folder, "HS_ZF_meth_values.Rdata"))

# all 
save(AC_meth_values, AS_meth_values, EH_meth_values, ZF_meth_values, file = paste0(save_folder, "all_meth_values.Rdata"))

message("methylation values saved")

### saving methylation SITES
## AC
write.csv(AC_methyl_sites, file = paste0(save_folder, "HS_AC_methyl_sites.csv") )
save(AC_methyl_sites, file = paste0(save_folder, "HS_AC_methyl_sites.Rdata"))

## AS
write.csv(AS_methyl_sites, file = paste0(save_folder, "HS_AS_methyl_sites.csv") )
save(AS_methyl_sites, file = paste0(save_folder, "HS_AS_methyl_sites.Rdata"))

## EH
write.csv(EH_methyl_sites, file = paste0(save_folder, "HS_EH_methyl_sites.csv") )
save(EH_methyl_sites, file = paste0(save_folder, "HS_EH_methyl_sites.Rdata"))

## ZF
write.csv(ZF_methyl_sites, file = paste0(save_folder, "HS_ZF_methyl_sites.csv") )
save(ZF_methyl_sites, file = paste0(save_folder, "HS_ZF_methyl_sites.Rdata"))

# all 
save(AC_methyl_sites, AS_methyl_sites, EH_methyl_sites, ZF_methyl_sites, file = paste0(save_folder, "all_methyl_sites.Rdata"))

message("methylation sites saved")

### saving age metadata
save(AC_age, AS_age, EH_age, ZF_age, file = paste0(save_folder, "HS_all_age.Rdata"))

message("age metadata saved")

message("methylation extraction completed!")