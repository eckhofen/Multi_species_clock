#### Overview ####

#Reducing overlapping ranges of all species and extract the shared (methylation) regions 
# adding the methylation data to the extracted methylation sites 

#### Preparation ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(GenomicAlignments)
library(tidyverse)
library(ggforce)

#### Settings ####
setwd("/workspace/cfngle/results-data/02_conserved_seq")

#### Data manipulation ####
# loading overlapping sequences for all species (only bt2 in this case) 
load("AC_AS_EH_JM_ZF_overlaps_bt2_v_2.RData")

load("JM_overlaps_bt2.RData")
overlap_seqs_bt2
overlap_seqs_bt2_JM

# converting GAlignments object to GRanges to be able to reduce them 
gr_overlap_seqs_bt2 <- lapply(overlap_seqs_bt2, function(x) granges(x))
gr_overlap_seqs_bt2_JM <- lapply(overlap_seqs_bt2_JM, function(x) granges(x))

### Shared methylation regions
## OPTION A)

# intersecting all sequences (option a!) This means that only the sequences which are intersecting will be selected
intersected_seqs_bt2 <- GenomicRanges::intersect(gr_overlap_seqs_bt2[[1]], gr_overlap_seqs_bt2[[2]]) %>% 
  GenomicRanges::intersect(gr_overlap_seqs_bt2[[3]])
boxplot(width(intersected_seqs_bt2))
names(intersected_seqs_bt2) <- sprintf("SMR_c_bt2_%03d", 1:length(intersected_seqs_bt2))
SMR_a_bt2 <- intersected_seqs_bt2

## OPTION B)
# reducing them is option b of how to define the shared methyl regions
# group them first
group_gr_overlap_bt2 <- c(gr_overlap_seqs_bt2[[1]], gr_overlap_seqs_bt2[[2]], gr_overlap_seqs_bt2[[3]], gr_overlap_seqs_bt2[[4]],gr_overlap_seqs_bt2[[5]])
SMR_b_bt2 <- GenomicRanges::reduce(group_gr_overlap_bt2)
names(SMR_b_bt2) <- sprintf("SMR_b_bt2_%03d", 1:length(SMR_b_bt2))

all_SMR_b_bt2 <- SMR_b_bt2

save(all_SMR_b_bt2, file = "/workspace/cfngle/results-data/04_SMRs/all_SMR_b_bt2_v_2.Rdata")

### FOR JM
group_gr_overlap_bt2_JM <- c(gr_overlap_seqs_bt2_JM[[1]], gr_overlap_seqs_bt2_JM[[2]], gr_overlap_seqs_bt2_JM[[3]], gr_overlap_seqs_bt2_JM[[4]], gr_overlap_seqs_bt2_JM[[5]])
SMR_b_bt2_JM <- GenomicRanges::reduce(group_gr_overlap_bt2_JM)
names(SMR_b_bt2_JM) <- sprintf("SMR_b_bt2_JM_%03d", 1:length(SMR_b_bt2_JM))

SMR_b_bt2_JM

save(SMR_b_bt2_JM, file = "/workspace/cfngle/results-data/04_SMRs/SMR_b_bt2_JM.RData")

#### Getting methylation values ####
### Loading sequences

## getting overlapped sequence data
load("/workspace/cfngle/results-data/02_conserved_seq/AC_AS_EH_JM_ZF_overlaps_bt2_v_2.RData")
# load("/workspace/cfngle/results-data/02_conserved_seq/AC_AS_EH_JM_overlaps_mini.R")

overlap_seqs_bt2
# overlap_seqs_mini

## BT2
overlap_AC_AC_bt2 <- overlap_seqs_bt2[[1]]
overlap_AC_AS_bt2 <- overlap_seqs_bt2[[2]]
overlap_AC_EH_bt2 <- overlap_seqs_bt2[[3]]
overlap_AC_JM_bt2 <- overlap_seqs_bt2[[4]]
overlap_AC_ZF_bt2 <- overlap_seqs_bt2[[5]]

## BT2
overlap_JM_AC_bt2 <- overlap_seqs_bt2_JM[[1]]
overlap_JM_AS_bt2 <- overlap_seqs_bt2_JM[[2]]
overlap_JM_EH_bt2 <- overlap_seqs_bt2_JM[[3]]
overlap_JM_JM_bt2 <- overlap_seqs_bt2_JM[[4]]
overlap_JM_ZF_bt2 <- overlap_seqs_bt2_JM[[5]]


## function to get methyl sites

get.methyl.sites <- function(seqs_aligned, species = "undefined", SMRs = "undefined") {
  ### A) extract information from CIGAR code
  cigar_sep <- cigar(seqs_aligned) %>% 
    sapply(function(x) regmatches(x, gregexpr("\\d+[A-Z]", x)))
  
  start_temp <- lapply(cigar_sep, function(x) if(grepl("S$", x[1])) {
    x[1] <- as.integer(strsplit(x[1], split = "S")[1])
  }) %>% 
    gsub('NULL', '0',.) %>% as.integer()
  
  # getting end of alignment of CIGAR
  end_ZF_AC <- lapply(cigar_sep, function(x)
    if(grepl("S$", x[length(x)])) {
      x[length(x)] <- as.integer(strsplit(x[length(x)], split = "S")[1])
    }) %>% 
    gsub('NULL', '0',.) %>% as.integer()
  
  seq_end <- qwidth(seqs_aligned) - end_ZF_AC 
  
  cigar_width_df <- data.frame(seq_names = names(seqs_aligned),start_temp, seq_end, seq_end-start_temp, width(seqs_aligned),width(seqs_aligned)-(seq_end-start_temp))
  names(cigar_width_df) <- c("seq_names", "start","end","width","width_align","diff")
  
  ### B) getting the start of sequence which was aligned
  seq_start_pos <- names(seqs_aligned) %>% 
    gsub(":.*$", "",.) %>% 
    gsub("^.+_", "",.) %>% as.integer()
  
  seq_chr_name <- names(seqs_aligned) %>% 
    gsub("ZF_", "",.) %>% 
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
  
  
  final_methyl_sites_ZF <- lapply(seq_along(mapped_methyl_sites), function(i) {
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
  df_final_methyl_sites_ZF <- bind_rows(final_methyl_sites_ZF)
  return(df_final_methyl_sites_ZF)
}

# getting methyl sites for all species

AC_methyl_sites <- get.methyl.sites(overlap_AC_AC_bt2, species = "AC", SMRs = SMR_b_bt2)
AS_methyl_sites <- get.methyl.sites(overlap_AC_AS_bt2, species = "AS", SMRs = SMR_b_bt2)
EH_methyl_sites <- get.methyl.sites(overlap_AC_EH_bt2, species = "EH", SMRs = SMR_b_bt2)
JM_methyl_sites <- get.methyl.sites(overlap_AC_JM_bt2, species = "JM", SMRs = SMR_b_bt2)
ZF_methyl_sites <- get.methyl.sites(overlap_AC_ZF_bt2, species = "ZF", SMRs = SMR_b_bt2)

AC_methyl_sites_JM <- get.methyl.sites(overlap_JM_AC_bt2, species = "AC", SMRs = SMR_b_bt2_JM)
AS_methyl_sites_JM <- get.methyl.sites(overlap_JM_AS_bt2, species = "AS", SMRs = SMR_b_bt2_JM)
EH_methyl_sites_JM <- get.methyl.sites(overlap_JM_EH_bt2, species = "EH", SMRs = SMR_b_bt2_JM)
JM_methyl_sites_JM <- get.methyl.sites(overlap_JM_JM_bt2, species = "JM", SMRs = SMR_b_bt2_JM)
ZF_methyl_sites_JM <- get.methyl.sites(overlap_JM_ZF_bt2, species = "ZF", SMRs = SMR_b_bt2_JM)

### load methylation data for samples as well as age vector for each species

##AC
xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/meth-corrected-batchcorrected-cod.Rdata")
# xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/Meth-complete-nobatchcorrection-cod.RData")
assign("AC_meth_data", get(xx))
AC_meth_data <- as.data.frame(AC_meth_data)
# tail(colnames(AC_meth_data))
AC_age <- AC_meth_data$age

##AS
xx <- load("/workspace/cfngle/raw-data/AS/zzz_methyl_data/Meth-complete-snapper.RData")
assign("AS_meth_data", get(xx))
# tail(colnames(ZF_meth_data))
AS_age <- AS_meth_data$age

##EH
xx <- load("/workspace/cfngle/raw-data/EH/zzz-methyl_data/Meth-complete-hake.RData")
assign("EH_meth_data", get(xx))
tail(colnames(EH_meth_data))
EH_age <- EH_meth_data$age

##JM
# JM_meth_data <- read.csv("/workspace/cfngle/raw-data/JM/zzz-methyldata/00_JM_methyldata_66079_CpGs.csv")
JM_meth_data <- load("/workspace/cfngle/raw-data/JM/zzz-methyldata/00_JM_methyldata_243285_CpGs.Rdata")
JM_meth_data <- JM_24_methyl_data

JM_age <- JM_meth_data$age

##ZF
ZF_meth_data <- load("/workspace/cfngle/raw-data/ZF/zzz_methyldata/ZF_methyldata_88.RData")
ZF_meth_data <- ZF_methyl_data
tail(colnames(ZF_meth_data))

ZF_age <- ZF_meth_data$age

#### extract methylation data for all samples ####
# Not all the datasets have the same naming structure, hence the steps are different and are done one by one

### FOR AC rgenome
##AC
meth_sites_names_tmp_AC <- paste0(AC_methyl_sites$Chr, ".", AC_methyl_sites$pos_align) %>% 
  gsub("AC_", "Chr", .)
AC_meth_data_test <- gsub("X", "Chr", colnames(AC_meth_data))
table(meth_sites_names_tmp_AC %in% AC_meth_data_test)
meth_columns_tmp <- sapply(meth_sites_names_tmp_AC[meth_sites_names_tmp_AC %in% AC_meth_data_test], function(x) grep(x, AC_meth_data_test)) %>% 
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

##ZF
meth_sites_names_tmp <- paste0(gsub("ZF_", "",ZF_methyl_sites$Chr), ":", ZF_methyl_sites$pos_rgenome)
meth_sites_names_tmp %in% colnames(ZF_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(ZF_meth_data))) %>% 
  as.vector()
ZF_meth_values <- ZF_meth_data[,meth_columns_tmp] 

### saving data
save_dir <- "/workspace/cfngle/results-data/05_shared_methyl_values/"

##AC
write.csv(AC_meth_values, file = paste0(save_dir, "AC_meth_values_v_2.csv") )
save(AC_meth_values, file = paste0(save_dir, "AC_meth_values_v_2.Rdata"))

##AS
write.csv(AS_meth_values, file = paste0(save_dir, "AS_meth_values_v_2.csv") )
save(AS_meth_values, file = paste0(save_dir, "AS_meth_values_v_2.Rdata"))

##EH
write.csv(EH_meth_values, file = paste0(save_dir, "EH_meth_values_v_2.csv") )
save(EH_meth_values, file = paste0(save_dir, "EH_meth_values_v_2.Rdata"))

##JM
write.csv(JM_meth_values, file = paste0(save_dir, "JM_meth_values_v_2.csv") )
save(JM_meth_values, file = paste0(save_dir, "JM_meth_values_v_2.Rdata"))

##ZF
write.csv(ZF_meth_values, file = paste0(save_dir, "ZF_meth_values_v_2.csv") )
save(ZF_meth_values, file = paste0(save_dir, "ZF_meth_values_v_2.Rdata"))

### FOR JM rgenome

##AC
meth_sites_names_tmp_AC <- paste0(AC_methyl_sites_JM$Chr, ".", AC_methyl_sites_JM$pos_rgenome) %>% 
  gsub("AC_", "Chr", .)
AC_meth_data_test <- gsub("X", "Chr", colnames(AC_meth_data))
table(meth_sites_names_tmp_AC %in% AC_meth_data_test)

meth_columns_tmp <- sapply(meth_sites_names_tmp_AC[meth_sites_names_tmp_AC %in% AC_meth_data_test], function(x) grep(x, AC_meth_data_test)) %>% 
  as.vector()

AC_meth_values_JM <- AC_meth_data[,meth_columns_tmp]

##AS
meth_sites_names_tmp <- paste0(AS_methyl_sites_JM$Chr, "-", AS_methyl_sites_JM$pos_rgenome)
meth_sites_names_tmp %in% colnames(AS_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(AS_meth_data))) %>% 
  as.vector()

AS_meth_values_JM <- AS_meth_data[,meth_columns_tmp] 

# EH_meth_values_t <- t(EH_meth_data[,meth_columns_tmp]) %>% 
#   cbind(., data.frame(SMR = EH_methyl_sites_JM$SMR))

##EH
meth_sites_names_tmp <- paste0(gsub("EH_", "",EH_methyl_sites_JM$Chr), ".", EH_methyl_sites_JM$pos_rgenome)
meth_sites_names_tmp %in% colnames(EH_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(EH_meth_data))) %>% 
  as.vector()

EH_meth_values_JM <- EH_meth_data[,meth_columns_tmp] 

# EH_meth_values_t <- t(EH_meth_data[,meth_columns_tmp]) %>% 
#   cbind(., data.frame(SMR = EH_methyl_sites_JM$SMR))

##JM
meth_sites_names_tmp <- paste0(gsub("JM_", "",JM_methyl_sites_JM$Chr), ":", JM_methyl_sites_JM$pos_rgenome)
meth_sites_names_tmp %in% colnames(JM_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(JM_meth_data))) %>% 
  as.vector()

JM_meth_values_JM <- JM_meth_data[,meth_columns_tmp] 

##ZF
meth_sites_names_tmp <- paste0(gsub("ZF_", "",ZF_methyl_sites_JM$Chr), ":", ZF_methyl_sites_JM$pos_rgenome)
table(meth_sites_names_tmp %in% colnames(ZF_meth_data))
meth_columns_tmp <- ZF_meth_data[meth_sites_names_tmp] %>% as.vector()
meth_columns_tmp <- unlist(meth_columns_tmp)
ZF_meth_values_JM <- ZF_meth_data[,meth_columns_tmp] 
ZF_meth_values_JM <- ZF_meth_data[meth_sites_names_tmp]

### saving data
save_dir <- "/workspace/cfngle/results-data/05_shared_methyl_values/"

##AC
write.csv(AC_meth_values_JM, file = paste0(save_dir, "AC_meth_values_JM.csv") )
save(AC_meth_values_JM, file = paste0(save_dir, "AC_meth_values_JM.Rdata"))

##AS
write.csv(AS_meth_values_JM, file = paste0(save_dir, "AS_meth_values_JM.csv") )
save(AS_meth_values_JM, file = paste0(save_dir, "AS_meth_values_JM.Rdata"))

##EH
write.csv(EH_meth_values_JM, file = paste0(save_dir, "EH_meth_values_JM.csv") )
save(EH_meth_values_JM, file = paste0(save_dir, "EH_meth_values_JM.Rdata"))

##JM
write.csv(JM_meth_values_JM, file = paste0(save_dir, "JM_meth_values_JM.csv") )
save(JM_meth_values_JM, file = paste0(save_dir, "JM_meth_values_JM.Rdata"))

##ZF
write.csv(ZF_meth_values_JM, file = paste0(save_dir, "ZF_meth_values_JM.csv") )
save(ZF_meth_values_JM, file = paste0(save_dir, "ZF_meth_values_JM.Rdata"))

#### PCA tests and stat tests ####

##AC
PCA_AC <- prcomp(AC_meth_values,scale = TRUE)
PCA_values_AC <- as.data.frame(PCA_AC$x)
PCA_values_AC$species <- "AC"

##ZF
PCA_AS <- prcomp(AS_meth_values,scale = TRUE)
PCA_values_AS <- as.data.frame(PCA_AS$x)
PCA_values_AS$species <- "AS"

##EH
PCA_EH <- prcomp(EH_meth_values,scale = TRUE)
PCA_values_EH <- as.data.frame(PCA_EH$x)
PCA_values_EH$species <- "EH"

##JM
PCA_JM <- prcomp(JM_meth_values,scale = TRUE)
PCA_values_JM <- as.data.frame(PCA_JM$x)
PCA_values_JM$species <- "JM"

##ZF
PCA_ZF <- prcomp(ZF_meth_values,scale = TRUE)
PCA_values_ZF <- as.data.frame(PCA_ZF$x)
PCA_values_ZF$species <- "ZF"


# PCA_values_all <- rbind(PCA_values_AC[,c(1:10, length(PCA_values_AC))],
#                         PCA_values_AS[,c(1:10, length(PCA_values_AS))],
#                         PCA_values_EH[,c(1:10, length(PCA_values_EH))])
# age_all <- c(AC_age/mean(AC_age), AS_age/mean(AS_age), EH_age/mean(EH_age))

#### plotting results ####

# AC_age_tmp <- AC_age[c(-71,-74,-27,-69,-2,-3,-17)]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##AC
ggplot(PCA_values_AC, aes(x = PC1, y = PC2, color = AC_age, shape = species)) +
  geom_point(aes(cex = 7)) +
  facet_wrap(~species) +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6]) +
  geom_text(aes(label = 1:nrow(AC_meth_data)), nudge_x = 0.2, nudge_y = 0) +
  theme_classic()

##AS
ggplot(PCA_values_AS, aes(x = PC1, y = PC2, color = AS_age, shape = species)) +
  geom_point(aes(cex = 7)) +
  facet_wrap(~species) +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6]) +
  geom_text(aes(label = AS_meth_data$id), nudge_x = 0.2, nudge_y = 0) +
  theme_classic()

##EH
ggplot(PCA_values_EH, aes(x = PC2, y = PC1, color = EH_age)) +
  geom_point(aes(cex = 7)) +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6]) +
  geom_text(aes(label = EH_meth_data$id), nudge_x = 0.2, nudge_y = 0) +
  theme_classic()

##ZF
ggplot(PCA_values_ZF, aes(x = PC1, y = PC2, color = ZF_age)) +
  geom_point(aes(cex = 7)) +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6]) +
  geom_text(aes(label = rownames(ZF_meth_data)), nudge_x = 0.2, nudge_y = 0) +
  theme_classic()

#### Transforming for plots ####

#run this if JM is rgenome 
AC_meth_values <- AC_meth_values_JM
AS_meth_values <- AS_meth_values_JM
EH_meth_values <- EH_meth_values_JM
JM_meth_values <- JM_meth_values_JM
ZF_meth_values <- ZF_meth_values_JM

AC_methyl_sites <- AC_methyl_sites_JM
AS_methyl_sites <- AS_methyl_sites_JM
EH_methyl_sites <- EH_methyl_sites_JM
JM_methyl_sites <- JM_methyl_sites_JM
ZF_methyl_sites <- ZF_methyl_sites_JM

AS_meth_values_long <- pivot_longer(AS_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AS_meth_values_long$age <- rep(AS_age, each = ncol(AS_meth_values))
AS_meth_values_long$SMR <- as.factor(rep(AS_methyl_sites$SMR, times = length(AS_age)))
AS_meth_values_long$Site_i <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.integer()
AS_meth_values_long$Site_f <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.factor()
AS_meth_values_long$species <- "AS"

AC_meth_values_long <- pivot_longer(AC_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AC_meth_values_long$age <- rep(AC_age, each = ncol(AC_meth_values))
AC_meth_values_long$SMR <- as.factor(rep(AC_methyl_sites$SMR[meth_sites_names_tmp_AC %in% AC_meth_data_test], times = length(AC_age))) # indexing is necessary because not all CpGs were able to be extracted from the shared sites due o batch correction
AC_meth_values_long$Site_i <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.integer()
AC_meth_values_long$Site_f <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.factor()
AC_meth_values_long$species <- "AC"

EH_meth_values_long <- pivot_longer(EH_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
EH_meth_values_long$age <- rep(EH_meth_data$age, each = ncol(EH_meth_values))
EH_meth_values_long$SMR <- as.factor(rep(EH_methyl_sites$SMR, times = length(EH_age)))
EH_meth_values_long$Site_i <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.integer()
EH_meth_values_long$Site_f <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.factor()
EH_meth_values_long$species <- "EH"
# EH_meth_values_long$id <- rep(paste0("EH_", EH_meth_data$id), times = ncol(EH_meth_values))

JM_meth_values_long <- pivot_longer(JM_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
JM_meth_values_long$age <- rep(JM_meth_data$age, each = ncol(JM_meth_values))
JM_meth_values_long$SMR <- as.factor(rep(JM_methyl_sites$SMR, times = length(JM_age)))
JM_meth_values_long$Site_i <- gsub(".*\\:", "", JM_meth_values_long$Site) %>% as.integer()
JM_meth_values_long$Site_f <- gsub(".*\\:", "", JM_meth_values_long$Site) %>% as.factor()
JM_meth_values_long$species <- "JM"
# C_meth_values_long$sample <- rep(paste0("sample_", 1:length(AC_meth_values)),each = ncol(AC_age))

ZF_meth_values_long <- pivot_longer(ZF_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
ZF_meth_values_long$age <- rep(ZF_meth_data$age, each = ncol(ZF_meth_values))
ZF_meth_values_long$SMR <- as.factor(rep(ZF_methyl_sites$SMR, times = length(ZF_age)))
ZF_meth_values_long$Site_i <- gsub(".*\\:", "", ZF_meth_values_long$Site) %>% as.integer()
ZF_meth_values_long$Site_f <- gsub(".*\\:", "", ZF_meth_values_long$Site) %>% as.factor()
ZF_meth_values_long$species <- "ZF"

all_meth_values_long <- rbind(AC_meth_values_long, AS_meth_values_long, EH_meth_values_long, JM_meth_values_long, ZF_meth_values_long)


##color palettes
colpalOI <- palette.colors(palette = "Okabe-Ito") %>% 
  as.vector()
colpal <- hcl.colors(7, "SunsetDark") 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


ggplot(AS_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = colpal[1], high = colpal[7], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values AS")

ggplot(AC_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values AC")

ggplot(EH_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values EH")

ggplot(JM_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, cex = 0.1,  alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = "#E69F00", high = "#0072B2", guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values JM")

## all boxplot and sina
ggplot(all_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina() +
  geom_boxplot(aes(group = Site_f, fill = species), alpha = 0.9) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_fill_manual(values = colpalOI[c(-1,-9)]) +
  # scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Cod (AC), Snapper (ZF), Hake (EH), Medaka (JM), Zebrafish (ZF)")

## all boxplot
ggplot(all_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  # geom_sina(aes(color = species)) +
  geom_boxplot(aes(group = Site_f, fill = species), alpha = 0.9) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_fill_manual(values = colpalOI[c(-1,-9)]) +
  # scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Cod (AC), Snapper (ZF), Hake (EH), Medaka (JM), Zebrafish (ZF)")

## all sina
ggplot(all_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = species)) +
  # geom_boxplot(aes(group = Site_f, fill = species), alpha = 0.9) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_fill_manual(values = colpalOI[c(-1,-9)]) +
  # scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Cod (AC), Snapper (ZF), Hake (EH), Medaka (JM), Zebrafish (ZF)")



#### testing correlation tests ####
# Initialize a list to store the results
library(tibble)
correlation_results_ZF <- list()

#### function ####

cor.test.age <- function(methyl_values, age, method = "pearson") {
  correlation_results <- list()
  print(paste0("Running correlation test against age with ", method, " method. Results are stored in tibble."))
  # Loop through each methylation site
  for (i in 1:ncol(methyl_values)) {
    site_name <- colnames(methyl_values)[i]
    # Perform correlation test with age
    test_result <- cor.test(methyl_values[,i], age, method = method) # Use "spearman" or "kendall" if more appropriate
    
    # Store the results
    correlation_results[[site_name]] <- list(
      correlation_coefficient = test_result$estimate,
      p_value = test_result$p.value
    )
  }
  
  # Optionally, convert the results list to a more convenient format like a dataframe
  correlation_summary <- tibble(
    Site = names(correlation_results),
    Correlation = sapply(correlation_results, function(x) x$correlation_coefficient),
    P_value = sapply(correlation_results, function(x) x$p_value)
  )
  return(correlation_summary)
}

cor.test.age.filter <- function(input, p_value = 0.05) {
  input[input$P_value <= p_value,]
}


AC_cor_age_pearson <- cor.test.age(AC_meth_values, AC_meth_data$age)
AC_cor_age_filtered_pearson <- cor.test.age.filter(AC_cor_age_pearson, 0.05)

AC_cor_age_spearman <- cor.test.age(AC_meth_values, AC_meth_data$age, method = "spearman")
AC_cor_age_filtered_spearman <- cor.test.age.filter(AC_cor_age_spearman, 0.05)

AC_cor_age_kendall <- cor.test.age(AC_meth_values, AC_meth_data$age, method = "kendall")
AC_cor_age_filtered_kendall <- cor.test.age.filter(AC_cor_age_kendall, 0.05)

AS_cor_age_pearson <- cor.test.age(AS_meth_values, AS_meth_data$age)
AS_cor_age_filtered_pearson <- cor.test.age.filter(AS_cor_age, 0.05)

EH_cor_age_pearson <- cor.test.age(EH_meth_values, EH_meth_data$age)
EH_cor_age_filtered_pearson <- cor.test.age.filter(EH_cor_age, 0.05)

JM_cor_age_pearson <- cor.test.age(JM_meth_values, JM_meth_data$age)
JM_cor_age_filtered_pearson <- cor.test.age.filter(JM_cor_age, 0.05)

ZF_cor_age_pearson <- cor.test.age(ZF_meth_values, ZF_meth_data$age)
ZF_cor_age_filtered_pearson <- cor.test.age.filter(ZF_cor_age, 0.05)

