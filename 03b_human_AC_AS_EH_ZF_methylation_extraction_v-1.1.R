#### Overview ####
# Extraction of methylation data based on overlapping sequences

#### Settings ####
setwd("/powerplant/workspace/cfngle")

#### Preparation ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(tidyverse)

## loading data 

load("results-data/02_conserved_seq/HS_AC_AS_EH_ZF_overlaps_bt2.Rdata")


HS_overlap_seqs_bt2


overlap_HS_AC_bt2 <- HS_overlap_seqs_bt2[[1]]
overlap_HS_AS_bt2 <- HS_overlap_seqs_bt2[[2]]
overlap_HS_EH_bt2 <- HS_overlap_seqs_bt2[[3]]
overlap_HS_ZF_bt2 <- HS_overlap_seqs_bt2[[4]]

## getting SMR 
HS_gr_overlap_seqs_bt2 <- lapply(HS_overlap_seqs_bt2, function(x) granges(x))

HS_group_gr_overlap_bt2 <-  c(HS_gr_overlap_seqs_bt2[[1]], HS_gr_overlap_seqs_bt2[[2]], HS_gr_overlap_seqs_bt2[[3]], HS_gr_overlap_seqs_bt2[[4]])
HS_SMR_b_bt2 <- GenomicRanges::reduce(HS_group_gr_overlap_bt2)
names(HS_SMR_b_bt2) <- sprintf("HS_SMR_b_bt2_%03d", 1:length(HS_SMR_b_bt2))

HS_AC_AS_EH_ZF_SMR_b_bt2 <- HS_SMR_b_bt2

save(HS_AC_AS_EH_ZF_SMR_b_bt2, file = "/workspace/cfngle/results-data/04_SMRs/HS_AC_AS_EH_ZF_SMR_b_bt2.Rdata")

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
ZF_methyl_df_bt2 <- get.methyl.sites(overlap_HS_ZF_bt2, species = "JM", SMRs = HS_SMR_b_bt2)

combined_df_bt2 <- bind_rows(AC_methyl_df_bt2, AS_methyl_df_bt2, EH_methyl_df_bt2, ZF_methyl_df_bt2)

HS_chr_names <- sort(unique(combined_df_bt2$chr_align))
# AC_chr_names_simple <-  paste0("Chr ", 1:23)

# combined_df_bt2 <- cbind(combined_df_bt2, data.frame(align_method = "Bowtie2")) 

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

#### get methyl sites and value ####
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
ZF_methyl_sites <- get.methyl.sites(overlap_HS_ZF_bt2, species = "ZF", SMRs = HS_SMR_b_bt2)

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
# tail(colnames(AS_meth_data))
AS_age <- AS_meth_data$age

##EH
xx <- load("/workspace/cfngle/raw-data/EH/zzz-methyl_data/Meth-complete-hake.RData")
assign("EH_meth_data", get(xx))
# tail(colnames(EH_meth_data))
EH_age <- EH_meth_data$age
EH_metadata_samples <- read.csv("/workspace/cfngle/raw-data/EH/zzz-methyl_data/hake-samples.txt", sep = "\t")
EH_sex <- EH_metadata_samples$sex


##ZF
ZF_meth_data <- load("/workspace/cfngle/raw-data/ZF/zzz_methyldata/ZF_methyldata_88.RData")
ZF_meth_data <- ZF_methyl_data

ZF_age <- ZF_meth_data$age/52

#### extract methylation data for all samples ####
# Not all the datasets have the same naming structure, hence the steps are different and are done one by one

##AC
meth_sites_names_tmp_AC <- paste0(AC_methyl_sites$Chr, ".", AC_methyl_sites$pos_rgenome) %>% 
  gsub("AC_", "Chr", .)
AC_meth_data_test <- gsub("X", "Chr", colnames(AC_meth_data))
table(meth_sites_names_tmp_AC %in% AC_meth_data_test)
meth_columns_tmp <- sapply(meth_sites_names_tmp_AC[meth_sites_names_tmp_AC %in% AC_meth_data_test], function(x) grep(x, AC_meth_data_test)) %>% 
  as.vector()

AC_meth_values <- AC_meth_data[ ,meth_columns_tmp]

##AS
meth_sites_names_tmp <- paste0(AS_methyl_sites$Chr, "-", AS_methyl_sites$pos_rgenome)
meth_sites_names_tmp %in% colnames(AS_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(AS_meth_data))) %>% 
  as.vector()

AS_meth_values <- AS_meth_data[,meth_columns_tmp] 

##EH
meth_sites_names_tmp <- paste0(gsub("EH_", "",EH_methyl_sites$Chr), ".", EH_methyl_sites$pos_rgenome)
meth_sites_names_tmp %in% colnames(EH_meth_data)
meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, colnames(EH_meth_data))) %>% 
  as.vector()

EH_meth_values <- EH_meth_data[,meth_columns_tmp] 
# if males ONLY
EH_meth_values <- EH_meth_values[EH_sex == "M",]
EH_age <- EH_age[EH_sex == "M"]

##ZF
meth_sites_names_tmp <- paste0(gsub("ZF_", "",ZF_methyl_sites$Chr), ":", ZF_methyl_sites$pos_rgenome)
table(meth_sites_names_tmp %in% colnames(ZF_meth_data))
meth_columns_tmp <- ZF_meth_data[meth_sites_names_tmp] %>% as.vector()
meth_columns_tmp <- unlist(meth_columns_tmp)
# ZF_meth_values_JM <- ZF_meth_data[,meth_columns_tmp] 
ZF_meth_values <- ZF_meth_data[meth_sites_names_tmp]

### saving data
save_dir <- "/workspace/cfngle/results-data/05_shared_methyl_values/"

##AC
write.csv(AC_meth_values, file = paste0(save_dir, "HS_AC_meth_values.csv") )
save(AC_meth_values, file = paste0(save_dir, "HS_AC_meth_values.Rdata"))

##AS
write.csv(AS_meth_values, file = paste0(save_dir, "HS_AS_meth_values.csv") )
save(AS_meth_values, file = paste0(save_dir, "HS_AS_meth_values.Rdata"))

##EH
write.csv(EH_meth_values, file = paste0(save_dir, "HS_EH_meth_values.csv") )
save(EH_meth_values, file = paste0(save_dir, "HS_EH_meth_values.Rdata"))

##ZF
write.csv(ZF_meth_values, file = paste0(save_dir, "HS_ZF_meth_values.csv") )
save(ZF_meth_values, file = paste0(save_dir, "HS_ZF_meth_values.Rdata"))

#### Imputation ####

# There are several ways of imputing missing values, here we present two of them. Always set.seed() for imputations.
#a) Method 1 using package “mice” (Multiple Imputation by Chained Equation)
library(mice)

set.seed(123)

init <- mice(ZF_meth_values, maxit=0)

m_method <- init$method

pred_matrix <- init$predictorMatrix

colnames(ZF_meth_values)
predM[, c("age")]=0
meth[c("age")]=""

imputed <- mice(ZF_meth_values, method = m_method, predictorMatrix = pred_matrix, m = 5)

#b) Method 2 using package “zoo” (Missing values replaced by the mean or other function of its group)
library(zoo)
set.seed(123)
ZF_meth_values_imputed <- na.aggregate(ZF_meth_values)

#### PCA ####

##AC
PCA_AC <- prcomp(AC_meth_values,scale = TRUE)
PCA_values_AC <- as.data.frame(PCA_AC$x)
PCA_values_AC$species <- "AC"

##AS
PCA_AS <- prcomp(AS_meth_values,scale = TRUE)
PCA_values_AS <- as.data.frame(PCA_AS$x)
PCA_values_AS$species <- "AS"

##EH
PCA_EH <- prcomp(EH_meth_values,scale = TRUE)
PCA_values_EH <- as.data.frame(PCA_EH$x)
PCA_values_EH$species <- "EH"

##ZF
PCA_ZF <- prcomp(ZF_meth_values_imputed,scale = TRUE)
PCA_values_ZF <- as.data.frame(PCA_ZF$x)
PCA_values_ZF$species <- "ZF"

### plotting PCA
library(patchwork)
# color palettes
colpalOI <- palette.colors(palette = "Okabe-Ito") %>% 
  as.vector() %>% 
  .[c(-1,-9)]
colpal <- hcl.colors(7, "SunsetDark") 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##AC
AC_pca_plot <- ggplot(PCA_values_AC, aes(x = PC1, y = PC2, color = AC_age)) +
  geom_point(cex = 3) +
  facet_wrap(~species) +
  scale_color_gradient(low = colpalOI[1], high = colpalOI[2]) +
  # geom_text(aes(label = 1:nrow(AC_meth_data)), nudge_x = 0.4, nudge_y = 0) +
  theme_classic()

##AS
AS_pca_plot <- ggplot(PCA_values_AS, aes(x = PC1, y = PC2, color = AS_age)) +
  geom_point(cex = 3) +
  facet_wrap(~species) +
  scale_color_gradient(low = colpalOI[3], high = colpalOI[4]) +
  # geom_text(aes(label = AS_meth_data$id), nudge_x = 0, nudge_y = 0.5) +
  theme_classic()

##EH
EH_pca_plot <- ggplot(PCA_values_EH, aes(x = PC2, y = PC1, color = EH_age)) +
  geom_point(cex = 3) +
  facet_wrap(~species) +
  scale_color_gradient(low = colpalOI[5], high = colpalOI[6]) +
  # geom_text(aes(label = EH_meth_data$id), nudge_x = 0.4, nudge_y = 0) +
  theme_classic()

##ZF
ZF_pca_plot <- ggplot(PCA_values_ZF, aes(x = PC1, y = PC2, color = ZF_age)) +
  geom_point(cex = 3) +
  facet_wrap(~species) +
  scale_color_gradient(low = colpalOI[1], high = colpalOI[5]) +
  # geom_text(aes(label = rownames(ZF_meth_data)), nudge_x = 0.4, nudge_y = 0) +
  theme_classic()

## plot all 
AC_pca_plot + AS_pca_plot + EH_pca_plot + ZF_pca_plot +
  plot_layout(nrow=2)


#### Methylation values ####

AS_meth_values_long <- pivot_longer(AS_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AS_meth_values_long$age <- rep(AS_age, each = ncol(AS_meth_values))
AS_meth_values_long$max_age <- 54
AS_meth_values_long$rel_age <- AS_meth_values_long$age / AS_meth_values_long$max_age
AS_meth_values_long$SMR <- as.factor(rep(AS_methyl_sites$SMR, times = length(AS_age)))
AS_meth_values_long$Site_i <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.integer()
AS_meth_values_long$Site_f <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.factor()
AS_meth_values_long$species <- "AS"

AC_meth_values_long <- pivot_longer(AC_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AC_meth_values_long$age <- rep(AC_age, each = ncol(AC_meth_values))
AC_meth_values_long$max_age <- 25
AC_meth_values_long$rel_age <- AC_meth_values_long$age / AC_meth_values_long$max_age
AC_meth_values_long$SMR <- as.factor(rep(AC_methyl_sites$SMR[meth_sites_names_tmp_AC %in% AC_meth_data_test], times = length(AC_age))) # indexing is necessary because not all CpGs were able to be extracted from the shared sites due o batch correction
AC_meth_values_long$Site_i <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.integer()
AC_meth_values_long$Site_f <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.factor()
AC_meth_values_long$species <- "AC"

EH_meth_values_long <- pivot_longer(EH_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
EH_meth_values_long$age <- rep(EH_meth_data$age, each = ncol(EH_meth_values))
EH_meth_values_long$max_age <- 20
EH_meth_values_long$rel_age <- EH_meth_values_long$age / EH_meth_values_long$max_age
EH_meth_values_long$SMR <- as.factor(rep(EH_methyl_sites$SMR, times = length(EH_age)))
EH_meth_values_long$Site_i <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.integer()
EH_meth_values_long$Site_f <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.factor()
EH_meth_values_long$species <- "EH"

ZF_meth_values_long <- pivot_longer(ZF_meth_values_imputed, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
ZF_meth_values_long$age <- rep(ZF_meth_data$age, each = ncol(ZF_meth_values)) / 52
ZF_meth_values_long$max_age <- 5
ZF_meth_values_long$rel_age <- ZF_meth_values_long$age / ZF_meth_values_long$max_age
ZF_meth_values_long$SMR <- as.factor(rep(ZF_methyl_sites$SMR, times = length(ZF_age)))
ZF_meth_values_long$Site_i <- gsub(".*\\:", "", ZF_meth_values_long$Site) %>% as.integer()
ZF_meth_values_long$Site_f <- gsub(".*\\:", "", ZF_meth_values_long$Site) %>% as.factor()
ZF_meth_values_long$species <- "ZF"

all_meth_values_long <- rbind(AC_meth_values_long, AS_meth_values_long, EH_meth_values_long, ZF_meth_values_long)

### plotting

ggplot(AS_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = colpal[1], high = colpal[7], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values AS (human rgenome)")

ggplot(AC_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values AC (human rgenome)")

ggplot(EH_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values EH (human rgenome)")

ggplot(ZF_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, cex = 0.1,  alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = "#E69F00", high = "#0072B2", guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values ZF (human rgenome)")

## all boxplot
ggplot(all_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_boxplot(aes(group = Site_f, fill = species), alpha = 0.9) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_fill_manual(values = colpalOI) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (ZF), European Hake (EH), Zebrafish (ZF) (human rgenome)")


#### correlation testing ####
## correlation test between the selected CpGs and to age as well

#### Preparation ####
library(tibble)
library(dplyr)
library(ggplot2)



#### functions ####

cor.test.age <- function(methyl_values, age, SMR = "not_defined", species = "undefined", method = "pearson") {
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
    P_value = sapply(correlation_results, function(x) x$p_value),
    SMR = SMR,
    species = species
  )
  return(correlation_summary)
}

cor.test.age.filter <- function(input, p_value = 0.05) {
  significant_vector <- as.vector(ifelse(input$P_value <= p_value, TRUE, FALSE))
  input$significant <- significant_vector
  return(input)
}


#### Correlation tests ####

AC_cor_age_pearson <- cor.test.age(AC_meth_values, AC_age, AC_methyl_sites$SMR[meth_sites_names_tmp_AC %in% AC_meth_data_test], species = "AC")
AC_cor_age_filtered_pearson <- cor.test.age.filter(AC_cor_age_pearson, 0.05)

AS_cor_age_pearson <- cor.test.age(AS_meth_values, AS_age, AS_methyl_sites$SMR, species = "AS")
AS_cor_age_filtered_pearson <- cor.test.age.filter(AS_cor_age_pearson, 0.05)

EH_cor_age_pearson <- cor.test.age(EH_meth_values, EH_age, EH_methyl_sites$SMR, species = "EH")
EH_cor_age_filtered_pearson <- cor.test.age.filter(EH_cor_age_pearson, 0.05)

ZF_cor_age_pearson <- cor.test.age(ZF_meth_values, ZF_age, ZF_methyl_sites$SMR, species = "ZF")
ZF_cor_age_filtered_pearson <- cor.test.age.filter(ZF_cor_age_pearson, 0.05)

cor_all <- rbind(AC_cor_age_filtered_pearson,
                 AS_cor_age_filtered_pearson,
                 EH_cor_age_filtered_pearson, 
                 ZF_cor_age_filtered_pearson)

## selecting CpGs

ncol(AC_meth_values[AC_cor_age_filtered_pearson$significant] == TRUE)
ncol(AS_meth_values[AS_cor_age_filtered_pearson$significant] == TRUE)
ncol(EH_meth_values[EH_cor_age_filtered_pearson$significant] == TRUE)
ncol(ZF_meth_values[ZF_cor_age_filtered_pearson$significant] == TRUE)

#### function to choose only the highest correlating CpGs per SMR

select.max.cor <- function(cor_tibble, filter_significant = FALSE) {
  filtered_data <- cor_tibble
  if(filter_significant == TRUE) {filtered_data <- filter(filtered_data, significant)}
   
  filtered_data <- filtered_data %>% 
    group_by(SMR) %>%
    # Add a temporary column for the absolute correlation values
    mutate(abs_correlation = abs(Correlation)) %>%
    # For each group, filter the row with the max absolute correlation
    filter(abs_correlation == max(abs_correlation)) %>%
    # Remove the temporary column
    select(-abs_correlation) %>%
    # Optionally, ensure only one row per group if there are ties
    slice(1)
  return(filtered_data)
}

select.max.cor(AC_cor_age_filtered_pearson[AC_cor_age_filtered_pearson$Correlation > 0,])
select.max.cor(AS_cor_age_filtered_pearson)
select.max.cor(EH_cor_age_filtered_pearson)
select.max.cor(ZF_cor_age_filtered_pearson)

AC_pos_cor_CpGs <- select.max.cor(AC_cor_age_filtered_pearson[AC_cor_age_filtered_pearson$Correlation > 0,])
AS_pos_cor_CpGs <- select.max.cor(AS_cor_age_filtered_pearson[AS_cor_age_filtered_pearson$Correlation > 0,])
EH_pos_cor_CpGs <- select.max.cor(EH_cor_age_filtered_pearson[EH_cor_age_filtered_pearson$Correlation > 0,])
ZF_pos_cor_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson[ZF_cor_age_filtered_pearson$Correlation > 0,])

all_pos_cor_CpG <- rbind(AC_pos_cor_CpGs, AS_pos_cor_CpGs, EH_pos_cor_CpGs, ZF_pos_cor_CpGs)
all_pos_cor_CpG_common  <- all_pos_cor_CpG %>% 
  group_by(SMR) %>% 
  filter(n() == 4) %>% 
  ungroup

AC_sig_CpGs <- select.max.cor(AC_cor_age_filtered_pearson, TRUE)
AS_sig_CpGs <- select.max.cor(AS_cor_age_filtered_pearson, TRUE)
EH_sig_CpGs <- select.max.cor(EH_cor_age_filtered_pearson, TRUE)
ZF_sig_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson, TRUE)

all_sig_CpGs <- rbind(AC_sig_CpGs, AS_sig_CpGs, EH_sig_CpGs, ZF_sig_CpGs)

all_sig_CpGs_common <- all_sig_CpGs %>% 
  group_by(SMR) %>% 
  filter(n() == 4) %>% 
  ungroup


#### model creation ####
library(glmnet)
AC_glm <- cv.glmnet(as.matrix(AC_meth_values), AC_age, alpha = 0.5)
coef(AC_glm, s=0.01)
AC_glm_prediction <- predict(AC_glm, newx = as.matrix(AC_meth_values))
plot(AC_glm)


### selecting only pos cor ones 
#AC
AC_meth_values_selected <- AC_meth_values[,colnames(AC_meth_values) %in% all_pos_cor_CpG_common$Site]
AC_name_index <- match(colnames(AC_meth_values_selected), all_pos_cor_CpG_common$Site)
colnames(AC_meth_values_selected) <- all_pos_cor_CpG_common$SMR[AC_name_index]
AC_meth_values_selected <- AC_meth_values_selected[, order(colnames(AC_meth_values_selected))]
AC_meth_values_selected$rel_age <- AC_age/25
AC_meth_values_selected$species <- "AC"

#AS
AS_meth_values_selected <- AS_meth_values[,colnames(AS_meth_values) %in% all_pos_cor_CpG_common$Site]
AS_name_index <- match(colnames(AS_meth_values_selected), all_pos_cor_CpG_common$Site)
colnames(AS_meth_values_selected) <- all_pos_cor_CpG_common$SMR[AS_name_index]
AS_meth_values_selected <- AS_meth_values_selected[, order(colnames(AS_meth_values_selected))]
AS_meth_values_selected$rel_age <- AS_age/54
AS_meth_values_selected$species <- "AS"

#EH
EH_meth_values_selected <- EH_meth_values[,colnames(EH_meth_values) %in% all_pos_cor_CpG_common$Site]
EH_name_index <- match(colnames(EH_meth_values_selected), all_pos_cor_CpG_common$Site)
colnames(EH_meth_values_selected) <- all_pos_cor_CpG_common$SMR[EH_name_index]
EH_meth_values_selected <- EH_meth_values_selected[, order(colnames(EH_meth_values_selected))]
EH_meth_values_selected$rel_age <- EH_age/20
EH_meth_values_selected$species <- "EH"

#ZF
ZF_meth_values_selected <- ZF_meth_values_imputed[,colnames(ZF_meth_values) %in% all_pos_cor_CpG_common$Site]
ZF_name_index <- match(colnames(ZF_meth_values_selected), all_pos_cor_CpG_common$Site)
colnames(ZF_meth_values_selected) <- all_pos_cor_CpG_common$SMR[ZF_name_index]
ZF_meth_values_selected <- ZF_meth_values_selected[, order(colnames(ZF_meth_values_selected))]
ZF_meth_values_selected$rel_age <- ZF_age/5
ZF_meth_values_selected$species <- "ZF"


all_meth_values_selected <- rbind(AC_meth_values_selected, AS_meth_values_selected, EH_meth_values_selected, ZF_meth_values_selected)

all_meth_values_selected_train <- all_meth_values_selected[-seq(1, nrow(all_meth_values_selected), 4),]
all_meth_values_selected_test <- all_meth_values_selected[seq(1, nrow(all_meth_values_selected), 4),]

X <- all_meth_values_selected_train %>% 
  select(-rel_age, -species) %>% 
  as.matrix()

Y <- all_meth_values_selected_train[,"rel_age"]

X_test <- all_meth_values_selected_test %>% 
  select(-rel_age, -species) %>% 
  as.matrix()
Y_test <- all_meth_values_selected_test[,"rel_age"]

set.seed(123)
glm_test <- cv.glmnet(X, Y, alpha = 0.5)
plot(glm_test)

coef(glm_test, s=glm_test$lambda.min)

predictions_test <- predict(glm_test, newx= X_test, s = glm_test$lambda.min)
result_df <- data.frame(predictions = predictions_test,
                        rel_age = all_meth_values_selected_test$rel_age,
                        species = all_meth_values_selected_test$species)
colnames(result_df) <- c("age_predicted", "age", "species")

ggplot(result_df, aes(color = species)) +
  geom_point(aes(x = age, y = age_predicted)) +
  scale_color_manual(values = colpalOI) +
  ylim(0,0.65) +
  xlim(0,0.65) +
  theme_minimal()


## evaluation

mse <- mean((predictions_test - all_meth_values_selected_test$rel_age)^2)
print(paste("MSE:", mse))

rmse <- sqrt(mse)
print(paste("RMSE:", rmse))

mae <- mean(abs(predictions_test - all_meth_values_selected_test$rel_age))
print(paste("MAE:", mae))

pearson_cor <- cor(predictions_test, all_meth_values_selected_test$rel_age, method = "pearson")


#### plotting ####
colpalOI <- palette.colors(palette = "Okabe-Ito") %>% 
  as.vector() %>%
  .[c(-1,-9)]
colpal <- hcl.colors(7, "SunsetDark") 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## all
ggplot(cor_all, aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species, alpha = significant)) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  # facet_row(~SMR) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggplot(cor_all, aes()) +
  geom_point(aes(x = Site, y = Correlation, color = species, alpha = significant)) +
  # geom_line(aes(x = c(-1,1), y = log2(0.05), color = "#CC79A7")) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  facet_wrap(~SMR) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## only significant
ggplot(subset(cor_all, significant == TRUE), aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species)) +
  # facet_wrap(~SMR) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## only selected 
ggplot(all_sig_CpGs_common, aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species)) +
  facet_row(~SMR) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## only cor positive 
ggplot(all_pos_cor_CpG, aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species, alpha = significant)) +
  facet_row(~SMR) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  theme(axis.ticks.x = element_blank())

## plotting SMR groups 24 and 28
selected_methyl_values <- subset(all_meth_values_long, Site %in% subset(all_sig_CpGs_common, SMR == "SMR_024" | SMR == "SMR_026")$Site)
selected_methyl_values <- subset(all_meth_values_long, Site %in% all_pos_cor_CpG$Site)

ggplot(selected_methyl_values, aes(x = Site, y = Methylation_Value)) +
  geom_boxplot(aes(group = Site_f, fill = species), alpha = 0.9) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_fill_manual(values = colpalOI) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (ZF), European Hake (EH), Zebrafish (ZF) (human rgenome)")

ggplot(selected_methyl_values, aes(x = species, y = Methylation_Value)) +
  geom_sina(aes(color = rel_age, shape = species)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_manual(aesthetics = "legend") +
  theme_classic() +
  # theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (ZF), European Hake (EH), Zebrafish (ZF) (human rgenome)")
