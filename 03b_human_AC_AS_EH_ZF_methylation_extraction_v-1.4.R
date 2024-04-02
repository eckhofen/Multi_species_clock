#### Overview ####
# Extraction of methylation data based on overlapping sequences

#### Settings ####
setwd("/powerplant/workspace/cfngle")

#### Preparation ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(tidyverse)
library(ggforce)

#### getting SHARED METHYLATION REGION

## loading data 
# loading overlapping alignment reads (rgenome: human)
load("results-data/02_conserved_seq/HS_AC_AS_EH_ZF_overlaps_bt2.Rdata")

# assign overlapping sequences for each species
overlap_HS_AC_bt2 <- HS_overlap_seqs_bt2[[1]]
overlap_HS_AS_bt2 <- HS_overlap_seqs_bt2[[2]]
overlap_HS_EH_bt2 <- HS_overlap_seqs_bt2[[3]]
overlap_HS_ZF_bt2 <- HS_overlap_seqs_bt2[[4]]

# transforming aligned reads into GRanges object
HS_gr_overlap_seqs_bt2 <- lapply(HS_overlap_seqs_bt2, function(x) granges(x))
# using the overlap of the sequences to get the SMRs
HS_group_gr_overlap_bt2 <-  c(HS_gr_overlap_seqs_bt2[[1]], HS_gr_overlap_seqs_bt2[[2]], HS_gr_overlap_seqs_bt2[[3]], HS_gr_overlap_seqs_bt2[[4]])
HS_SMR_b_bt2 <- GenomicRanges::reduce(HS_group_gr_overlap_bt2)
names(HS_SMR_b_bt2) <- sprintf("HS_SMR_b_bt2_%03d", 1:length(HS_SMR_b_bt2))
# renaming
HS_AC_AS_EH_ZF_SMR_b_bt2 <- HS_SMR_b_bt2

save(HS_AC_AS_EH_ZF_SMR_b_bt2, file = "/workspace/cfngle/results-data/04_SMRs/HS_AC_AS_EH_ZF_SMR_b_bt2.Rdata")

### function for methylation extraction

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

AC_methyl_df_bt2 <- get.methyl.sites(overlap_HS_AC_bt2, species = "AC", SMRs = HS_SMR_b_bt2)
AS_methyl_df_bt2 <- get.methyl.sites(overlap_HS_AS_bt2, species = "AS", SMRs = HS_SMR_b_bt2)
EH_methyl_df_bt2 <- get.methyl.sites(overlap_HS_EH_bt2, species = "EH", SMRs = HS_SMR_b_bt2)
ZF_methyl_df_bt2 <- get.methyl.sites(overlap_HS_ZF_bt2, species = "ZF", SMRs = HS_SMR_b_bt2)

combined_df_bt2 <- bind_rows(AC_methyl_df_bt2, AS_methyl_df_bt2, EH_methyl_df_bt2, ZF_methyl_df_bt2)

HS_chr_names <- sort(unique(combined_df_bt2$chr_align))

### plotting data 

library(ggplot2)
# library(ggpattern)

# adding color palette
colpalOI <- palette.colors(palette = "Okabe-Ito") %>% 
  as.vector() %>%
  .[c(-1,-9)]

# this is visualising the genomic locations and the frequency of CpGs per species
ggplot(combined_df_bt2, aes(x = pos_align, fill = species)) +
  geom_histogram(position = "dodge") +  
  facet_wrap(~ chr_align, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "All CpGs") +
  scale_color_manual(values = colpalOI) + 
  labs(title = "Number of CpGs and their genomic position") +
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0), )  

#### getting methyl sites for all species

AC_methyl_sites <- AC_methyl_df_bt2
AS_methyl_sites <- AS_methyl_df_bt2
EH_methyl_sites <- EH_methyl_df_bt2
ZF_methyl_sites <- ZF_methyl_df_bt2

### load methylation data for samples as well as age vector for each species

##AC
xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/meth-corrected-batchcorrected-cod.Rdata")
# xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/Meth-complete-nobatchcorrection-cod.RData")
assign("AC_meth_data", get(xx))
AC_meth_data <- as.data.frame(AC_meth_data)
# tail(colnames(AC_meth_data))
AC_age <- AC_meth_data$age
AC_age[AC_age == 0] <- 0.01

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
EH_age <- EH_metadata_samples$age

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
AC_methyl_sites <- AC_methyl_sites[meth_sites_names_tmp_AC %in% AC_meth_data_test,]

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

# EH_meth_values <- EH_meth_data[,meth_columns_tmp]
EH_meth_values <- EH_meth_data[meth_sites_names_tmp]
EH_age <- EH_meth_data$age

# # if males ONLY
# EH_meth_values <- EH_meth_values[EH_sex == "M",]
# EH_age <- EH_age[EH_sex == "M"]

# if only fish younger than 0.25 relative age
EH_meth_values <- EH_meth_values[EH_age <= 5,]
EH_sex <- EH_metadata_samples$sex[EH_age <= 5]
EH_age <- EH_age[EH_age <= 5]

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

### saving methyl_sites
##AC
write.csv(AC_methyl_sites, file = paste0(save_dir, "HS_AC_methyl_sites.csv") )
save(AC_methyl_sites, file = paste0(save_dir, "HS_AC_methyl_sites.Rdata"))

##AS
write.csv(AS_methyl_sites, file = paste0(save_dir, "HS_AS_methyl_sites.csv") )
save(AS_methyl_sites, file = paste0(save_dir, "HS_AS_methyl_sites.Rdata"))

##EH
write.csv(EH_methyl_sites, file = paste0(save_dir, "HS_EH_methyl_sites.csv") )
save(EH_methyl_sites, file = paste0(save_dir, "HS_EH_methyl_sites.Rdata"))

##ZF
write.csv(ZF_methyl_sites, file = paste0(save_dir, "HS_ZF_methyl_sites.csv") )
save(ZF_methyl_sites, file = paste0(save_dir, "HS_ZF_methyl_sites.Rdata"))

### age metadata
save(AC_age, AS_age, EH_age, ZF_age, file = paste0(save_dir, "HS_all_age.Rdata"))
#### Imputation ####

# There are several ways of imputing missing values, here we present two of them. Always set.seed() for imputations.
#a) Method 1 using package “mice” (Multiple Imputation by Chained Equation)
# library(mice)
# 
# set.seed(123)
# 
# init <- mice(ZF_meth_values, maxit=0)
# 
# m_method <- init$method
# 
# pred_matrix <- init$predictorMatrix
# 
# colnames(ZF_meth_values)
# predM[, c("age")]=0
# meth[c("age")]=""
# 
# imputed <- mice(ZF_meth_values, method = m_method, predictorMatrix = pred_matrix, m = 5)

#b) Method 2 using package “zoo” (Missing values replaced by the mean or other function of its group)
library(zoo)
set.seed(123)
ZF_meth_values_imputed <- na.aggregate(ZF_meth_values)

##ZF
write.csv(ZF_meth_values_imputed, file = paste0(save_dir, "HS_ZF_meth_values_imputed.csv"))
save(ZF_meth_values_imputed, file = paste0(save_dir, "HS_ZF_meth_values_imputed.Rdata"))

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
PCA_EH <- prcomp(EH_meth_values, scale = TRUE)
PCA_values_EH <- as.data.frame(PCA_EH$x)
PCA_values_EH$species <- "EH"
PCA_values_EH$sex <- EH_sex

##ZF
PCA_ZF <- prcomp(ZF_meth_values_imputed,scale = TRUE)
PCA_values_ZF <- as.data.frame(PCA_ZF$x)
PCA_values_ZF$species <- "ZF"
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
  geom_point(cex = 3, aes(shape = sex)) +
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


#### Plotting Methylation values ####
## max age span modifier
AC_max_age_mod <- 1.3 #30% more

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
EH_meth_values_long$age <- rep(EH_age, each = ncol(EH_meth_values))
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
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values EH (human rgenome)")

ggplot(ZF_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
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