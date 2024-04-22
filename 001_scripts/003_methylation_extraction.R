#### Overview ####
# Extraction of methylation data based on overlapping sequences

#### Settings ####
setwd("/powerplant/workspace/cfngle/script_GH/Multi_species_clock/")
setwd("/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/Multi_species_clock/")
data_folder <- paste0(getwd(), "/000_data/")
save_folder <- paste0(data_folder, "003_SMR/") # folder where extracted sequences will be saved

#### Preparation ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(Biostrings) # https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(tidyverse)
library(ggplot2)
library(ggforce)

## loading data 
# loading overlapping alignment reads (rgenome: human) (see script "002....R")
load("000_data/002_conserved_seq/HS_AC_AS_EH_ZF_overlaps.Rdata")

# assign overlapping sequences for each species
overlap_HS_AC <- HS_overlap_seqs[[1]]
overlap_HS_AS <- HS_overlap_seqs[[2]]
overlap_HS_EH <- HS_overlap_seqs[[3]]
overlap_HS_ZF <- HS_overlap_seqs[[4]]

# transforming aligned reads into GRanges object
HS_gr_overlap_seqs <- lapply(HS_overlap_seqs, function(x) granges(x))

# using the overlap of the sequences to get the shared methylation regions (SMRs)
HS_group_gr_overlap <- do.call(c, HS_gr_overlap_seqs)
HS_SMR_b <- GenomicRanges::reduce(HS_group_gr_overlap)
names(HS_SMR_b) <- sprintf("HS_SMR_b_bt2_%03d", 1:length(HS_SMR_b_bt2))

# Saving file
save(HS_SMR_b, file = paste0(save_folder, "HS_SMR.RData"))

#### Methylation extraction ####

## function for extraction
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
#### Getting all species transformed ####

AC_methyl_sites <- get.methyl.sites(overlap_HS_AC, species = "AC", SMRs = HS_SMR_b)
AS_methyl_sites <- get.methyl.sites(overlap_HS_AS, species = "AS", SMRs = HS_SMR_b)
EH_methyl_sites <- get.methyl.sites(overlap_HS_EH, species = "EH", SMRs = HS_SMR_b)
ZF_methyl_sites <- get.methyl.sites(overlap_HS_ZF, species = "ZF", SMRs = HS_SMR_b)

methyl_sites_combined <- bind_rows(AC_methyl_sites, AS_methyl_sites, EH_methyl_sites, ZF_methyl_sites)
methyl_sites_combined$species <- as.factor(methyl_sites_combined$species)
HS_chr_names <- sort(unique(methyl_sites_combined$chr_align))

#### plotting data ####
# adding color palette
colpal_CB_c <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

color_species_df <- data.frame(species = as.factor(c("AC","AS","EH","JM","ZF")), color = colpal_CB_c[c(1, 5, 3, 7, 8)])
color_species <- setNames(color_species_df$color, color_species_df$species)

# visualising the genomic locations and the frequency of CpGs per species
ggplot(methyl_sites_combined) +
  geom_histogram(position = "dodge", aes(x = pos_align, fill = species)) +  
  facet_wrap(~ chr_align, scales = "free_x") +
  labs(x = "Position", y = "CpGs", title = "Number of CpGs and their genomic position") +
  scale_color_manual(values = color_species) + 
  theme_minimal() +
  theme(strip.text.y = element_text(angle = 0), )  

#### Methylation metadata ####
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
EH_age <- EH_meth_data$age # this is the average age and was calculated
EH_metadata_samples <- read.csv("/workspace/cfngle/raw-data/EH/zzz-methyl_data/hake-samples.txt", sep = "\t")
EH_sex <- EH_metadata_samples$sex
EH_age <- EH_metadata_samples$age

##ZF
ZF_meth_data <- load("/workspace/cfngle/raw-data/ZF/zzz_methyldata/ZF_methyldata_88.RData")
ZF_meth_data <- ZF_methyl_data

ZF_age <- ZF_meth_data$age/52

#### Extract methylation data for all samples ####
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

# # if males ONLY
# EH_meth_values <- EH_meth_values[EH_sex == "M",]
# EH_age <- EH_age[EH_sex == "M"]

# # if only fish younger than 0.25 relative age
# EH_meth_values <- EH_meth_values[EH_age <= 5,]
# EH_sex <- EH_metadata_samples$sex[EH_age <= 5]
# EH_age <- EH_age[EH_age <= 5]

##ZF
meth_sites_names_tmp <- paste0(gsub("ZF_", "",ZF_methyl_sites$Chr), ":", ZF_methyl_sites$pos_rgenome)
table(meth_sites_names_tmp %in% colnames(ZF_meth_data))
meth_columns_tmp <- ZF_meth_data[meth_sites_names_tmp] %>% as.vector()
meth_columns_tmp <- unlist(meth_columns_tmp)
# ZF_meth_values_JM <- ZF_meth_data[,meth_columns_tmp] 
ZF_meth_values <- ZF_meth_data[meth_sites_names_tmp]

## saving data
save_folder <- paste0(data_folder, "004_methyl_values/")

### saving methylation VALUES
##AC
write.csv(AC_meth_values, file = paste0(save_folder, "HS_AC_meth_values.csv") )
save(AC_meth_values, file = paste0(save_folder, "HS_AC_meth_values.Rdata"))

##AS
write.csv(AS_meth_values, file = paste0(save_folder, "HS_AS_meth_values.csv") )
save(AS_meth_values, file = paste0(save_folder, "HS_AS_meth_values.Rdata"))

##EH
write.csv(EH_meth_values, file = paste0(save_folder, "HS_EH_meth_values.csv") )
save(EH_meth_values, file = paste0(save_folder, "HS_EH_meth_values.Rdata"))

##ZF CONTAIN NA values, check below for imputation
write.csv(ZF_meth_values, file = paste0(save_folder, "HS_ZF_meth_values.csv") )
save(ZF_meth_values, file = paste0(save_folder, "HS_ZF_meth_values.Rdata"))

### saving methylation SITES
##AC
write.csv(AC_methyl_sites, file = paste0(save_folder, "HS_AC_methyl_sites.csv") )
save(AC_methyl_sites, file = paste0(save_folder, "HS_AC_methyl_sites.Rdata"))

##AS
write.csv(AS_methyl_sites, file = paste0(save_folder, "HS_AS_methyl_sites.csv") )
save(AS_methyl_sites, file = paste0(save_folder, "HS_AS_methyl_sites.Rdata"))

##EH
write.csv(EH_methyl_sites, file = paste0(save_folder, "HS_EH_methyl_sites.csv") )
save(EH_methyl_sites, file = paste0(save_folder, "HS_EH_methyl_sites.Rdata"))

##ZF
write.csv(ZF_methyl_sites, file = paste0(save_folder, "HS_ZF_methyl_sites.csv") )
save(ZF_methyl_sites, file = paste0(save_folder, "HS_ZF_methyl_sites.Rdata"))

### saving age metadata
save(AC_age, AS_age, EH_age, ZF_age, file = paste0(save_folder, "HS_all_age.Rdata"))

#### Imputation ####
library(zoo)
set.seed(123)
ZF_meth_values_imputed <- na.aggregate(ZF_meth_values)

##ZF
write.csv(ZF_meth_values_imputed, file = paste0(save_folder, "HS_ZF_meth_values_imputed.csv"))
save(ZF_meth_values_imputed, file = paste0(save_folder, "HS_ZF_meth_values_imputed.Rdata"))

#### Loading data ####

load("000_data/004_methyl_values/HS_AC_meth_values.Rdata")
load("000_data/004_methyl_values/HS_AS_meth_values.Rdata")
load("000_data/004_methyl_values/HS_EH_meth_values.Rdata")
load("000_data/004_methyl_values/HS_ZF_meth_values_imputed.Rdata")
load("000_data/004_methyl_values/HS_all_age.Rdata")
load("000_data/004_methyl_values/")

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


### plotting PCA
library(patchwork)
# color palettes
colpal_CB_c <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

color_species_df <- data.frame(species = as.factor(c("AC","AS","EH","JM","ZF")), color = colpal_CB_c[c(1, 5, 3, 7, 8)])
color_species <- setNames(color_species_df$color, color_species_df$species)

##AC
AC_pca_plot <- ggplot(PCA_values_AC, aes(x = PC1, y = PC2, color = AC_age)) +
  geom_point(cex = 3) +
  facet_wrap(~species) +
  scale_color_gradient(low = colpal_CB_c[1], high = colpal_CB_c[5]) +
  # geom_text(aes(label = 1:nrow(AC_meth_data)), nudge_x = 0.4, nudge_y = 0) +
  theme_classic()

##AS
AS_pca_plot <- ggplot(PCA_values_AS, aes(x = PC1, y = PC2, color = AS_age)) +
  geom_point(cex = 3) +
  facet_wrap(~species) +
  scale_color_gradient(low = colpal_CB_c[2], high = colpal_CB_c[6]) +
  # geom_text(aes(label = AS_meth_data$id), nudge_x = 0, nudge_y = 0.5) +
  theme_classic()

##EH
EH_pca_plot <- ggplot(PCA_values_EH, aes(x = PC2, y = PC1, color = EH_age)) +
  geom_point(cex = 3) +
  facet_wrap(~species) +
  scale_color_gradient(low = colpal_CB_c[3], high = colpal_CB_c[7]) +
  # geom_text(aes(label = EH_meth_data$id), nudge_x = 0.4, nudge_y = 0) +
  theme_classic()

##ZF
ZF_pca_plot <- ggplot(PCA_values_ZF, aes(x = PC1, y = PC2, color = ZF_age)) +
  geom_point(cex = 3) +
  facet_wrap(~species) +
  scale_color_gradient(low = colpal_CB_c[4], high = colpal_CB_c[8]) +
  # geom_text(aes(label = rownames(ZF_meth_data)), nudge_x = 0.4, nudge_y = 0) +
  theme_classic()

## plot all 
PCA_plot_all <- (AC_pca_plot + AS_pca_plot + EH_pca_plot + ZF_pca_plot) +
  plot_layout(nrow=2)

ggsave(filename = "002_plots/003_PCA_all.pdf", plot = PCA_plot_all, width = 8, height = 7)

#### Plotting Methylation values ####
## max age span modifier
# AC_max_age_mod <- 1.3 #30% more

# transforming all the values into a plot-frindly dataframe for ggplot2
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
# AC_meth_values_long$SMR <- as.factor(rep(AC_methyl_sites$SMR[meth_sites_names_tmp_AC %in% AC_meth_data_test], times = length(AC_age))) # indexing is necessary because not all CpGs were able to be extracted from the shared sites due o batch correction
AC_meth_values_long$SMR <- as.factor(rep(AC_methyl_sites$SMR, times = length(AC_age))) # indexing is necessary because not all CpGs were able to be extracted from the shared sites due o batch correction
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
save(all_meth_values_long, file = paste0(save_folder, "all_meth_values_long.RData"))

### plotting

ggplot(AS_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[1], high = colpal_CB_c[5], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values AS (human rgenome)")

ggplot(AC_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[2], high = colpal_CB_c[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values AC (human rgenome)")

ggplot(EH_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[3], high = colpal_CB_c[7], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values EH (human rgenome)")

ggplot(ZF_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[4], high = colpal_CB_c[8], guide = "legend") +
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
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (AS), European Hake (EH), Zebrafish (ZF) (human rgenome)")

