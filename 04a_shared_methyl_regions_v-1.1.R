#### Reducing overlapping ranges of all species and extract the shared (methylation) regions 
#### adding the methylation data to the extracted methylation sites 

library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(GenomicAlignments)
library(tidyverse)

setwd("/workspace/cfngle/results-data/02_conserved_seq")
load("AC_AS_EH_JM_overlaps_bt2.R")
overlap_seqs_bt2

# converting GAlignments object to GRanges to be able to reduce them 
gr_overlap_seqs_bt2 <- lapply(overlap_seqs_bt2, function(x) granges(x))

# intersecting all sequences (option a!) This means that only the sequences which are intersecting will be selected
intersected_seqs_bt2 <- GenomicRanges::intersect(gr_overlap_seqs_bt2[[1]], gr_overlap_seqs_bt2[[2]]) %>% 
  GenomicRanges::intersect(gr_overlap_seqs_bt2[[3]])
boxplot(width(intersected_seqs_bt2))
names(intersected_seqs_bt2) <- sprintf("SMR_c_bt2_%03d", 1:length(intersected_seqs_bt2))
SMR_a_bt2 <- intersected_seqs_bt2

# reducing them is option b of how to define the shared methyl regions
#group them first
group_gr_overlap_bt2 <-  c(gr_overlap_seqs_bt2[[1]], gr_overlap_seqs_bt2[[2]], gr_overlap_seqs_bt2[[3]], gr_overlap_seqs_bt2[[4]])
SMR_b_bt2 <- GenomicRanges::reduce(group_gr_overlap_bt2)
names(SMR_b_bt2) <- sprintf("SMR_b_bt2_%03d", 1:length(SMR_b_bt2))

AC_AS_EH_JM_SMR_b_bt2 <- SMR_b_bt2

save(AC_AS_EH_JM_SMR_b_bt2, file = "/workspace/cfngle/results-data/04_SMRs/AC_AS_EH_JM_SMR_b_bt2.Rdata")

# check which CpGs are in which SMR
hist(subjectHits(findOverlaps(SMR_b_bt2, gr_overlap_seqs_bt2[[3]])))

length(subjectHits(findOverlaps(intersected_seqs_bt2, gr_overlap_seqs_bt2[[3]])))

## getting overlapped sequence data
load("/workspace/cfngle/results-data/02_conserved_seq/AC_AS_EH_JM_overlaps_bt2.R")
load("/workspace/cfngle/results-data/02_conserved_seq/AC_AS_EH_JM_overlaps_mini.R")

overlap_seqs_bt2
overlap_seqs_mini

overlap_AC_AC_bt2 <- overlap_seqs_bt2[[1]]
overlap_AC_AS_bt2 <- overlap_seqs_bt2[[2]]
overlap_AC_EH_bt2 <- overlap_seqs_bt2[[3]]
overlap_AC_JM_bt2 <- overlap_seqs_bt2[[4]]

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
  
  ## finding overlaps between sequences and SMRs
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
                     SMR = paste0("SMR_", SMR_index[i]),
                     species = species)
    df$Chr <- seq_chr_name[i]
    df <- cbind(df[length(df)], df[-length(df)])
    return(df)
  })
  df_final_methyl_sites_AS <- bind_rows(final_methyl_sites_AS)
  return(df_final_methyl_sites_AS)
}

AC_methyl_sites <- get.methyl.sites(overlap_AC_AC_bt2, species = "AC", SMRs = SMR_b_bt2)
AS_methyl_sites <- get.methyl.sites(overlap_AC_AS_bt2, species = "AS", SMRs = SMR_b_bt2)
EH_methyl_sites <- get.methyl.sites(overlap_AC_EH_bt2, species = "EH", SMRs = SMR_b_bt2)
JM_methyl_sites <- get.methyl.sites(overlap_AC_JM_bt2, species = "JM", SMRs = SMR_b_bt2)

SMR_b_bt2

### load methylation data for samples

##AC
xx <- load("/workspace/cfngle/raw-data/AC/zzz_methyl_data/Meth-complete-nobatchcorrection-cod.RData")
assign("AC_meth_data", get(xx))
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
tail(colnames(EH_meth_data))
EH_age <- EH_meth_data$age

##JM
# JM_meth_data <- read.csv("/workspace/cfngle/raw-data/JM/zzz-methyldata/00_JM_methyldata_66079_CpGs.csv")
JM_meth_data <- load("/workspace/cfngle/raw-data/JM/zzz-methyldata/00_JM_methyldata_179818_CpGs.Rdata")
JM_meth_data <- JM_24_methyl_data

tail(colnames(JM_meth_data))
JM_age <- JM_meth_data$age

##ZF ##NOT DONE!!
ZF_meth_data <- read.csv("/workspace/cfngle/raw-data/ZF/")
tail(colnames(ZF_meth_data))
ZF_age <- ZF_meth_data$age

#### extract methylation data for all samples ####

# paste0(gsub("JM_", "",JM_methyl_sites$Chr), ".", JM_methyl_sites$pos_rgenome)

##AC
meth_sites_names_tmp <- paste0(AC_methyl_sites$Chr, ".", AC_methyl_sites$pos_align) %>% 
  gsub("AC_", "Chr", .)
AC_meth_data_test <- gsub("X", "Chr", colnames(AC_meth_data))
meth_sites_names_tmp %in% AC_meth_data_test

# meth_sites_names_tmp %in% colnames(AC_meth_data)

meth_columns_tmp <- sapply(meth_sites_names_tmp, function(x) grep(x, AC_meth_data_test)) %>% 
  as.vector()

AC_meth_values <- AC_meth_data[,meth_columns_tmp]

# AS_meth_values_t <- t(AS_meth_data[,meth_columns_tmp]) %>% 
#   cbind(., data.frame(SMR = AS_methyl_sites$SMR))

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





#### PCA tests and stat tests ####

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

##JM
PCA_JM <- prcomp(JM_meth_values,scale = TRUE)
PCA_values_JM <- as.data.frame(PCA_JM$x)
PCA_values_JM$species <- "JM"


PCA_values_all <- rbind(PCA_values_AC[,c(1:10, length(PCA_values_AC))],
                        PCA_values_AS[,c(1:10, length(PCA_values_AS))],
                        PCA_values_EH[,c(1:10, length(PCA_values_EH))])
age_all <- c(AC_age/mean(AC_age), AS_age/mean(AS_age), EH_age/mean(EH_age))

#### plotting results ####

##AC
ggplot(PCA_values_all, aes(x = PC5, y = PC1, color = age_all, shape = species)) +
  geom_point(aes(cex = 7)) +
  facet_wrap(~species) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_classic()

##EH
ggplot(PCA_values_EH, aes(x = PC1, y = PC2, color = EH_age)) +
  geom_point(aes(cex = 7)) +
  scale_color_gradient(low = "blue", high = "red") +
  theme_classic()


AS_meth_values_long <- pivot_longer(AS_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AS_meth_values_long$age <- rep(AS_age, times = ncol(AS_meth_values))
AS_meth_values_long$SMR <- as.factor(rep(AS_methyl_sites$SMR, times = length(AS_age)))
AS_meth_values_long$Site_i <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.integer()
AS_meth_values_long$Site_f <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.factor()
AS_meth_values_long$species <- "AS"

AC_meth_values_long <- pivot_longer(AC_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AC_meth_values_long$age <- rep(AC_meth_data$age, times = ncol(AC_meth_values))
AC_meth_values_long$SMR <- as.factor(rep(AC_methyl_sites$SMR, times = length(AC_age)))
AC_meth_values_long$Site_i <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.integer()
AC_meth_values_long$Site_f <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.factor()
AC_meth_values_long$species <- "AC"

EH_meth_values_long <- pivot_longer(EH_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
EH_meth_values_long$age <- rep(EH_meth_data$age, times = ncol(EH_meth_values))
EH_meth_values_long$SMR <- as.factor(rep(EH_methyl_sites$SMR, times = length(EH_age)))
EH_meth_values_long$Site_i <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.integer()
EH_meth_values_long$Site_f <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.factor()
EH_meth_values_long$species <- "EH"

JM_meth_values_long <- pivot_longer(JM_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
JM_meth_values_long$age <- rep(JM_meth_data$age, times = ncol(JM_meth_values))
JM_meth_values_long$SMR <- as.factor(rep(JM_methyl_sites$SMR, times = length(JM_age)))
JM_meth_values_long$Site_i <- gsub(".*\\.", "", JM_meth_values_long$Site) %>% as.integer()
JM_meth_values_long$Site_f <- gsub(".*\\.", "", JM_meth_values_long$Site) %>% as.factor()
JM_meth_values_long$species <- "JM"
# C_meth_values_long$sample <- rep(paste0("sample_", 1:length(AC_meth_values)),each = ncol(AC_age))

all_meth_values_long <- rbind(AC_meth_values_long, AS_meth_values_long, EH_meth_values_long, JM_meth_values_long)


ggplot(AS_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_point(aes(color = age, cex = 1, alpha = 0.5)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = "green", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values Snapper")
  
ggplot(AC_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_point(aes(color = age, size = 0.2, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = "green", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values Cod")

ggplot(EH_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_point(aes(color = age, size = 0.2, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = "green", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values EH")

ggplot(JM_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_point(aes(color = age, size = 0.2, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = "green", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values JM")

ggplot(all_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_point(aes(alpha=0.5)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5, fill = species)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = "green", high = "red") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values Cod (AC), Snapper (AS), Hake (EH), Medaka (JM)")



#### testing correlation tests ####
# Initialize a list to store the results
correlation_results <- list()

# Loop through each methylation site
for (i in 1:ncol(AS_meth_values)) {
  site_name <- colnames(AS_meth_values)[i]
  # Perform correlation test with age
  test_result <- cor.test(AS_meth_values[,i], AS_meth_data$age, method = "pearson") # Use "spearman" or "kendall" if more appropriate
  
  # Store the results
  correlation_results[[site_name]] <- list(
    correlation_coefficient = test_result$estimate,
    p_value = test_result$p.value
  )
}

# Optionally, convert the results list to a more convenient format like a dataframe
library(tibble)
correlation_summary <- tibble(
  Site = names(correlation_results),
  Correlation = sapply(correlation_results, function(x) x$correlation_coefficient),
  P_value = sapply(correlation_results, function(x) x$p_value)
)
