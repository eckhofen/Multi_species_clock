#### Overview ####

# plotting methylation values 

#### Preparation ####
library(GenomicRanges) # https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html
library(GenomicAlignments)
library(tidyverse)
library(ggforce)

##AC
load("/workspace/cfngle/results-data/05_shared_methyl_values/AC_meth_values_v_2.Rdata")

##AS
load("/workspace/cfngle/results-data/05_shared_methyl_values/AS_meth_values_v_2.Rdata")

##EH
load("/workspace/cfngle/results-data/05_shared_methyl_values/EH_meth_values_v_2.Rdata")

##JM
load("/workspace/cfngle/results-data/05_shared_methyl_values/JM_meth_values_v_2.Rdata")

##ZF
load("/workspace/cfngle/results-data/05_shared_methyl_values/ZF_meth_values_v_2.Rdata")


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

##ZF
PCA_ZF <- prcomp(ZF_meth_values,scale = TRUE)
PCA_values_ZF <- as.data.frame(PCA_ZF$x)
PCA_values_ZF$species <- "ZF"

#### plotting results ####

## color palettes
colpalOI <- palette.colors(palette = "Okabe-Ito") %>% 
  as.vector()
colpal <- hcl.colors(7, "SunsetDark") 
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

##AC
ggplot(PCA_values_AC, aes(x = PC1, y = PC2, color = AC_age, shape = species)) +
  geom_point(aes(cex = 7)) +
  facet_wrap(~species) +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6]) +
  geom_text(aes(label = 1:nrow(AC_meth_data)), nudge_x = 0.2, nudge_y = 0) +
  theme_classic()

##AS
ggplot(PCA_values_AS, aes(x = PC3, y = PC2, color = AS_age, shape = species)) +
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

AS_meth_values_long <- pivot_longer(AS_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AS_meth_values_long$age <- rep(AS_age, times = ncol(AS_meth_values))
AS_meth_values_long$SMR <- as.factor(rep(AS_methyl_sites$SMR, times = length(AS_age)))
AS_meth_values_long$Site_i <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.integer()
AS_meth_values_long$Site_f <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.factor()
AS_meth_values_long$pos_rgenome <- rep(AS_methyl_sites$pos_rgenome, times = length(AS_age))
AS_meth_values_long$species <- "AS"

AC_meth_values_long <- pivot_longer(AC_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AC_meth_values_long$age <- rep(AC_age, times = ncol(AC_meth_values))
AC_meth_values_long$SMR <- as.factor(rep(AC_methyl_sites$SMR[meth_sites_names_tmp_AC %in% AC_meth_data_test], times = length(AC_age))) # indexing is necessary because not all CpGs were able to be extracted from the shared sites due o batch correction
AC_meth_values_long$Site_i <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.integer()
AC_meth_values_long$Site_f <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.factor()
AC_meth_values_long$pos_rgenome <- rep(AC_methyl_sites$pos_rgenome[meth_sites_names_tmp_AC %in% AC_meth_data_test], times = length(AC_age))
AC_meth_values_long$species <- "AC"

EH_meth_values_long <- pivot_longer(EH_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
EH_meth_values_long$age <- rep(EH_meth_data$age, times = ncol(EH_meth_values))
EH_meth_values_long$SMR <- as.factor(rep(EH_methyl_sites$SMR, times = length(EH_age)))
EH_meth_values_long$Site_i <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.integer()
EH_meth_values_long$Site_f <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.factor()
EH_meth_values_long$pos_rgenome <- rep(EH_methyl_sites$pos_rgenome, times = length(EH_age))
EH_meth_values_long$species <- "EH"
# EH_meth_values_long$id <- rep(paste0("EH_", EH_meth_data$id), times = ncol(EH_meth_values))

JM_meth_values_long <- pivot_longer(JM_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
JM_meth_values_long$age <- rep(JM_meth_data$age, times = ncol(JM_meth_values))
JM_meth_values_long$SMR <- as.factor(rep(JM_methyl_sites$SMR, times = length(JM_age)))
JM_meth_values_long$Site_i <- gsub(".*\\:", "", JM_meth_values_long$Site) %>% as.integer()
JM_meth_values_long$Site_f <- gsub(".*\\:", "", JM_meth_values_long$Site) %>% as.factor()
JM_meth_values_long$pos_rgenome <- rep(JM_methyl_sites$pos_rgenome, times = length(JM_age))
JM_meth_values_long$species <- "JM"
# C_meth_values_long$sample <- rep(paste0("sample_", 1:length(AC_meth_values)),each = ncol(AC_age))

ZF_meth_values_long <- pivot_longer(ZF_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
ZF_meth_values_long$age <- rep(ZF_meth_data$age, times = ncol(ZF_meth_values))
ZF_meth_values_long$SMR <- as.factor(rep(ZF_methyl_sites$SMR, times = length(ZF_age)))
ZF_meth_values_long$Site_i <- gsub(".*\\:", "", ZF_meth_values_long$Site) %>% as.integer()
ZF_meth_values_long$Site_f <- gsub(".*\\:", "", ZF_meth_values_long$Site) %>% as.factor()
ZF_meth_values_long$pos_rgenome <- rep(ZF_methyl_sites$pos_rgenome, times = length(ZF_age))
ZF_meth_values_long$species <- "ZF"

all_meth_values_long <- rbind(AC_meth_values_long, AS_meth_values_long, EH_meth_values_long, JM_meth_values_long, ZF_meth_values_long)

AC_meth_positions <- data.frame(species = "AC",
                                SMR = AC_methyl_sites$SMR[meth_sites_names_tmp_AC %in% AC_meth_data_test],
                                position = AC_methyl_sites$pos_align[meth_sites_names_tmp_AC %in% AC_meth_data_test])

AS_meth_positions <- data.frame(species = "AS",
                                SMR = AS_methyl_sites$SMR,
                                position = AS_methyl_sites$pos_align)

EH_meth_positions <- data.frame(species = "EH",
                                SMR = EH_methyl_sites$SMR,
                                position = EH_methyl_sites$pos_align)

JM_meth_positions <- data.frame(species = "JM",
                                SMR = JM_methyl_sites$SMR,
                                position = JM_methyl_sites$pos_align)

ZF_meth_positions <- data.frame(species = "ZF",
                                SMR = ZF_methyl_sites$SMR,
                                position = ZF_methyl_sites$pos_align)


all_meth_positions <- rbind(AC_meth_positions, AS_meth_positions, EH_meth_positions, JM_meth_positions, ZF_meth_positions)

## plotting methyl values

ggplot(ZF_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = colpal[1], high = colpal[7], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values AS")
  
ggplot(AC_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age, cex = 0.5, alpha = 0.7)) +
  geom_boxplot(aes(group = Site_f, alpha = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_gradient(low = cbbPalette[2], high = cbbPalette[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values JM")

## all

ggplot(all_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina() +
  geom_boxplot(aes(group = Site_f, fill = species), alpha = 0.9) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_fill_manual(values = colpalOI[c(-1,-9)]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values Cod (AC), Snapper (AS), Hake (EH), Medaka (JM), Zebrafish (ZF)")

# with rgenome position
ggplot(all_meth_values_long, aes(y = Methylation_Value)) +
  geom_sina(aes(x = Site, color = species)) +
  # geom_boxplot(aes(x = Site, group = Site_f, fill = species), alpha = 0.9) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_fill_manual(values = colpalOI[c(-1,-9)]) +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Methylation values Cod (AC), Snapper (AS), Hake (EH), Medaka (JM), Zebrafish (ZF)")

## testing

ggplot(all_meth_positions) +
  geom_point(aes(x = position, y = species, color = species, cex = 0.5)) +
  facet_wrap(~SMR, scale = "free_x") +
  scale_color_manual(values = colpalOI[c(-1,-9)]) +
  theme_classic() +
  labs(title = "Methylation positions Cod (AC), Snapper (AS), Hake (EH), Medaka (JM), Zebrafish (ZF)") +
  xlab("CpG position on rgenome")

