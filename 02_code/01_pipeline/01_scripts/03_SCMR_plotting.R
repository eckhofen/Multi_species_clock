# Metadata ----------------------------------------------------------------
# Project: Piscine Multispecies Epigenetic Clock
# Description: SCMR plotting
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026-02

#### Preparation ####
library(tidyverse)
library(patchwork)

#### Loading data ####

load("01_data/03_intermediate/04_methyl_values/all_meth_values.Rdata")
load("01_data/03_intermediate/04_methyl_values/HS_all_age.Rdata")

# check how many CpGs were captured in the SMRs for each species
length(AC_meth_values) #518
length(AS_meth_values) #281
length(EH_meth_values) #1268
length(ZF_meth_values) #1127

# in total 
length(c(AC_meth_values, AS_meth_values, EH_meth_values, ZF_meth_values)) # 3194

#### PCA ####
message("computing PCA...")
##AC
PCA_AC <- prcomp(AC_meth_values,scale = TRUE)
PCA_values_AC <- as.data.frame(PCA_AC$x)
PCA_values_AC$species <- "AC"
PCA_AC_sum <- summary(PCA_AC)

##AS
PCA_AS <- prcomp(AS_meth_values,scale = TRUE)
PCA_values_AS <- as.data.frame(PCA_AS$x)
PCA_values_AS$species <- "AS"
PCA_AS_sum <- summary(PCA_AS)

##EH
PCA_EH <- prcomp(EH_meth_values, scale = TRUE)
PCA_values_EH <- as.data.frame(PCA_EH$x)
PCA_values_EH$species <- "EH"
PCA_values_EH$sex <- EH_sex
PCA_EH_sum <- summary(PCA_EH)

##ZF
PCA_ZF <- prcomp(ZF_meth_values, scale = TRUE)
PCA_values_ZF <- as.data.frame(PCA_ZF$x)
PCA_values_ZF$species <- "ZF"
PCA_ZF_sum <- summary(PCA_ZF)


### plotting PCA
# custom theme
theme_custom <- theme_minimal() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(color = "black", linewidth = 1),
            plot.title = element_text(face = "bold", hjust = 0.5), 
            axis.ticks = element_line(color = "black"),
            axis.ticks.length = unit(0.1, "cm"), 
            legend.ticks = element_blank())

# setting theme 
theme_set(theme_custom) 

message("plotting PCAs...")
# color palettes
colpal_CB_c <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

##AC
AC_pca_plot <- ggplot(PCA_values_AC, aes(x = PC1, y = PC2, color = AC_age)) +
  geom_point(cex = 2.5) +
  scale_color_gradient(low = colpal_CB_c[1], high = colpal_CB_c[5]) +
  labs(title = "AC", color = "Age", y = paste0("PC2 (", round(PCA_AC_sum$importance[2,2],4)*100, "%)"), x = paste0("PC1 (", round(PCA_AC_sum$importance[2,1],4)*100, "%)"))

##AS
AS_pca_plot <- ggplot(PCA_values_AS, aes(x = PC1, y = PC2, color = AS_age)) +
  geom_point(cex = 2.5) +
  scale_color_gradient(low = colpal_CB_c[2], high = colpal_CB_c[6]) +
  labs(title = "AS", color = "Age", y = paste0("PC2 (", round(PCA_AS_sum$importance[2,2],4)*100, "%)"), x = paste0("PC1 (", round(PCA_AS_sum$importance[2,1],4)*100, "%)"))

##EH
EH_pca_plot <- ggplot(PCA_values_EH, aes(y = PC2, x = PC1, color = EH_age)) +
  geom_point(cex = 2.5) +
  scale_color_gradient(low = colpal_CB_c[3], high = colpal_CB_c[7]) +
  labs(title = "EH", color = "Age", y = paste0("PC2 (", round(PCA_EH_sum$importance[2,2],4)*100, "%)"), x = paste0("PC1 (", round(PCA_EH_sum$importance[2,1],4)*100, "%)"))

##ZF
ZF_pca_plot <- ggplot(PCA_values_ZF, aes(x = PC1, y = PC2, color = ZF_age)) +
  geom_point(cex = 2.5) +
  scale_color_gradient(low = colpal_CB_c[4], high = colpal_CB_c[8]) +
  labs(title = "ZF", color = "Age", y = paste0("PC2 (", round(PCA_ZF_sum$importance[2,2],4)*100, "%)"), x = paste0("PC1 (", round(PCA_ZF_sum$importance[2,1],4)*100, "%)"))

## plot all 
PCA_plot_all <- (AC_pca_plot + AS_pca_plot + EH_pca_plot + ZF_pca_plot) +
  plot_layout(nrow=2)

ggsave(filename = paste0(results_folder, "03_PCA_all.pdf"), plot = PCA_plot_all, width = 8, height = 7)
ggsave(filename = paste0(results_folder, "03_PCA_all.png"), plot = PCA_plot_all, width = 8, height = 7)

#### Plotting Methylation values ####
message("plotting methylation values...")

# transforming all the values into a plot-friendly dataframe for ggplot2
AS_meth_values_long <- pivot_longer(AS_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AS_meth_values_long$age <- rep(AS_age, each = ncol(AS_meth_values))
AS_meth_values_long$max_age <- 54
AS_meth_values_long$rel_age <- AS_meth_values_long$age / AS_meth_values_long$max_age
AS_region_id <- if ("SCMR" %in% colnames(AS_methyl_sites)) AS_methyl_sites$SCMR else AS_methyl_sites$SMR
AS_meth_values_long$SCMR <- as.factor(rep(AS_region_id, times = length(AS_age)))
AS_meth_values_long$Site_i <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.integer()
AS_meth_values_long$Site_f <- gsub(".*-", "", AS_meth_values_long$Site) %>% as.factor()
AS_meth_values_long$species <- "AS"

AC_meth_values_long <- pivot_longer(AC_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
AC_meth_values_long$age <- rep(AC_age, each = ncol(AC_meth_values))
AC_meth_values_long$max_age <- 25
AC_meth_values_long$rel_age <- AC_meth_values_long$age / AC_meth_values_long$max_age
# AC_meth_values_long$SMR <- as.factor(rep(AC_methyl_sites$SMR[meth_sites_names_tmp_AC %in% AC_meth_data_test], times = length(AC_age))) # indexing is necessary because not all CpGs were able to be extracted from the shared sites due o batch correction
AC_region_id <- if ("SCMR" %in% colnames(AC_methyl_sites)) AC_methyl_sites$SCMR else AC_methyl_sites$SMR
AC_meth_values_long$SCMR <- as.factor(rep(AC_region_id, times = length(AC_age))) # indexing is necessary because not all CpGs were able to be extracted from the shared sites due o batch correction
AC_meth_values_long$Site_i <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.integer()
AC_meth_values_long$Site_f <- gsub(".*\\.", "", AC_meth_values_long$Site) %>% as.factor()
AC_meth_values_long$species <- "AC"

EH_meth_values_long <- pivot_longer(EH_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
EH_meth_values_long$age <- rep(EH_age, each = ncol(EH_meth_values))
EH_meth_values_long$max_age <- 20
EH_meth_values_long$rel_age <- EH_meth_values_long$age / EH_meth_values_long$max_age
EH_region_id <- if ("SCMR" %in% colnames(EH_methyl_sites)) EH_methyl_sites$SCMR else EH_methyl_sites$SMR
EH_meth_values_long$SCMR <- as.factor(rep(EH_region_id, times = length(EH_age)))
EH_meth_values_long$Site_i <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.integer()
EH_meth_values_long$Site_f <- gsub(".*\\.", "", EH_meth_values_long$Site) %>% as.factor()
EH_meth_values_long$species <- "EH"

ZF_meth_values_long <- pivot_longer(ZF_meth_values, cols = everything(), names_to = "Site", values_to = "Methylation_Value")
ZF_meth_values_long$age <- rep(ZF_meth_data$age, each = ncol(ZF_meth_values)) / 52
ZF_meth_values_long$max_age <- 5
ZF_meth_values_long$rel_age <- ZF_meth_values_long$age / ZF_meth_values_long$max_age
ZF_region_id <- if ("SCMR" %in% colnames(ZF_methyl_sites)) ZF_methyl_sites$SCMR else ZF_methyl_sites$SMR
ZF_meth_values_long$SCMR <- as.factor(rep(ZF_region_id, times = length(ZF_age)))
ZF_meth_values_long$Site_i <- gsub(".*\\:", "", ZF_meth_values_long$Site) %>% as.integer()
ZF_meth_values_long$Site_f <- gsub(".*\\:", "", ZF_meth_values_long$Site) %>% as.factor()
ZF_meth_values_long$species <- "ZF"

all_meth_values_long <- rbind(AC_meth_values_long, AS_meth_values_long, EH_meth_values_long, ZF_meth_values_long)
save(all_meth_values_long, AC_meth_values_long, AS_meth_values_long, EH_meth_values_long, ZF_meth_values_long,  file = paste0(save_folder, "all_meth_values_long.RData"))

all_meth_values_long <- all_meth_values_long %>% 
  group_by(SCMR) %>% 
  mutate(pos_nor = Site_i - min(Site_i)) %>% 
  mutate(pos_nor_kb = (Site_i - min(Site_i))/1000) %>% 
  ungroup()

### plotting
save_folder <- paste0(data_folder, "03_intermediate/04_methyl_values/")
load(paste0(save_folder,"all_meth_values_long.RData"))

ggplot(AS_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[1], high = colpal_CB_c[5], guide = "legend") +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values AS (human rgenome)", color = "Age")

ggplot(AC_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[2], high = colpal_CB_c[6], guide = "legend") +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values AC (human rgenome)", color = "Age")

ggplot(EH_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[3], high = colpal_CB_c[7], guide = "legend") +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values EH (human rgenome)", color = "Age")

ggplot(ZF_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[4], high = colpal_CB_c[8], guide = "legend") +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values ZF (human rgenome)", color = "Age")

## all boxplot
ggplot(all_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_boxplot(aes(group = Site_f, fill = species), color = NA, alpha = 0.9, outlier.shape = "") +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_fill_manual(values = color_species) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (AS), European Hake (EH), Zebrafish (ZF) (human rgenome)")

message("plotting completed")