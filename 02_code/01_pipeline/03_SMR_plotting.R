# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: SMR plotting
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026-02

#### Preparation ####
library(tidyverse)
library(patchwork)

data_folder <- "01_data/"
save_folder <- "01_data/03_intermediate/03_SMR/"
results_folder <- "03_results/01_figures/"

#### Loading data ####

load("01_data/03_intermediate/04_methyl_values/all_meth_values.Rdata")
load("01_data/03_intermediate/04_methyl_values/HS_all_age.Rdata")
load("01_data/03_intermediate/04_methyl_values/all_methyl_sites.Rdata")

load("01_data/04_metadata/color_palettes.RData")

# plot helper
source("02_code/02_helpers/plot_style.R")

# check how many CpGs were captured in the SMRs for each species
length(AC_meth_values) #518
length(AS_meth_values) #281
length(EH_meth_values) #1268
length(ZF_meth_values) #1127

nrow(AC_methyl_sites) #443 (due to batch correction)
nrow(AS_methyl_sites) #281
nrow(EH_methyl_sites) #1268
nrow(ZF_methyl_sites) #1127

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
PCA_EH_sum <- summary(PCA_EH)

##ZF
PCA_ZF <- prcomp(ZF_meth_values, scale = TRUE)
PCA_values_ZF <- as.data.frame(PCA_ZF$x)
PCA_values_ZF$species <- "ZF"
PCA_ZF_sum <- summary(PCA_ZF)


### plotting PCA

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

message("plotting completed")
message("Pipeline completed!")