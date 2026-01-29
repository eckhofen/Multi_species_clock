# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: December 2024

# DISCLAIMER:
# none

#### Overview ####
## Further analysis of SMR data

#### Settings ####
setwd("/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/Multi_species_clock/")
data_folder <- "/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/Multi_species_clock/000_data/"
save_folder <- paste0(data_folder, "009_further_analysis/") # folder where extracted sequences will be saved

## color palette for plotting
# for species comparison
colpal_CB_c <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

color_species_df <- data.frame(species = as.factor(c("AC","AS","EH","JM","ZF")), color = colpal_CB_c[c(1, 5, 3, 7, 8)])
color_species <- setNames(color_species_df$color, color_species_df$species)
# for comparison of two groups
color_compare <- c("#005AB5", "#DC3220")

#### Preparation ####
library(tidyverse)
library(tidymodels)
library(tibble)
library(dplyr)
library(ggplot2)
library(patchwork)

#### Loading data ####
## Methylation values for all species
# all SMRs
load("000_data/006_model_creation/all_meth_values_all_SMR.RData")
# data verification (age is saved as age relative to maximum lifespan of each species)
AC_meth_values_selected_all
AS_meth_values_selected_all
EH_meth_values_selected_all
ZF_meth_values_selected_all
all_meth_values_all_SMR

# selected SMRs (only the ones with same correlation trajectories in all 4 species)
load("000_data/006_model_creation/all_meth_values_selected.RData")
AC_meth_values_selected
AS_meth_values_selected
EH_meth_values_selected
ZF_meth_values_selected

#### Data manipulation ####

## function for extracting slope from a linear model for multiple xs in a dataframe
slope <- function(df, v_dependent, group) {
  
  results <- data.frame(SMR = character(), slope = numeric(), species = character(), stringsAsFactors = FALSE)
  
  groups <- split(df, df[, group])
  
  for (i in names(groups)) {
    group_data <- groups[[i]]
    
    for (column in sort(names(group_data))) {
      if (column != v_dependent & column != group) {
        model <- lm(group_data[, column] ~ group_data[, v_dependent]) 
        df_temp <- data.frame(SMR = as.character(column), slope = coef(model)[2], species = i)
        results <- rbind(results, df_temp)
      }
    }
  }
  return(results)
}  

# getting slopes for all 
slopes_all_SMR <- slope(all_meth_values_all_SMR, "rel_age", "species")

# adding column describing pos. or neg. slope
slopes_all_SMR$cor <- ifelse(slopes_all_SMR$slope > 0, "positive", "negative")

abs_diff_slopes <- slopes_all_SMR
abs_diff_slopes$slope <- abs(abs_diff_slopes[,"slope"])

#### Plotting ####
p_diff_all <- ggplot(slopes_all_SMR, aes(x = species, y = slope, fill = species, group = interaction(species, cor))) +
  geom_boxplot(color = "black", position = position_dodge(width = 0)) +
  scale_fill_manual(values = color_species) +
  labs(x = "Species", y = "Methylation difference (slope)") +
  theme_minimal()
p_diff_all

ggsave("002_plots/008_diff_all.pdf", p_diff_all, width = 8, height = 6)

# absolute difference
p_abs_diff_all <- ggplot(abs_diff_slopes, aes(x = species, y = slope, color = species, fill = species)) +
  geom_boxplot(color = "black") +
  # scale_color_manual(values = color_species) +
  scale_fill_manual(values = color_species) +
  labs(x = "Species", y = "Absolute methylation difference") +
  theme_minimal()
p_abs_diff_all

ggsave("002_plots/008_abs_diff_all.pdf", p_abs_diff_all, width = 8, height = 6)

# plotting values per SMR
p_diff_all_bars <- ggplot(slopes_all_SMR, aes(y = SMR, x = slope, fill = species)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = color_species) +
  labs(x = "Species", y = "Methylation difference (slope)") +
  theme_minimal()
p_diff_all_bars

ggsave("002_plots/008_diff_all_bars.pdf", p_diff_all_bars, width = 6, height = 12)


all_meth_values_all_SMR_long <- all_meth_values_all_SMR %>% 
  pivot_longer(cols = c(1:53), 
               names_to = "SMR", 
               values_to = "methylation")

p_slope_all <-  ggplot(all_meth_values_all_SMR_long, aes(x = rel_age, y = methylation, color = species, fill = species)) +
  geom_smooth(method = lm) +
  geom_point() +
  facet_wrap(vars(SMR)) +
  scale_color_manual(values = color_species) +
  scale_fill_manual(values = color_species) +
  labs(x = "Relative age", y = "CpG methylation (%)") +
  theme_minimal()
p_slope_all

ggsave("002_plots/008_slope_all.pdf", p_slope_all, width = 10, height = 9)

#### Gene Ontology ####
# tables retrieved with obtained gene list on NCBI with DAVID 

go_all_cc <- read.csv("000_data/008_annotation/GO_gene_all_GOterm_CC.csv", sep = "\t")
go_sel_bp <- read.csv("000_data/008_annotation/GO_gene_sel_GOterm_BP.csv", sep = "\t")

go_sel_mf <- read.csv("000_data/008_annotation/GO_gene_sel_GOterm_MF.csv", sep = "\t") %>% as.list()
go_sel_mf$GOTERM_MF_DIRECT <- lapply(go_sel_mf$GOTERM_MF_DIRECT, function(x) unlist(strsplit(x, ","))) %>% 
  lapply(., function(x) unlist(strsplit(x, "~")))





coef(lm(AS_meth_values_selected_all$SMR_050_neg ~ AS_meth_values_selected_all$rel_age))
coef(lm(AC_meth_values_selected_all$rel_age ~ AC_meth_values_selected_all$SMR_050_neg))
coef(lm(EH_meth_values_selected_all$rel_age ~ EH_meth_values_selected_all$SMR_050_neg))
coef(lm(ZF_meth_values_selected_all$rel_age ~ ZF_meth_values_selected_all$SMR_001_neg))

AS_meth_values_selected_all == all_meth_values_all_SMR_long %>% 
  filter(species == "AS", SMR == "SMR_050_neg")


