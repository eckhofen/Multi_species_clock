#### Overview ####
# Data comparison

#### Settings ####
# change working directory accordingly extension
# setwd("/powerplant/workspace/cfngle/script_GH/Multi_species_clock/")
setwd("/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/Multi_species_clock/")
extension <- ".pdf"

# setting up color palette 
colpal_CB_c <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

color_species_df <- data.frame(species = as.factor(c("AC","AS","EH","JM","ZF")), color = colpal_CB_c[c(1, 5, 3, 7, 8)])
color_species <- setNames(color_species_df$color, color_species_df$species)

color_compare <- c("#005AB5", "#DC3220")

#### Preparation ####
library(tidyverse)
library(tidymodels)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggpattern)
library(patchwork)
library(caret)
library(glmnet)
library(svglite)

#### load data ####
load("000_data/007_data_comparison/mlm_age_summary.Rdata")
load("000_data/007_data_comparison/mlm_rel_age_summary.Rdata")


#### data manipulation ####
# define significance threshold
significance <- 0.05

##select significant regression coefficients 
mlm_age_sign_coef <- mlm_age_summary$coefficients[mlm_age_summary$coefficients[,4] < significance,c(1,4)] %>% 
  .[-1,] #getting rid of intercept
mlm_rel_age_sign_coef <- mlm_rel_age_summary$coefficients[mlm_rel_age_summary$coefficients[,4] < significance,c(1,4)] %>% 
  .[-1,] #getting rid of intercept

df_mlm_age_sign_coef <- data.frame(model = "Chronological age", 
                                   reg_coef = mlm_age_sign_coef[,1], 
                                   reg_coef_scaled = mlm_age_sign_coef[,1]/max(abs(mlm_age_sign_coef[,1])), 
                                   p_value = mlm_age_sign_coef[,2], 
                                   SMR = rownames(mlm_age_sign_coef))
df_mlm_rel_age_sign_coef <- data.frame(model = "Relative age", 
                                       reg_coef = mlm_rel_age_sign_coef[,1], 
                                       reg_coef_scaled = mlm_rel_age_sign_coef[,1]/max(abs(mlm_rel_age_sign_coef[,1])),
                                       p_value = mlm_rel_age_sign_coef[,2], 
                                       SMR = rownames(mlm_rel_age_sign_coef))

# only significant
df_SMR_comparison <- rbind(df_mlm_age_sign_coef[-1,], df_mlm_rel_age_sign_coef[-1,])
# add significance
df_SMR_comparison <- df_SMR_comparison %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ "'",
    TRUE ~ ""
  ))

### all SMRs
df_mlm_age <- data.frame(model = "Chronological age", 
                         reg_coef = mlm_age_summary$coefficients[-1,1],
                         reg_coef_scaled = mlm_age_summary$coefficients[-1,1]/max(abs(mlm_age_summary$coefficients[-1,1])), 
                         p_value = mlm_age_summary$coefficients[-1,4], 
                         SMR = rownames(mlm_age_summary$coefficients[-1,]))
df_mlm_rel_age <- data.frame(model = "Relative age", 
                             reg_coef = mlm_rel_age_summary$coefficients[-1,1],
                             reg_coef_scaled = mlm_rel_age_summary$coefficients[-1,1]/max(abs(mlm_rel_age_summary$coefficients[-1,1])), 
                             p_value = mlm_rel_age_summary$coefficients[-1,4], 
                             SMR = rownames(mlm_rel_age_summary$coefficients[-1,]))

# only significant
df_SMR_comparison_all <- rbind(df_mlm_age, df_mlm_rel_age)
# colnames(df_SMR_comparison_all) <-c("reg_coef", "p_value", "SMR")

# add significance
df_SMR_comparison_all <- df_SMR_comparison_all %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ "'",
    TRUE ~ ""
  ))


## only significant
plot_SMR_comparison <- ggplot(df_SMR_comparison, aes(x = reg_coef_scaled, y = SMR, fill = model)) +
  geom_point(aes(color = model), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model), position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", y = "Shared methylation regions", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.15, , color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.8) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare) +
  theme_bw() +
  theme(legend.position = "bottom")
plot_SMR_comparison

ggsave(file = paste0("002_plots/006_SMR_comparison",extension), plot_SMR_comparison, width = 5, height = 5)

## all
plot_SMR_comparison_all <- ggplot(df_SMR_comparison_all, aes(x = reg_coef_scaled, y = SMR, fill = model)) +
  geom_point(aes(color = model), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model), position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", y = "Shared methylation regions", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.1, , color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.8) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare) +
  theme_bw() +
  theme(legend.position = "bottom")
plot_SMR_comparison_all

ggsave(file = paste0("002_plots/006_SMR_comparison_all",extension), plot_SMR_comparison_all, width = 5, height = 5)