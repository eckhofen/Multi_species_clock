# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: Data comparisons and combined plots
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026 05

#### Settings ####

extension <- ".pdf"
data_folder <- "01_data/"
results_folder <- "03_results/01_figures/"
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

load("01_data/04_metadata/color_palettes.RData")

#### Preparation ####
library(tidyverse)
library(tidymodels)
library(ggplot2)
library(patchwork)

#### load data ####
load("01_data/03_intermediate/06_model_creation/mlm_chron_summary.Rdata")
load("01_data/03_intermediate/06_model_creation/mlm_rel_summary.Rdata")
load("01_data/03_intermediate/07_data_comparison/methyl_sites_combined_nor.Rdata")
load("01_data/03_intermediate/05_correlation_data/all_mix_cor_CpG_common.RData")

# plot helper
source("02_code/02_helpers/plot_style.R")

#### data manipulation ####
# define significance threshold
significance <- 0.05

##select significant regression coefficients 
mlm_age_sign_coef <- mlm_chron_summary$coefficients[mlm_chron_summary$coefficients[, 4] < significance, c(1, 4), drop = FALSE]
mlm_rel_age_sign_coef <- mlm_rel_summary$coefficients[mlm_rel_summary$coefficients[, 4] < significance, c(1, 4), drop = FALSE]

mlm_age_sign_coef <- mlm_age_sign_coef[rownames(mlm_age_sign_coef) != "(Intercept)", , drop = FALSE]
mlm_rel_age_sign_coef <- mlm_rel_age_sign_coef[rownames(mlm_rel_age_sign_coef) != "(Intercept)", , drop = FALSE]

max_abs_age_sign <- if (nrow(mlm_age_sign_coef) == 0) NA_real_ else max(abs(mlm_age_sign_coef[, 1]), na.rm = TRUE)
max_abs_rel_age_sign <- if (nrow(mlm_rel_age_sign_coef) == 0) NA_real_ else max(abs(mlm_rel_age_sign_coef[, 1]), na.rm = TRUE)

df_mlm_age_sign_coef <- data.frame(
  model = "Chronological age",
  reg_coef = mlm_age_sign_coef[, 1],
  reg_coef_scaled = ifelse(is.na(max_abs_age_sign) || max_abs_age_sign == 0, NA_real_, mlm_age_sign_coef[, 1] / max_abs_age_sign),
  p_value = mlm_age_sign_coef[, 2],
  SMR = gsub("^SMR_", "SMR_", rownames(mlm_age_sign_coef))
)
df_mlm_rel_age_sign_coef <- data.frame(
  model = "Relative age",
  reg_coef = mlm_rel_age_sign_coef[, 1],
  reg_coef_scaled = ifelse(is.na(max_abs_rel_age_sign) || max_abs_rel_age_sign == 0, NA_real_, mlm_rel_age_sign_coef[, 1] / max_abs_rel_age_sign),
  p_value = mlm_rel_age_sign_coef[, 2],
  SMR = gsub("^SMR_", "SMR_", rownames(mlm_rel_age_sign_coef))
)

# only significant
df_SMR_comparison <- rbind(df_mlm_age_sign_coef, df_mlm_rel_age_sign_coef)
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
                         reg_coef = mlm_chron_summary$coefficients[-1,1],
                         reg_coef_scaled = mlm_chron_summary$coefficients[-1,1]/max(abs(mlm_chron_summary$coefficients[-1,1])), 
                         p_value = mlm_chron_summary$coefficients[-1,4], 
                         SMR = gsub("^SMR_", "SMR_", rownames(mlm_chron_summary$coefficients[-1,])))
df_mlm_rel_age <- data.frame(model = "Relative age", 
                             reg_coef = mlm_rel_summary$coefficients[-1,1],
                             reg_coef_scaled = mlm_rel_summary$coefficients[-1,1]/max(abs(mlm_rel_summary$coefficients[-1,1])), 
                             p_value = mlm_rel_summary$coefficients[-1,4], 
                             SMR = gsub("^SMR_", "SMR_", rownames(mlm_rel_summary$coefficients[-1,])))

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

#### Plots ####
theme_custom <- theme_minimal() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(color = "black", linewidth = 1),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 12), 
            axis.ticks = element_line(color = "black"),
            axis.ticks.length = unit(0.1, "cm"), 
            axis.title = element_text(size = 10), 
            legend.ticks = element_blank())

theme_set(theme_custom)

### SMR comparison
## only significant
plot_SMR_comparison <- ggplot(df_SMR_comparison, aes(x = reg_coef_scaled, y = SMR, fill = model)) +
  geom_point(aes(color = model), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model), position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", y = "Shared methylation region", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.15, color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.8) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare) 
plot_SMR_comparison

save_plot(file = paste0(results_folder, "06_SMR_comparison",extension), plot_SMR_comparison, width = 5, height = 5)

## all
plot_SMR_comparison_all <- ggplot(df_SMR_comparison_all, aes(x = reg_coef_scaled, y = SMR, fill = model)) +
  geom_point(aes(color = model), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model), position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", y = "Shared methylation region", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.1, color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.75) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare)
plot_SMR_comparison_all

save_plot(file = paste0(results_folder, "06_SMR_comparison_all",extension), plot_SMR_comparison_all, width = 10, height = 7)

### CpG location for SMRs
## CpG position for all 
plot_SMR_position_all <- ggplot(methyl_sites_combined_nor) +
  geom_segment(aes(x = pos_nor_kb, xend = pos_nor_kb+.02, y = SMR, color = species), linewidth = 2) +
  labs(x = "CpG position (kb)", y = "Shared methylation region", color = "Species") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_color_manual(values = color_species)

plot_SMR_position_all

save_plot(filename = paste0(results_folder, "06_SMR_position_all.pdf"), plot = plot_SMR_position_all, width = 10, height = 7)

# distribution of CpG positions 
plot_SMR_distribtution <- ggplot(methyl_sites_combined_nor) +
  geom_segment(aes(x = pos_nor_kb, xend = pos_nor_kb+.02, y = SMR, color = species), linewidth = 2) +
  labs(x = "CpG position (kb)", y = "Shared methylation region", color = "Species") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_color_manual(values = color_species)

plot_SMR_position_all

save_plot(filename = paste0(results_folder, "06_SMR_position_all.pdf"), plot = plot_SMR_position_all, width = 10, height = 7)

## only significant
# subsetting for significant SMRs only
methyl_sites_significant <- subset(methyl_sites_combined_nor, SMR %in% df_SMR_comparison$SMR)
 
if (nrow(methyl_sites_significant) > 0) {
  plot_SMR_position_significant <- ggplot(methyl_sites_significant) +
    geom_segment(aes(x = pos_nor_kb, xend = pos_nor_kb+.03, y = SMR, color = species), linewidth = 5) +
    labs(x = "CpG position (kb)", y = "Shared co-methylation region", color = "Species") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_color_manual(values = color_species)
 
  save_plot(filename = paste0(results_folder, "06_SMR_position_significant.pdf"), plot = plot_SMR_position_significant, width = 6, height = 5)
} else {
  plot_SMR_position_significant <- patchwork::plot_spacer()
}

## combined

plot_SMR_combined <- plot_SMR_comparison + plot_SMR_position_significant +
  plot_layout(guides = "collect", axes = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

save_plot(filename = paste0(results_folder, "06_SMR_combined.pdf"), plot = plot_SMR_combined, width = 10, height = 6)

## combined all together

plot_SMR_combined_all <- plot_SMR_comparison_all + plot_SMR_position_all +
  plot_layout(guides = "collect", axes = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

save_plot(filename = paste0(results_folder, "06_SMR_combined_all.pdf"), plot = plot_SMR_combined_all, width = 10, height = 6)
