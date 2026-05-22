# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: SMR analysis and position plots
# eckhofen@icm.csic.es
# Date: 2026-02

#### Settings ####
data_folder <- "01_data/"
intermediate_folder <- "01_data/03_intermediate/06_model_creation/"
dir.create(intermediate_folder, recursive = TRUE, showWarnings = FALSE)
results_folder <- "03_results/01_figures/"

extension <- ".png"
# color palette
load(paste0(data_folder, "04_metadata/color_palettes.RData"))

#### Preparation ####
library(tidyverse)
library(tidymodels)
library(patchwork)

# plot settings
source("02_code/02_helpers/plot_style.R")

#### Loading data ####

load(file = paste0(intermediate_folder, "all_data.Rdata"))
load(file = paste0(intermediate_folder, "models.Rdata"))
load(file = paste0(intermediate_folder, "models_summary.Rdata"))
load("01_data/03_intermediate/06_model_creation/all_meth_values_selected.RData")
load("01_data/03_intermediate/07_data_comparison/methyl_sites_combined_nor.Rdata")
load("01_data/03_intermediate/05_correlation_data/cor_all.RData")

#### Key parameters ####
message("Extracting key parameters...")

significance <- 0.05

# Extract coefficients for chronological age model
mlm_age_coef <- mlm_chron_summary$coefficients[-1, , drop = FALSE]  # Remove intercept
mlm_age_sign_coef <- mlm_age_coef[mlm_age_coef[, 4] < significance, c(1, 4), drop = FALSE]

max_abs_age_sign <- max(abs(mlm_age_sign_coef[, 1]), na.rm = TRUE)
max_abs_age_all <- max(abs(mlm_age_coef[, 1]), na.rm = TRUE)

df_mlm_age_sign_coef <- data.frame(
  model = "Chronological age",
  reg_coef = mlm_age_sign_coef[, 1],
  reg_coef_scaled = mlm_age_sign_coef[, 1] / max_abs_age_sign,
  p_value = mlm_age_sign_coef[, 2],
  SMR = gsub("^SMR_", "SMR_", rownames(mlm_age_sign_coef))
)

df_mlm_age_all <- data.frame(
  model = "Chronological age",
  reg_coef = mlm_age_coef[, 1],
  reg_coef_scaled = mlm_age_coef[, 1] / max_abs_age_all,
  p_value = mlm_age_coef[, 4],
  SMR = gsub("^SMR_", "SMR_", rownames(mlm_age_coef))
)

# Extract coefficients for relative age model
mlm_rel_coef <- mlm_rel_summary$coefficients[-1, , drop = FALSE]  # Remove intercept
mlm_rel_age_sign_coef <- mlm_rel_coef[mlm_rel_coef[, 4] < significance, c(1, 4), drop = FALSE]

max_abs_rel_age_sign <- if (nrow(mlm_rel_age_sign_coef) == 0) NA_real_ else max(abs(mlm_rel_age_sign_coef[, 1]), na.rm = TRUE)
max_abs_rel_age_all <- max(abs(mlm_rel_coef[, 1]), na.rm = TRUE)

df_mlm_rel_age_sign_coef <- data.frame(
  model = "Relative age",
  reg_coef = mlm_rel_age_sign_coef[, 1],
  reg_coef_scaled = mlm_rel_age_sign_coef[, 1] / max_abs_rel_age_sign,
  p_value = mlm_rel_age_sign_coef[, 2],
  SMR = gsub("^SMR_", "SMR_", rownames(mlm_rel_age_sign_coef))
)

df_mlm_rel_age_all <- data.frame(
  model = "Relative age",
  reg_coef = mlm_rel_coef[, 1],
  reg_coef_scaled = mlm_rel_coef[, 1] / max_abs_rel_age_all,
  p_value = mlm_rel_coef[, 4],
  SMR = gsub("^SMR_", "SMR_", rownames(mlm_rel_coef))
)

# Combine significant SMRs
df_SMR_coefficient <- rbind(df_mlm_age_sign_coef, df_mlm_rel_age_sign_coef)
df_SMR_coefficient <- df_SMR_coefficient %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ "'",
    TRUE ~ ""
  ), 
        significant = ifelse(p_value < .05, TRUE, FALSE))

# Combine all SMRs
df_SMR_coefficient_all <- rbind(df_mlm_age_all, df_mlm_rel_age_all)
df_SMR_coefficient_all <- df_SMR_coefficient_all %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ "'",
    TRUE ~ ""
  ),
        significant = ifelse(p_value < .05, TRUE, FALSE))

# Save key parameters
save(file = paste0(intermediate_folder, "df_SMR_coefficient.Rdata"), df_SMR_coefficient)
save(file = paste0(intermediate_folder, "df_SMR_coefficient_all.Rdata"), df_SMR_coefficient_all)

message("Key parameters extracted and saved.")

# Adding metadata to methyl sites
methyl_sites_combined_nor <- 
  inner_join(methyl_sites_combined_nor, cor_all, keep = FALSE) %>%
  na.omit()


#### SMR Comparison Plots ####
message("Creating plots...")

# Check for color_compare palette, create default if not available
if (!exists("color_compare")) {
  color_compare <- c("#005AB5", "#DC3220")
}

### SMR coefficient plots
# Only significant SMRs
p_SMR_coefficient_significant <- ggplot(df_SMR_coefficient, aes(x = reg_coef_scaled, y = SMR, fill = model)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "solid") +
  geom_point(aes(color = model),
    position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model),
    position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", 
       y = "", 
       fill = "Model", 
       color = "Model") +
  geom_text(aes(label = significance, x = ifelse(reg_coef_scaled > 0, reg_coef_scaled + 0.02, reg_coef_scaled - 0.02), color = model, hjust = ifelse(reg_coef_scaled > 0, 0, 1)),
    show.legend = FALSE, position = position_dodge(0.9), vjust = 0.75) +
  scale_fill_manual(values = color_compare) +
  scale_x_continuous(limits = c(-1.1, .6)) +
  scale_color_manual(values = color_compare) +
  theme(legend.position = c(0.15, 0.17), 
        legend.background = element_rect(fill = "grey95", color = NA), 
        legend.title = element_text(face = "bold", hjust = .5))

save_plot(file = paste0(results_folder, "06_SMR_coefficient_significant", extension), p_SMR_coefficient_significant, width = 5, height = 5)

# All SMRs
p_SMR_coefficient_all <- ggplot(df_SMR_coefficient_all, aes(x = reg_coef_scaled, y = SMR, fill = model)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "solid") +
  geom_point(aes(color = model, alpha = significant),
    position = position_dodge(width = 0.9), show.legend = c(alpha = FALSE)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model, alpha = significant),
    position = position_dodge(width = 0.9), show.legend = c(alpha = FALSE)) +
  labs(x = "Normalized regression coefficient", 
       y = "", 
       fill = "Model", 
       color = "Model") +
  geom_text(aes(label = significance, x = ifelse(reg_coef_scaled > 0, reg_coef_scaled + 0.02, reg_coef_scaled - 0.02), color = model, hjust = ifelse(reg_coef_scaled > 0, 0, 1)),
    show.legend = FALSE, position = position_dodge(0.9), vjust = 0.75) +
  scale_fill_manual(values = color_compare) +
  scale_x_continuous(limits = c(-1.1, .6)) +
  scale_color_manual(values = color_compare) +
  theme(legend.position = c(0.15, 0.17), 
        legend.background = element_rect(fill = "grey95", color = NA), 
        legend.title = element_text(face = "bold", hjust = .5))

save_plot(file = paste0(results_folder, "06_SMR_coefficient_all", extension), p_SMR_coefficient_all, width = 8, height = 6)

# boxplot for comparing coefficients between models
p_boxplot <- ggplot(df_SMR_coefficient_all, aes(y = model, x = reg_coef_scaled)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "solid") +
  geom_boxplot(fill = "grey95", outlier.colour = NA, show.legend = FALSE) +
  geom_jitter(aes(color = model, alpha = significant),
                inherit.aes = TRUE, 
                cex = 1.5,
                position = position_jitter(width = 0),
                show.legend = FALSE) +
  labs(y = "", x = "Normalized regression coefficient") +
  scale_color_manual(values = color_compare) +
  scale_x_continuous(limits = c(-1.1, .6)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p_comparison_full <- p_SMR_coefficient_all + p_boxplot +
  plot_layout(nrow = 2, height = c(5,1), axes = "collect") +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

save_plot(file = paste0(results_folder, "06_SMR_comparison_full", extension), p_comparison_full, width = 8, height = 8)

message("Plotting CpG location on SMRs...")

### CpG location for SMRs
# only locations which were used for the model
site_selected <- methyl_sites_combined_nor %>%
  filter(SMR %in% df_SMR_coefficient_all$SMR)

## CpG position for all 
p_SMR_position_all <- ggplot(site_selected) +
  geom_segment(aes(x = pos_nor_kb, xend = pos_nor_kb+.02, y = SMR, color = species, linewidth = selected, alpha = selected), show.legend = c(alpha = FALSE)) +
  labs(x = "Relative length (kb)", y = "", color = "Species") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,7)) +
  scale_color_manual(values = color_species) +
  scale_linewidth_manual(values = c("FALSE" = 3, "TRUE" = 4)) +
  scale_alpha_manual(values = c("FALSE" = .2, "TRUE" = 1)) +
  theme(legend.position = "top", legend.background = element_rect(fill = "grey95", color = NA))

save_plot(filename = paste0(results_folder, "06_SMR_position_all.pdf"), plot = p_SMR_position_all, width = 8, height = 9)

# distribution of CpG positions 
p_SMR_distribtution_all <- ggplot(site_selected) +
  geom_boxplot(aes(x = pos_nor_kb, y = species, fill = species, alpha = selected), 
    show.legend = FALSE, outliers = FALSE) +
  labs(x = "Relative length (kb)", y = "Species") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1), limits = c(0,7)) +
  scale_fill_manual(values = color_species) +
  scale_alpha_manual(values = c("FALSE" = .2, "TRUE" = 1))

save_plot(filename = paste0(results_folder, "06_SMR_distribution_all.pdf"), plot = p_SMR_distribtution_all, width = 7, height = 4)

p_SMR_dist_combined <- p_SMR_position_all + p_SMR_distribtution_all +
  plot_layout(nrow = 2, height = c(4, 1), axes = "collect") &
  theme(legend.position = "none")

save_plot(filename = paste0(results_folder, "06_SMR_distribution_combined.pdf"), plot = p_SMR_dist_combined, width = 6, height = 8)

## only significant
# subsetting for significant SMRs only
site_selected_significant <- methyl_sites_combined_nor %>%
  filter(SMR %in% df_SMR_coefficient$SMR)

p_SMR_position_significant <- ggplot(site_selected_significant) +
  geom_segment(aes(x = pos_nor_kb, xend = pos_nor_kb+.03, y = SMR, color = species, linewidth = selected, alpha = selected), show.legend = FALSE) +
  labs(x = "Relative length (kb)", y = "", color = "Species") +
  scale_x_continuous(limits = c(0, 6.2), breaks = seq(0, 6, 1)) +
  scale_color_manual(values = color_species) +
  scale_linewidth_manual(values = c("FALSE" = 4, "TRUE" = 6)) +
  scale_alpha_manual(values = c("FALSE" = .2, "TRUE" = 1))
 
save_plot(filename = paste0(results_folder, "06_SMR_position_significant.pdf"), plot = p_SMR_position_significant, width = 6, height = 5)

p_SMR_distribtution_significant <- ggplot(site_selected_significant) +
  geom_boxplot(aes(x = pos_nor_kb, y = species, fill = species, alpha = selected), 
    show.legend = FALSE, outliers = FALSE) +
  labs(x = "Relative length (kb)", y = "Species") +
  scale_x_continuous(limits = c(0, 6.2), breaks = seq(0, 6, 1)) +
  scale_fill_manual(values = color_species) +
  scale_alpha_manual(values = c("FALSE" = .2, "TRUE" = 1))

save_plot(filename = paste0(results_folder, "06_SMR_distribution_significant.pdf"), plot = p_SMR_distribtution_significant, width = 7, height = 4)
message("Plotting combined plots...")

# combined
p_SMR_combined <- p_SMR_coefficient_significant + p_SMR_position_significant + p_boxplot + p_SMR_distribtution_significant +
  plot_layout(nrow = 2, height = c(4, 1, 4, 1), axes = "collect_y", guides = "collect", axis_titles = "collect") +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"), 
        legend.position = "bottom", 
        legend.direction = "horizontal")

save_plot(filename = paste0(results_folder, "06_SMR_combined.pdf"), plot = p_SMR_combined, width = 10, height = 8)

## both all together

p_SMR_combined_all <- p_SMR_coefficient_all + p_SMR_position_all +
  plot_layout(guides = "collect", axes = "collect") + 
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

save_plot(filename = paste0(results_folder, "06_SMR_combined_all.pdf"), plot = p_SMR_combined_all, width = 10, height = 6)

#### Distance analysis ####
message("Distance analysis...")

# calcualte relative distance of selected sites
cluster_distance <- methyl_sites_combined_nor %>% group_by(selected, significant, SMR) %>% 
  mutate(cluster_distance = max(pos_nor_kb) - min(pos_nor_kb), 
         mean_correlation = mean(Correlation)) 

p_cluster_distance <- ggplot(cluster_distance) +
  geom_boxplot(aes(y = selected, x = cluster_distance), 
    show.legend = FALSE, outliers = FALSE, fill = "grey95") +
  geom_jitter(aes(y = selected, x = cluster_distance), 
    show.legend = FALSE, cex = 1, shape = 16, position = position_jitter(width = 0)) +
  labs(y = "Selected for model", x = "Cluster distance (kb)", fill = "Significant age correlation") +
  scale_fill_manual(values = color_compare)

save_plot(filename = paste0(results_folder, "06_cluster_distance.pdf"), plot = p_cluster_distance, width = 7, height = 4)

message("SMR comparison plots created and saved.")

#### SMR length distribution ####

# helper functions
summarise_smr <- function(df) {
  df %>%
    group_by(SMR) %>%
    summarise(
      length_bp = max(pos_align) - min(pos_align),
      .groups = "drop"
    ) %>%
    summarise(
      n_SMRs = n(),
      min_length = min(length_bp),
      max_length = max(length_bp),
      mean_length = mean(length_bp),
      median_length = median(length_bp),
      sd_length = sd(length_bp)
    )
}

# All in SMR
smr_lengths_all <- summarise_smr(methyl_sites_combined_nor)

# Selected in SMR
smr_lengths_selected <- summarise_smr(site_selected %>% filter (selected == TRUE))
