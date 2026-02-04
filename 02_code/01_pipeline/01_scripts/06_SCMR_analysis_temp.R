# Metadata ----------------------------------------------------------------
# Project: Piscine Multispecies Epigenetic Clock
# Description: Epigenetic model creation and testing
# Author: Gabriel Ecker-Eckhofen
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

#### Loading data ####

load(file = paste0(intermediate_folder, "all_data.Rdata"))
load(file = paste0(intermediate_folder, "models.Rdata"))
load(file = paste0(intermediate_folder, "models_summary.Rdata"))
load("01_data/03_intermediate/06_model_creation/all_meth_values_selected.RData")
load("01_data/03_intermediate/07_data_comparison/methyl_sites_combined_nor.Rdata")
load("01_data/03_intermediate/05_correlation_data/all_mix_cor_CpG_common.RData")

#### Key parameters ####
message("Extracting key parameters...")

significance <- 0.05

# Extract coefficients for chronological age model
mlm_age_coef <- mlm_chron_summary$coefficients[-1, , drop = FALSE]  # Remove intercept
mlm_age_sign_coef <- mlm_age_coef[mlm_age_coef[, 4] < significance, c(1, 4), drop = FALSE]

max_abs_age_sign <- if (nrow(mlm_age_sign_coef) == 0) NA_real_ else max(abs(mlm_age_sign_coef[, 1]), na.rm = TRUE)
max_abs_age_all <- max(abs(mlm_age_coef[, 1]), na.rm = TRUE)

df_mlm_age_sign_coef <- data.frame(
  model = "Chronological age",
  reg_coef = mlm_age_sign_coef[, 1],
  reg_coef_scaled = ifelse(is.na(max_abs_age_sign) || max_abs_age_sign == 0, NA_real_, mlm_age_sign_coef[, 1] / max_abs_age_sign),
  p_value = mlm_age_sign_coef[, 2],
  SCMR = gsub("^SCMR_", "SCMR_", rownames(mlm_age_sign_coef))
)

df_mlm_age_all <- data.frame(
  model = "Chronological age",
  reg_coef = mlm_age_coef[, 1],
  reg_coef_scaled = mlm_age_coef[, 1] / max_abs_age_all,
  p_value = mlm_age_coef[, 4],
  SCMR = gsub("^SCMR_", "SCMR_", rownames(mlm_age_coef))
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
  SCMR = gsub("^SCMR_", "SCMR_", rownames(mlm_rel_age_sign_coef))
)

df_mlm_rel_age_all <- data.frame(
  model = "Relative age",
  reg_coef = mlm_rel_coef[, 1],
  reg_coef_scaled = mlm_rel_coef[, 1] / max_abs_rel_age_all,
  p_value = mlm_rel_coef[, 4],
  SCMR = gsub("^SCMR_", "SCMR_", rownames(mlm_rel_coef))
)

# Combine significant SCMRs
df_SCMR_comparison <- rbind(df_mlm_age_sign_coef, df_mlm_rel_age_sign_coef)
df_SCMR_comparison <- df_SCMR_comparison %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ "'",
    TRUE ~ ""
  ), 
        significant = ifelse(p_value < .05, TRUE, FALSE))

# Combine all SCMRs
df_SCMR_comparison_all <- rbind(df_mlm_age_all, df_mlm_rel_age_all)
df_SCMR_comparison_all <- df_SCMR_comparison_all %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ "'",
    TRUE ~ ""
  ),
        significant = ifelse(p_value < .05, TRUE, FALSE))

# Save key parameters
save(file = paste0(intermediate_folder, "df_SCMR_comparison.Rdata"), df_SCMR_comparison)
save(file = paste0(intermediate_folder, "df_SCMR_comparison_all.Rdata"), df_SCMR_comparison_all)

message("Key parameters extracted and saved.")
#-------------------------
# Adding SCMR selection
methyl_sites_combined_nor$selected <- 
  ifelse(methyl_sites_combined_nor$seq_names %in% df_SCMR_comparison$SCMR, TRUE, FALSE)

#-------------------------
#### SCMR Comparison Plots ####
message("Creating plots...")

# Check for color_compare palette, create default if not available
if (!exists("color_compare")) {
  color_compare <- c("#005AB5", "#DC3220")
}

### SCMR comparison plots
# Only significant SCMRs
p_SCMR_comparison <- ggplot(df_SCMR_comparison, aes(x = reg_coef_scaled, y = SCMR, fill = model)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "solid") +
  geom_point(aes(color = model), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model), position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", y = "Shared co-methylation region", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.15, color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.8) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare)

ggsave(file = paste0(results_folder, "06_SCMR_comparison_significant", extension), p_SCMR_comparison, width = 5, height = 5)

# All SCMRs
p_SCMR_comparison_all <- ggplot(df_SCMR_comparison_all, aes(x = reg_coef_scaled, y = SCMR, fill = model)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "solid") +
  geom_point(aes(color = model, alpha = significant), position = position_dodge(width = 0.9), show.legend = c(alpha = FALSE)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model, alpha = significant), position = position_dodge(width = 0.9), show.legend = c(alpha = FALSE)) +
  labs(x = "Normalized regression coefficient", y = "Shared co-methylation region", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.1, color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.75) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare) +
  theme(legend.position = c(0.15, 0.17), 
        legend.background = element_rect(fill = "grey95", color = NA), 
        legend.title = element_text(face = "bold", hjust = .5))

ggsave(file = paste0(results_folder, "06_SCMR_comparison_all", extension), p_SCMR_comparison_all, width = 8, height = 6)

# boxplot for comparing coefficients between models
p_boxplot <- ggplot(df_SCMR_comparison_all, aes(y = model, x = reg_coef_scaled)) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "solid") +
  geom_boxplot(fill = "grey95", outlier.colour = NA, show.legend = FALSE) +
  geom_jitter(aes(color = model, alpha = significant),
                inherit.aes = TRUE, 
                cex = 1.5,
                position = position_jitter(width = .05),
                show.legend = FALSE) +
  labs(y = "", x = "Normalized regression coefficient") +
  scale_color_manual(values = color_compare) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

p_comparison_full <- p_SCMR_comparison_all + p_boxplot +
  plot_layout(nrow = 2, height = c(4,1), axes = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(file = paste0(results_folder, "06_SCMR_comparison_full", extension), p_comparison_full, width = 8, height = 8)

message("Plotting CpG location on SCMRs...")

### CpG location for SCMRs
## CpG position for all 
p_SCMR_position_all <- ggplot(methyl_sites_combined_nor) +
  geom_segment(aes(x = pos_nor_kb, xend = pos_nor_kb+.02, y = SCMR, color = species), linewidth = 2) +
  labs(x = "CpG position (kb)", y = "Shared co-methylation region", color = "Species") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  scale_color_manual(values = color_species)

p_SCMR_position_all

ggsave(filename = paste0(results_folder, "06_SCMR_position_all.pdf"), plot = p_SCMR_position_all, width = 10, height = 7)

# distribution of CpG positions 
p_SCMR_distribtution <- ggplot(methyl_sites_combined_nor) +
  geom_jitter(aes(x = pos_nor_kb, y = species, alpha = selected), 
    show.legend = FALSE, cex = 1, shape = 16) +
  geom_violin(aes(x = pos_nor_kb, y = species, fill = species), 
    alpha = .4, show.legend = FALSE, outliers = FALSE) +
  labs(x = "CpG position (kb)", y = "Species") +
  scale_fill_manual(values = color_species) +
  scale_color_manual(values = color_species) +
  scale_alpha_manual(values = c(1, .2))

ggsave(filename = paste0(results_folder, "06_SCMR_distribution_all.pdf"), plot = p_SCMR_distribtution, width = 7, height = 4)

## only significant
# subsetting for significant SCMRs only
methyl_sites_significant <- subset(methyl_sites_combined_nor, SCMR %in% df_SCMR_comparison$SCMR)
 

p_SCMR_position_significant <- ggplot(methyl_sites_significant) +
    geom_segment(aes(x = pos_nor_kb, xend = pos_nor_kb+.03, y = SCMR, color = species), linewidth = 5) +
    labs(x = "CpG position (kb)", y = "Shared co-methylation region", color = "Species") +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_color_manual(values = color_species)
 
ggsave(filename = paste0(results_folder, "06_SCMR_position_significant.pdf"), plot = p_SCMR_position_significant, width = 6, height = 5)

message("Plotting combined plots...")
# combined
p_SCMR_combined <- p_SCMR_comparison + p_SCMR_position_significant +
  plot_layout(guides = "collect", axes = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(filename = paste0(results_folder, "06_SCMR_combined.pdf"), plot = p_SCMR_combined, width = 10, height = 6)

## both all together

p_SCMR_combined_all <- p_SCMR_comparison_all + p_SCMR_position_all +
  plot_layout(guides = "collect", axes = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(filename = paste0(results_folder, "06_SCMR_combined_all.pdf"), plot = p_SCMR_combined_all, width = 10, height = 6)



message("SCMR comparison plots created and saved.")

