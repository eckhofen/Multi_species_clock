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

# color palette
load(paste0(data_folder, "04_metadata/color_palettes.RData"))

#### Preparation ####
library(tidyverse)
library(tidymodels)

#### Loading data ####
load("01_data/03_intermediate/06_model_creation/all_meth_values_selected.RData")

#### Data preparation ####
message("Data preparation...")
# split breaks
ds_breaks <- 3
# training/testing set proportion
ds_prop <- 5/6

# stratified split based on age
set.seed(123)
AC_split <- initial_split(AC_meth_values_selected, strata = "age", breaks = ds_breaks, prop = ds_prop)

set.seed(123)
AS_split <- initial_split(AS_meth_values_selected, strata = "age", breaks = ds_breaks, prop = ds_prop)

set.seed(123)
EH_split <- initial_split(EH_meth_values_selected, strata = "age", breaks = ds_breaks, prop = ds_prop)

set.seed(123)
ZF_split <- initial_split(ZF_meth_values_selected, strata = "age", breaks = ds_breaks, prop = ds_prop)

# combining data into training and testing sets
meth_train <- rbind(training(AC_split), training(AS_split), training(EH_split), training(ZF_split))
meth_test <- rbind(testing(AC_split), testing(AS_split), testing(EH_split), testing(ZF_split))

metadata_train <- meth_train %>% select(age, rel_age, species)
metadata_test <- meth_test %>% select(age, rel_age, species)

data_train <- meth_train %>% select(-age, -rel_age, -species)
data_test <- meth_test %>% select(-age, -rel_age, -species)

# saving
save(file = paste0(intermediate_folder, "all_data.Rdata"), data_train, data_test, metadata_train, metadata_test)

# model creation
message("Model creation...")
mlm_chron <- lm(metadata_train$age ~ ., data = data_train)
mlm_chron_summary <- summary(mlm_chron)

save(file = paste0(intermediate_folder, "mlm_chron_summary.Rdata"), mlm_chron_summary)

mlm_rel <- lm(metadata_train$rel_age ~ ., data = data_train)
mlm_rel_summary <- summary(mlm_rel)

save(file = paste0(intermediate_folder, "mlm_rel_summary.Rdata"), mlm_rel_summary)

# log models
mlm_chron_log <- lm(log(metadata_train$age) ~ ., data = data_train)
mlm_rel_log <- lm(-log(-log(metadata_train$rel_age)) ~ ., data = data_train)

mlm_chron_log_summary <- summary(mlm_chron_log)
mlm_rel_log_summary <- summary(mlm_rel_log)

# save models
save(file = paste0(intermediate_folder, "models.Rdata"), mlm_chron, mlm_rel, mlm_chron_log, mlm_rel_log)

# key parameters
message(paste("MLM (chronological age), R2: ", round(mlm_chron_summary$r.squared, 3)))
message(paste("MLM (chronological age), R: ", round(sqrt(mlm_chron_summary$r.squared), 3)))

message(paste("MLM (relative age), R2: ", round(mlm_rel_summary$r.squared, 3)))
message(paste("MLM (relative age), R: ", round(sqrt(mlm_rel_summary$r.squared), 3)))

message(paste("MLM (chronological age, log), R2: ", round(mlm_chron_log_summary$r.squared, 3)))
message(paste("MLM (chronological age, log), R: ", round(sqrt(mlm_chron_log_summary$r.squared), 3)))

message(paste("MLM (relative age, log), R2: ", round(mlm_rel_log_summary$r.squared, 3)))
message(paste("MLM (relative age, log), R: ", round(sqrt(mlm_rel_log_summary$r.squared), 3)))

message("Models created and saved.")

# get key parameters
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
  SCMR = gsub("^SMR_", "SCMR_", rownames(mlm_age_sign_coef))
)

df_mlm_age_all <- data.frame(
  model = "Chronological age",
  reg_coef = mlm_age_coef[, 1],
  reg_coef_scaled = mlm_age_coef[, 1] / max_abs_age_all,
  p_value = mlm_age_coef[, 4],
  SCMR = gsub("^SMR_", "SCMR_", rownames(mlm_age_coef))
)

# Extract coefficients for relative age model
mlm_rel_coef <- mlm_rel_summary$coefficients[-1, , drop = FALSE]  # Remove intercept
mlm_rel_age_sign_coef <- mlm_rel_coef[mlm_rel_coef[, 4] < significance, c(1, 4), drop = FALSE]

max_abs_rel_age_sign <- if (nrow(mlm_rel_age_sign_coef) == 0) NA_real_ else max(abs(mlm_rel_age_sign_coef[, 1]), na.rm = TRUE)
max_abs_rel_age_all <- max(abs(mlm_rel_coef[, 1]), na.rm = TRUE)

df_mlm_rel_age_sign_coef <- data.frame(
  model = "Relative age",
  reg_coef = mlm_rel_age_sign_coef[, 1],
  reg_coef_scaled = ifelse(is.na(max_abs_rel_age_sign) || max_abs_rel_age_sign == 0, NA_real_, mlm_rel_age_sign_coef[, 1] / max_abs_rel_age_sign),
  p_value = mlm_rel_age_sign_coef[, 2],
  SCMR = gsub("^SMR_", "SCMR_", rownames(mlm_rel_age_sign_coef))
)

df_mlm_rel_age_all <- data.frame(
  model = "Relative age",
  reg_coef = mlm_rel_coef[, 1],
  reg_coef_scaled = mlm_rel_coef[, 1] / max_abs_rel_age_all,
  p_value = mlm_rel_coef[, 4],
  SCMR = gsub("^SMR_", "SCMR_", rownames(mlm_rel_coef))
)

# Combine significant SMRs
df_SMR_comparison <- rbind(df_mlm_age_sign_coef, df_mlm_rel_age_sign_coef)
df_SMR_comparison <- df_SMR_comparison %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ "'",
    TRUE ~ ""
  ))

# Combine all SMRs
df_SMR_comparison_all <- rbind(df_mlm_age_all, df_mlm_rel_age_all)
df_SMR_comparison_all <- df_SMR_comparison_all %>%
  mutate(significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    p_value < 0.1 ~ "'",
    TRUE ~ ""
  ))

# Save key parameters
save(file = paste0(intermediate_folder, "df_SMR_comparison.Rdata"), df_SMR_comparison)
save(file = paste0(intermediate_folder, "df_SMR_comparison_all.Rdata"), df_SMR_comparison_all)

message("Key parameters extracted and saved.")

#### SMR Comparison Plots ####
message("Creating SMR comparison plots...")

results_folder <- "03_results/01_figures/"
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

extension <- ".pdf"

# Check for color_compare palette, create default if not available
if (!exists("color_compare")) {
  color_compare <- c("#005AB5", "#DC3220")
}

### SMR comparison plots
## Only significant SMRs
plot_SMR_comparison <- ggplot(df_SMR_comparison, aes(x = reg_coef_scaled, y = SCMR, fill = model)) +
  geom_point(aes(color = model), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model), position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", y = "Shared methylation region", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.15, color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.8) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare) +
  theme_bw()

ggsave(file = paste0(results_folder, "05a_SCMR_comparison_significant", extension), plot_SMR_comparison, width = 5, height = 5)

## All SMRs
plot_SMR_comparison_all <- ggplot(df_SMR_comparison_all, aes(x = reg_coef_scaled, y = SCMR, fill = model)) +
  geom_point(aes(color = model), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model), position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", y = "Shared methylation region", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.1, color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.75) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare) +
  theme_bw()

ggsave(file = paste0(results_folder, "05a_SCMR_comparison_all", extension), plot_SMR_comparison_all, width = 10, height = 7)

message("SMR comparison plots created and saved.")
message("05a_model_creation.R completed successfully.")
