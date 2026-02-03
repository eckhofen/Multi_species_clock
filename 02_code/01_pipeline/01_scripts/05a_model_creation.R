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