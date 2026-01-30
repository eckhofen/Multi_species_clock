#### Overview ####
# multivariate linear model (MLM) for chronological age

#### Settings ####
data_folder <- "01_data/"
intermediate_folder <- "01_data/03_intermediate/07_data_comparison/"
dir.create(intermediate_folder, recursive = TRUE, showWarnings = FALSE)

#### Preparation ####
library(tidyverse)
library(tidymodels)

#### Loading data ####
required_input <- "01_data/03_intermediate/06_model_creation/all_meth_values_selected.RData"
if (!file.exists(required_input)) {
  message(
    "Missing required intermediate input: ", required_input, ".\n",
    "This typically means script 04 did not generate model inputs because script 03 could not generate methylation values (private raw inputs unavailable).\n",
    "Skipping chronological-age model fitting."
  )
  quit(status = 0)
}

load(required_input)

#### modifying age (optional) ####
# AC_meth_values_selected$rel_age <- AC_meth_values_selected$rel_age / 1.4
# AS_meth_values_selected$rel_age
# EH_meth_values_selected$rel_age <- EH_meth_values_selected$rel_age / 1
# ZF_meth_values_selected$rel_age

#### testing chronological age
AC_meth_values_selected$rel_age <- AC_meth_values_selected$rel_age * 25
AS_meth_values_selected$rel_age <- AS_meth_values_selected$rel_age * 54
EH_meth_values_selected$rel_age <- EH_meth_values_selected$rel_age * 20
ZF_meth_values_selected$rel_age <- ZF_meth_values_selected$rel_age * 5

all_meth_values_selected <- rbind(AC_meth_values_selected,AS_meth_values_selected,EH_meth_values_selected,ZF_meth_values_selected)
#### Data splitting ####
# defining arguments
ds_breaks <- 3
ds_prop <- 5/6
# using a stratified splitting technique
set.seed(123)
AC_split <- initial_split(AC_meth_values_selected, strata = "rel_age", breaks = ds_breaks, prop = ds_prop)

set.seed(123)
AS_split <- initial_split(AS_meth_values_selected, strata = "rel_age", breaks = ds_breaks, prop = ds_prop)

set.seed(123)
EH_split <- initial_split(EH_meth_values_selected, strata = "rel_age", breaks = ds_breaks, prop = ds_prop)

set.seed(123)
ZF_split <- initial_split(ZF_meth_values_selected, strata = "rel_age", breaks = ds_breaks, prop = ds_prop)

# combining data into training and testing sets
meth_train <- rbind(training(AC_split), training(AS_split), training(EH_split), training(ZF_split))
meth_test <- rbind(testing(AC_split), testing(AS_split), testing(EH_split), testing(ZF_split))

### defining data
# training data
X <- meth_train %>% select(-rel_age, -species)

Y <- meth_train[,"rel_age"]

# testing data
X_test <- meth_test %>% select(-rel_age, -species)

Y_test <- meth_test[,"rel_age"]

#### Multivariate linear model (MLM) ####
mlm_test <- lm(Y ~ ., data = X)
mlm_age_summary <- summary(mlm_test)
save(file = paste0(intermediate_folder, "mlm_age_summary.Rdata"), mlm_age_summary)
