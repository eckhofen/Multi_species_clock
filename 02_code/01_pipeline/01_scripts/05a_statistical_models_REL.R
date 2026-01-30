#### Overview ####
# multivariate linear model (MLM) for relative age

#### Settings ####
data_folder <- "01_data/"
intermediate_folder <- "01_data/03_intermediate/07_data_comparison/"
dir.create(intermediate_folder, recursive = TRUE, showWarnings = FALSE)

# setting up color palette 
colpal_CB <- c("#c06d00", "#f9cf6e", "#6a5d00", "#44a02b", "#008649", "#12ebf0", "#65a9ff", "#004588", "#660077", "#ff98f7", "#954674", "#630041")
colpal_CB_01 <- colpal_CB[c(TRUE, FALSE)]
colpal_CB_02 <- colpal_CB[c(FALSE, TRUE)]

colpal_CB_a <- c("#f8cbb1","#006786","#182057","#6b6300","#ff8ab9","#f1aaff","#bb005a","#013aa8","#01ef9a","#fa8200","#ee0028","#26c100")
colpal_CB_a_01 <- colpal_CB_a[1:6]
colpal_CB_a_02 <- colpal_CB_a[7:12]

#### Preparation ####
library(tidyverse)
library(tidymodels)

#### Loading data ####
required_input <- "01_data/03_intermediate/06_model_creation/all_meth_values_selected.RData"
if (!file.exists(required_input)) {
  message(
    "Missing required intermediate input: ", required_input, ".\n",
    "This typically means script 04 did not generate model inputs because script 03 could not generate methylation values (private raw inputs unavailable).\n",
    "Skipping relative-age model fitting."
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
# AC_meth_values_selected$rel_age <- AC_meth_values_selected$rel_age * 25
# AS_meth_values_selected$rel_age <- AS_meth_values_selected$rel_age * 54
# EH_meth_values_selected$rel_age <- EH_meth_values_selected$rel_age * 20
# ZF_meth_values_selected$rel_age <- ZF_meth_values_selected$rel_age * 5

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
mlm_rel_age_summary <- summary(mlm_test)
save(file = paste0(intermediate_folder, "mlm_rel_age_summary.Rdata"), mlm_rel_age_summary)
