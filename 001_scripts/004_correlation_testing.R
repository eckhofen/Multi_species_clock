#### Overview ####
## correlation testing between the selected CpGs and to age as well

#### Settings ####
setwd("/powerplant/workspace/cfngle/script_GH/Multi_species_clock/")
data_folder <- "/powerplant/workspace/cfngle/script_GH/Multi_species_clock/000_data/"
save_folder <- paste0(data_folder, "005_correlation_data/") # folder where extracted sequences will be saved

#### Preparation ####
library(tidyverse)
library(tidymodels)
library(tibble)
library(dplyr)
library(ggplot2)
library(patchwork)
library(caret)

## loading data 
# methyl values
load("000_data/004_methyl_values/HS_AC_meth_values.Rdata")
load("000_data/004_methyl_values/HS_AS_meth_values.Rdata")
load("000_data/004_methyl_values/HS_EH_meth_values.Rdata")
load("000_data/004_methyl_values/HS_ZF_meth_values_imputed.Rdata")
# shared methylation sites
load("000_data/004_methyl_values/HS_AC_methyl_sites.Rdata")
load("000_data/004_methyl_values/HS_AS_methyl_sites.Rdata")
load("000_data/004_methyl_values/HS_EH_methyl_sites.Rdata")
load("000_data/004_methyl_values/HS_ZF_methyl_sites.Rdata")
# age metadata
load("000_data/004_methyl_values/HS_all_age.Rdata")

#### functions ####
# function to run correlation tests for all the CpGs 
cor.test.age <- function(methyl_values, age, SMR = "not_defined", species = "undefined", method = "pearson") {
  correlation_results <- list()
  print(paste0("Running correlation test against age with ", method, " method. Results are stored in tibble."))
  
  for (i in 1:ncol(methyl_values)) {
    site_name <- colnames(methyl_values)[i]
    # perform correlation test with age
    test_result <- cor.test(methyl_values[,i], age, method = method)
    
    # store the results
    correlation_results[[site_name]] <- list(
      correlation_coefficient = test_result$estimate,
      p_value = test_result$p.value
    )
  }
  
  # convert the results list to a more convenient format
  correlation_summary <- tibble(
    Site = names(correlation_results),
    Correlation = sapply(correlation_results, function(x) x$correlation_coefficient),
    P_value = sapply(correlation_results, function(x) x$p_value),
    SMR = SMR,
    species = species
  )
  return(correlation_summary)
}

# checks whether the correlation is significant or not (p-value should be defined)
cor.test.age.filter <- function(input, p_value = 0.05) {
  significant_vector <- as.vector(ifelse(input$P_value <= p_value, TRUE, FALSE))
  input$significant <- significant_vector
  return(input)
}

#### Correlation tests ####

AC_cor_age_pearson <- cor.test.age(AC_meth_values, AC_age, AC_methyl_sites$SMR, species = "AC")
AC_cor_age_filtered_pearson <- cor.test.age.filter(AC_cor_age_pearson, 0.05)

AS_cor_age_pearson <- cor.test.age(AS_meth_values, AS_age, AS_methyl_sites$SMR, species = "AS")
AS_cor_age_filtered_pearson <- cor.test.age.filter(AS_cor_age_pearson, 0.05)

EH_cor_age_pearson <- cor.test.age(EH_meth_values, EH_age, EH_methyl_sites$SMR, species = "EH")
EH_cor_age_filtered_pearson <- cor.test.age.filter(EH_cor_age_pearson, 0.05)

ZF_cor_age_pearson <- cor.test.age(ZF_meth_values_imputed, ZF_age, ZF_methyl_sites$SMR, species = "ZF")
ZF_cor_age_filtered_pearson <- cor.test.age.filter(ZF_cor_age_pearson, 0.05)

cor_all <- rbind(AC_cor_age_filtered_pearson,
                 AS_cor_age_filtered_pearson,
                 EH_cor_age_filtered_pearson, 
                 ZF_cor_age_filtered_pearson)

## check number of selected CpGs which are significant
ncol(AC_meth_values[AC_cor_age_filtered_pearson$significant] == TRUE)
ncol(AS_meth_values[AS_cor_age_filtered_pearson$significant] == TRUE)
ncol(EH_meth_values[EH_cor_age_filtered_pearson$significant] == TRUE)
ncol(ZF_meth_values_imputed[ZF_cor_age_filtered_pearson$significant] == TRUE)
### 154, 120, 224, 64

#### function to choose only the highest correlating CpGs per SMR ####

select.max.cor <- function(cor_tibble, filter_significant = FALSE) {
  filtered_data <- cor_tibble
  if(filter_significant == TRUE) {filtered_data <- filter(filtered_data, significant)}
  
  filtered_data <- filtered_data %>% 
    group_by(SMR) %>%
    # Add a temporary column for the absolute correlation values
    mutate(abs_correlation = abs(Correlation)) %>%
    # For each group, filter the row with the max absolute correlation
    filter(abs_correlation == max(abs_correlation)) %>%
    # Remove the temporary column
    select(-abs_correlation) %>%
    # Optionally, ensure only one row per group if there are ties
    slice(1)
  return(filtered_data)
}

## selecting only positively correlating samples
AC_pos_cor_CpGs <- select.max.cor(AC_cor_age_filtered_pearson[AC_cor_age_filtered_pearson$Correlation > 0,])
AS_pos_cor_CpGs <- select.max.cor(AS_cor_age_filtered_pearson[AS_cor_age_filtered_pearson$Correlation > 0,])
EH_pos_cor_CpGs <- select.max.cor(EH_cor_age_filtered_pearson[EH_cor_age_filtered_pearson$Correlation > 0,])
ZF_pos_cor_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson[ZF_cor_age_filtered_pearson$Correlation > 0,])

all_pos_cor_CpG <- rbind(AC_pos_cor_CpGs, AS_pos_cor_CpGs, EH_pos_cor_CpGs, ZF_pos_cor_CpGs)

# keeping only CpGs when occurring in all species 
all_pos_cor_CpG_common  <- all_pos_cor_CpG %>% 
  group_by(SMR) %>% 
  filter(n() == 4) %>% 
  ungroup

## selecting only negatively correlating samples
AC_neg_cor_CpGs <- select.max.cor(AC_cor_age_filtered_pearson[AC_cor_age_filtered_pearson$Correlation < 0,])
AS_neg_cor_CpGs <- select.max.cor(AS_cor_age_filtered_pearson[AS_cor_age_filtered_pearson$Correlation < 0,])
EH_neg_cor_CpGs <- select.max.cor(EH_cor_age_filtered_pearson[EH_cor_age_filtered_pearson$Correlation < 0,])
ZF_neg_cor_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson[ZF_cor_age_filtered_pearson$Correlation < 0,])

all_neg_cor_CpG <- rbind(AC_neg_cor_CpGs, AS_neg_cor_CpGs, EH_neg_cor_CpGs, ZF_neg_cor_CpGs)

# keeping only CpGs when occurring in all species 
all_neg_cor_CpG_common  <- all_neg_cor_CpG %>% 
  group_by(SMR) %>% 
  filter(n() == 4) %>% 
  ungroup

## selecting a mixture of both
temp_index_vec <- (all_pos_cor_CpG_common$SMR %in% all_neg_cor_CpG_common$SMR) == FALSE
all_mix_cor_CpG_common  <- rbind(all_neg_cor_CpG_common, all_pos_cor_CpG_common[temp_index_vec,]) %>% 
  group_by(SMR) %>% 
  filter(n() == 4) %>% 
  ungroup

save_path <- paste0(save_folder, "all_mix_cor_CpG_common.RData")
save(all_mix_cor_CpG_common, file = save_path)

## selecting significant ones
# AC_sig_CpGs <- select.max.cor(AC_cor_age_filtered_pearson, TRUE)
# AS_sig_CpGs <- select.max.cor(AS_cor_age_filtered_pearson, TRUE)
# EH_sig_CpGs <- select.max.cor(EH_cor_age_filtered_pearson, TRUE)
# ZF_sig_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson, TRUE)
# 
# all_sig_CpGs <- rbind(AC_sig_CpGs, AS_sig_CpGs, EH_sig_CpGs, ZF_sig_CpGs)
# 
# # only keeping SMR which are present in all species (n = 4)
# all_sig_CpGs_common <- all_sig_CpGs %>% 
#   group_by(SMR) %>% 
#   filter(n() == 4) %>% 
#   ungroup

#### Retrieving methylation values ####
# selecting a mixture of pos and neg correlating CpGs

#AC
AC_meth_values_selected <- AC_meth_values[,colnames(AC_meth_values) %in% all_mix_cor_CpG_common$Site]
AC_name_index <- match(colnames(AC_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(AC_meth_values_selected) <- all_mix_cor_CpG_common$SMR[AC_name_index]
AC_meth_values_selected <- AC_meth_values_selected[, order(colnames(AC_meth_values_selected))]
AC_meth_values_selected$rel_age <- AC_age/25
AC_meth_values_selected$species <- "AC"

#AS
AS_meth_values_selected <- AS_meth_values[,colnames(AS_meth_values) %in% all_mix_cor_CpG_common$Site]
AS_name_index <- match(colnames(AS_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(AS_meth_values_selected) <- all_mix_cor_CpG_common$SMR[AS_name_index]
AS_meth_values_selected <- AS_meth_values_selected[, order(colnames(AS_meth_values_selected))]
AS_meth_values_selected$rel_age <- AS_age/54
AS_meth_values_selected$species <- "AS"

#EH
EH_meth_values_selected <- EH_meth_values[,colnames(EH_meth_values) %in% all_mix_cor_CpG_common$Site]
EH_name_index <- match(colnames(EH_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(EH_meth_values_selected) <- all_mix_cor_CpG_common$SMR[EH_name_index]
EH_meth_values_selected <- EH_meth_values_selected[, order(colnames(EH_meth_values_selected))]
EH_meth_values_selected$rel_age <- EH_age/20
EH_meth_values_selected$species <- "EH"

#ZF
ZF_meth_values_selected <- ZF_meth_values_imputed[,colnames(ZF_meth_values_imputed) %in% all_mix_cor_CpG_common$Site]
ZF_name_index <- match(colnames(ZF_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(ZF_meth_values_selected) <- all_mix_cor_CpG_common$SMR[ZF_name_index]
ZF_meth_values_selected <- ZF_meth_values_selected[, order(colnames(ZF_meth_values_selected))]
ZF_meth_values_selected$rel_age <- ZF_age/5
ZF_meth_values_selected$species <- "ZF"

# combining all data 
all_meth_values_selected <- rbind(AC_meth_values_selected, AS_meth_values_selected, EH_meth_values_selected, ZF_meth_values_selected)

# saving selected values
save_folder <- paste0(data_folder, "006_model_creation/")
save_path <- paste0(save_folder, "all_meth_values_selected.RData")
save(all_meth_values_selected, file = save_path)
