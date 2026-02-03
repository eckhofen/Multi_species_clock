# Metadata ----------------------------------------------------------------
# Project: Piscine Multispecies Epigenetic Clock
# Description: Correlation testing
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026-02

# DISCLAIMER:
# by this time and date the sequencing data is not available yet. Start here to reproduce the results

#### Overview ####
## correlation testing between the selected CpGs and to age as well

#### Settings ####
data_folder <- "01_data/"
save_folder <- "01_data/03_intermediate/05_correlation_data/"
results_folder <- "03_results/01_figures/"
dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

load("01_data/04_metadata/color_palettes.RData")

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
load("01_data/03_intermediate/04_methyl_values/all_meth_values.Rdata")

# shared methylation sites
load("01_data/03_intermediate/04_methyl_values/all_methyl_sites.Rdata")

# age metadata
load("01_data/03_intermediate/04_methyl_values/HS_all_age.Rdata")

#### helper functions ####
# function to run correlation tests for all the CpGs 
cor.test.age <- function(methyl_values, age, SCMR = "not_defined", species = "undefined", method = "pearson") {
  correlation_results <- list()
  print(paste0("Running correlation test against age with ", method, " method. Results are stored in tibble."))
  
  for (i in seq_len(ncol(methyl_values))) {
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
    SCMR = SCMR,
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

# function to choose only the highest correlating CpGs per SMR
select.max.cor <- function(cor_tibble, filter_significant = FALSE) {
  filtered_data <- cor_tibble
  if(filter_significant == TRUE) {filtered_data <- filter(filtered_data, significant)}
  
  filtered_data <- filtered_data %>% 
    group_by(SCMR) %>%
    # add a temporary column for the absolute correlation values
    mutate(abs_correlation = abs(Correlation)) %>%
    # for each group, filter the row with the max absolute correlation
    filter(abs_correlation == max(abs_correlation)) %>%
    # remove the temporary column
    select(-abs_correlation) %>%
    # ensure only one row per group if there are ties
    slice(1)
  return(filtered_data)
}

#### Correlation tests ####

AC_region_id <- AC_methyl_sites$SCMR
AC_cor_age_pearson <- cor.test.age(AC_meth_values, AC_age, AC_region_id, species = "AC")
AC_cor_age_filtered_pearson <- cor.test.age.filter(AC_cor_age_pearson, 0.05)

AS_region_id <- AS_methyl_sites$SCMR
AS_cor_age_pearson <- cor.test.age(AS_meth_values, AS_age, AS_region_id, species = "AS")
AS_cor_age_filtered_pearson <- cor.test.age.filter(AS_cor_age_pearson, 0.05)

EH_region_id <- EH_methyl_sites$SCMR
EH_cor_age_pearson <- cor.test.age(EH_meth_values, EH_age, EH_region_id, species = "EH")
EH_cor_age_filtered_pearson <- cor.test.age.filter(EH_cor_age_pearson, 0.05)

ZF_region_id <- ZF_methyl_sites$SCMR
ZF_cor_age_pearson <- cor.test.age(ZF_meth_values, ZF_age, ZF_region_id, species = "ZF")
ZF_cor_age_filtered_pearson <- cor.test.age.filter(ZF_cor_age_pearson, 0.05)

cor_all <- rbind(AC_cor_age_filtered_pearson,
                 AS_cor_age_filtered_pearson,
                 EH_cor_age_filtered_pearson, 
                 ZF_cor_age_filtered_pearson)

save(cor_all, file = paste0(save_folder, "cor_all.RData"))

## check number of selected CpGs which are significant
ncol(AC_meth_values[AC_cor_age_filtered_pearson$significant] == TRUE) # 154
ncol(AS_meth_values[AS_cor_age_filtered_pearson$significant] == TRUE) # 120
ncol(EH_meth_values[EH_cor_age_filtered_pearson$significant] == TRUE) # 224
ncol(ZF_meth_values[ZF_cor_age_filtered_pearson$significant] == TRUE) # 64

## selecting only positively correlating samples
AC_pos_cor_CpGs <- select.max.cor(AC_cor_age_filtered_pearson[AC_cor_age_filtered_pearson$Correlation > 0,])
AS_pos_cor_CpGs <- select.max.cor(AS_cor_age_filtered_pearson[AS_cor_age_filtered_pearson$Correlation > 0,])
EH_pos_cor_CpGs <- select.max.cor(EH_cor_age_filtered_pearson[EH_cor_age_filtered_pearson$Correlation > 0,])
ZF_pos_cor_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson[ZF_cor_age_filtered_pearson$Correlation > 0,])

all_pos_cor_CpG <- rbind(AC_pos_cor_CpGs, AS_pos_cor_CpGs, EH_pos_cor_CpGs, ZF_pos_cor_CpGs)

# keeping only CpGs when occurring in all species
all_pos_cor_CpG_common <- all_pos_cor_CpG %>%
  group_by(SCMR) %>%
  filter(n() == 4) %>%
  ungroup

## selecting only negatively correlating samples
AC_neg_cor_CpGs <- select.max.cor(AC_cor_age_filtered_pearson[AC_cor_age_filtered_pearson$Correlation < 0,])
AS_neg_cor_CpGs <- select.max.cor(AS_cor_age_filtered_pearson[AS_cor_age_filtered_pearson$Correlation < 0,])
EH_neg_cor_CpGs <- select.max.cor(EH_cor_age_filtered_pearson[EH_cor_age_filtered_pearson$Correlation < 0,])
ZF_neg_cor_CpGs <- select.max.cor(ZF_cor_age_filtered_pearson[ZF_cor_age_filtered_pearson$Correlation < 0,])

all_neg_cor_CpG <- rbind(AC_neg_cor_CpGs, AS_neg_cor_CpGs, EH_neg_cor_CpGs, ZF_neg_cor_CpGs)

# keeping only CpGs when occurring in all species
all_neg_cor_CpG_common <- all_neg_cor_CpG %>%
  group_by(SCMR) %>%
  filter(n() == 4) %>%
  ungroup

## selecting a mixture of both
temp_index_vec <- (all_pos_cor_CpG_common$SCMR %in% all_neg_cor_CpG_common$SCMR) == FALSE
all_mix_cor_CpG_common <- rbind(all_neg_cor_CpG_common, all_pos_cor_CpG_common[temp_index_vec,]) %>%
  group_by(SCMR) %>%
  filter(n() == 4) %>%
  ungroup

save_path <- paste0(save_folder, "all_mix_cor_CpG_common.RData")
save(all_mix_cor_CpG_common, file = save_path)

## selecting all and renaming the SMRs before
all_pos_cor_CpG_common$SCMR_cor <- paste0(all_pos_cor_CpG_common$SCMR, "_pos")
all_neg_cor_CpG_common$SCMR_cor <- paste0(all_neg_cor_CpG_common$SCMR, "_neg")

all_cor_CpG <- rbind(all_pos_cor_CpG_common, all_neg_cor_CpG_common)

# overview of significant values 
overview_CpGs <- rbind(summary(all_pos_cor_CpG_common$P_value), summary(all_neg_cor_CpG_common$P_value), summary(all_pos_cor_CpG_common$Correlation), summary(abs(all_neg_cor_CpG_common$Correlation)))
row.names(overview_CpGs) <- c("pos_cor p-value", "neg_cor p-value", "pos_cor correlation", "neg_cor correlation")
overview_CpGs

#AC
AC_meth_values_selected <- AC_meth_values[,colnames(AC_meth_values) %in% all_mix_cor_CpG_common$Site]
AC_name_index <- match(colnames(AC_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(AC_meth_values_selected) <- all_mix_cor_CpG_common$SCMR[AC_name_index]
AC_meth_values_selected <- AC_meth_values_selected[, order(colnames(AC_meth_values_selected))]
AC_meth_values_selected$rel_age <- AC_age/25
AC_meth_values_selected$age <- AC_age
AC_meth_values_selected$species <- "AC"

#AS
AS_meth_values_selected <- AS_meth_values[,colnames(AS_meth_values) %in% all_mix_cor_CpG_common$Site]
AS_name_index <- match(colnames(AS_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(AS_meth_values_selected) <- all_mix_cor_CpG_common$SCMR[AS_name_index]
AS_meth_values_selected <- AS_meth_values_selected[, order(colnames(AS_meth_values_selected))]
AS_meth_values_selected$rel_age <- AS_age/54
AS_meth_values_selected$age <- AS_age
AS_meth_values_selected$species <- "AS"

#EH
EH_meth_values_selected <- EH_meth_values[,colnames(EH_meth_values) %in% all_mix_cor_CpG_common$Site]
EH_name_index <- match(colnames(EH_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(EH_meth_values_selected) <- all_mix_cor_CpG_common$SCMR[EH_name_index]
EH_meth_values_selected <- EH_meth_values_selected[, order(colnames(EH_meth_values_selected))]
EH_meth_values_selected$rel_age <- EH_age/20
EH_meth_values_selected$age <- EH_age
EH_meth_values_selected$species <- "EH"

#ZF
ZF_meth_values_selected <- ZF_meth_values_imputed[,colnames(ZF_meth_values_imputed) %in% all_mix_cor_CpG_common$Site]
ZF_name_index <- match(colnames(ZF_meth_values_selected), all_mix_cor_CpG_common$Site)
colnames(ZF_meth_values_selected) <- all_mix_cor_CpG_common$SCMR[ZF_name_index]
ZF_meth_values_selected <- ZF_meth_values_selected[, order(colnames(ZF_meth_values_selected))]
ZF_meth_values_selected$rel_age <- ZF_age/5
ZF_meth_values_selected$age <- ZF_age
ZF_meth_values_selected$species <- "ZF"

# combining all data 
all_meth_values_selected <- rbind(AC_meth_values_selected, AS_meth_values_selected, EH_meth_values_selected, ZF_meth_values_selected)

# saving selected values
save_folder <- paste0(data_folder, "03_intermediate/06_model_creation/")
dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)
save_path <- paste0(save_folder, "all_meth_values_selected.RData")
save(AC_meth_values_selected, AS_meth_values_selected, EH_meth_values_selected, ZF_meth_values_selected, all_meth_values_selected, file = save_path)

## extracting all methylvalues for ALL SMRs

#AC
AC_meth_values_selected_all <- AC_meth_values[,colnames(AC_meth_values) %in% all_cor_CpG$Site]
AC_name_index <- match(colnames(AC_meth_values_selected_all), all_cor_CpG$Site)
colnames(AC_meth_values_selected_all) <- all_cor_CpG$SCMR_cor[AC_name_index]
AC_meth_values_selected_all <- AC_meth_values_selected_all[, order(colnames(AC_meth_values_selected_all))]
AC_meth_values_selected_all$rel_age <- AC_age/25
AC_meth_values_selected_all$age <- AC_age
AC_meth_values_selected_all$species <- "AC"

#AS
AS_meth_values_selected_all <- AS_meth_values[,colnames(AS_meth_values) %in% all_cor_CpG$Site]
AS_name_index <- match(colnames(AS_meth_values_selected_all), all_cor_CpG$Site)
colnames(AS_meth_values_selected_all) <- all_cor_CpG$SCMR_cor[AS_name_index]
AS_meth_values_selected_all <- AS_meth_values_selected_all[, order(colnames(AS_meth_values_selected_all))]
AS_meth_values_selected_all$rel_age <- AS_age/54
AS_meth_values_selected_all$age <- AS_age
AS_meth_values_selected_all$species <- "AS"

#EH
EH_meth_values_selected_all <- EH_meth_values[,colnames(EH_meth_values) %in% all_cor_CpG$Site]
EH_name_index <- match(colnames(EH_meth_values_selected_all), all_cor_CpG$Site)
colnames(EH_meth_values_selected_all) <- all_cor_CpG$SCMR_cor[EH_name_index]
EH_meth_values_selected_all <- EH_meth_values_selected_all[, order(colnames(EH_meth_values_selected_all))]
EH_meth_values_selected_all$rel_age <- EH_age/20
EH_meth_values_selected_all$age <- EH_age
EH_meth_values_selected_all$species <- "EH"

#ZF
ZF_meth_values_selected_all <- ZF_meth_values_imputed[,colnames(ZF_meth_values_imputed) %in% all_cor_CpG$Site]
ZF_name_index <- match(colnames(ZF_meth_values_selected_all), all_cor_CpG$Site)
colnames(ZF_meth_values_selected_all) <- all_cor_CpG$SCMR_cor[ZF_name_index]
ZF_meth_values_selected_all <- ZF_meth_values_selected_all[, order(colnames(ZF_meth_values_selected_all))]
ZF_meth_values_selected_all$rel_age <- ZF_age/5
ZF_meth_values_selected_all$age <- ZF_age
ZF_meth_values_selected_all$species <- "ZF"

# combining all data 
all_meth_values_all_SMR <- rbind(AC_meth_values_selected_all, AS_meth_values_selected_all, EH_meth_values_selected_all, ZF_meth_values_selected_all)

# saving selected values
save_folder <- paste0(data_folder, "03_intermediate/06_model_creation/")
dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)
save_path <- paste0(save_folder, "all_meth_values_all_SMR.RData")
save(AC_meth_values_selected_all, AS_meth_values_selected_all, EH_meth_values_selected_all, ZF_meth_values_selected_all, all_meth_values_all_SMR, file = save_path)

#### plotting ####
## all
ggplot(cor_all, aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species, alpha = significant)) +
  scale_color_manual(values = color_species) +
  # facet_row(~SMR) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggplot(cor_all, aes()) +
  geom_point(aes(x = Site, y = Correlation, color = species, alpha = significant)) +
  # geom_line(aes(x = c(-1,1), y = log2(0.05), color = "#CC79A7")) +
  scale_color_manual(values = color_species) +
  facet_wrap(~SCMR) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## only significant
ggplot(subset(cor_all, significant == TRUE), aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species)) +
  # facet_wrap(~SMR) +
  scale_color_manual(values = color_species) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## only selected
if (exists("all_sig_CpGs_common")) {
  ggplot(all_sig_CpGs_common, aes()) +
    geom_point(aes(y = Correlation, x = Site, color = species)) +
    facet_row(~SCMR) +
    scale_color_manual(values = color_species) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

## only cor positive
ggplot(all_pos_cor_CpG, aes()) +
  geom_point(aes(y = Correlation, x = Site, color = species, alpha = significant)) +
  facet_wrap(~SCMR) +
  scale_color_manual(values = color_species) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# plotting SMR
ggplot(selected_methyl_values, aes(x = Site, y = Methylation_Value)) +
  geom_boxplot(aes(group = Site_f, fill = species, color = species), alpha = 0.9, outlier.size = 0.1) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_fill_manual(values = color_species) +
  scale_color_manual(values = color_species) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (ZF), European Hake (EH), Zebrafish (ZF) (human rgenome)",
       subtitle = "Selected values are correlating with age")


ggplot(selected_methyl_values, aes(x = species, y = Methylation_Value)) +
  ggforce::geom_sina(aes(color = rel_age, shape = species)) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_color_manual(aesthetics = "legend") +
  theme_classic() +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (ZF), European Hake (EH), Zebrafish (ZF) (human rgenome)")

