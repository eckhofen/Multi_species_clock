#### Overview ####
# Data comparison

#### Settings ####
# change working directory accordingly extension
# setwd("/powerplant/workspace/cfngle/script_GH/Multi_species_clock/")
extension <- ".pdf"
data_folder <- "01_data/"
results_folder <- "03_results/01_figures/"
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

load("01_data/04_metadata/color_palettes.RData")

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
load("01_data/03_intermediate/007_data_comparison/mlm_age_summary.Rdata")
load("01_data/03_intermediate/007_data_comparison/mlm_rel_age_summary.Rdata")
load("01_data/03_intermediate/007_data_comparison/methyl_sites_combined_nor.Rdata")
load("01_data/03_intermediate/005_correlation_data/all_mix_cor_CpG_common.RData")

if (exists("methyl_sites_combined_nor") && ("SMR" %in% colnames(methyl_sites_combined_nor)) && !("SCMR" %in% colnames(methyl_sites_combined_nor))) {
  methyl_sites_combined_nor <- dplyr::rename(methyl_sites_combined_nor, SCMR = SMR)
}


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
                                   SCMR = rownames(mlm_age_sign_coef))
df_mlm_rel_age_sign_coef <- data.frame(model = "Relative age", 
                                       reg_coef = mlm_rel_age_sign_coef[,1], 
                                       reg_coef_scaled = mlm_rel_age_sign_coef[,1]/max(abs(mlm_rel_age_sign_coef[,1])),
                                       p_value = mlm_rel_age_sign_coef[,2], 
                                       SCMR = rownames(mlm_rel_age_sign_coef))

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
                         SCMR = rownames(mlm_age_summary$coefficients[-1,]))
df_mlm_rel_age <- data.frame(model = "Relative age", 
                             reg_coef = mlm_rel_age_summary$coefficients[-1,1],
                             reg_coef_scaled = mlm_rel_age_summary$coefficients[-1,1]/max(abs(mlm_rel_age_summary$coefficients[-1,1])), 
                             p_value = mlm_rel_age_summary$coefficients[-1,4], 
                             SCMR = rownames(mlm_rel_age_summary$coefficients[-1,]))

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

#### AE max life changes ####
## load AE for LF modification
load("01_data/03_intermediate/007_data_comparison/df_AE_t_max_LS_AC_p_30")
df_AE_t_AC_p_30 <- df_AE_t
load("01_data/03_intermediate/007_data_comparison/df_AE_max_LS_AC_p_30")
df_AE_AC_p_30 <- df_AE

load("01_data/03_intermediate/007_data_comparison/df_AE_REL")
df_AE_REL <- df_AE
load("01_data/03_intermediate/007_data_comparison/df_AE_t_REL")
df_AE_t_REL <- df_AE_t

## Students t test

t.test(df_AE_REL$ENR, df_AE_AC_p_30$ENR)

#### Plots ####

### SMR comparison
## only significant
plot_SMR_comparison <- ggplot(df_SMR_comparison, aes(x = reg_coef_scaled, y = SCMR, fill = model)) +
  geom_point(aes(color = model), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model), position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", y = "Shared methylation region", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.15, , color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.8) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare) +
  theme_bw()
plot_SMR_comparison

ggsave(file = paste0(results_folder, "06_SCMR_comparison",extension), plot_SMR_comparison, width = 5, height = 5)

## all
plot_SMR_comparison_all <- ggplot(df_SMR_comparison_all, aes(x = reg_coef_scaled, y = SCMR, fill = model)) +
  geom_point(aes(color = model), position = position_dodge(width = 0.9)) +
  geom_linerange(aes(xmin = 0, xmax = reg_coef_scaled, color = model), position = position_dodge(width = 0.9)) +
  labs(x = "Normalized regression coefficient", y = "Shared methylation region", fill = "Model", color = "Model") +
  geom_text(aes(label = significance, x = reg_coef_scaled * 1.1, , color = model), show.legend = FALSE, position = position_dodge(0.9), vjust = 0.75) +
  scale_fill_manual(values = color_compare) +
  scale_color_manual(values = color_compare) +
  theme_bw()
plot_SMR_comparison_all

ggsave(file = paste0(results_folder, "06_SCMR_comparison_all",extension), plot_SMR_comparison_all, width = 10, height = 7)

### CpG location for SMRs
## CpG position for all 
plot_SMR_position_all <- ggplot(methyl_sites_combined_nor) +
  geom_segment(aes(x = pos_nor_kb, xend = pos_nor_kb+.02, y = SCMR, color = species), linewidth = 2) +
  labs(x = "CpG position (kb)", y = "Shared methylation region", color = "Species") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  # scale_y_discrete(limits = rev(levels(methyl_sites_combined_nor$species))) +
  scale_color_manual(values = color_species) + 
  theme_bw()
# theme(axis.text.x = element_text(angle = -20, hjust = -.05))
plot_SMR_position_all

ggsave(filename = paste0(results_folder, "06_SCMR_position_all.pdf"), plot = plot_SMR_position_all, width = 10, height = 7)

## only significant
# subsetting for significant SMRs only
methyl_sites_significant <- subset(methyl_sites_combined_nor, SCMR %in% df_SMR_comparison$SCMR)


plot_SMR_position_significant <- ggplot(methyl_sites_significant) +
  geom_segment(aes(x = pos_nor_kb, xend = pos_nor_kb+.03, y = SCMR, color = species), linewidth = 5) +
  labs(x = "CpG position (kb)", y = "Shared methylation region", color = "Species") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
  # scale_y_discrete(limits = rev(levels(methyl_sites_combined_nor$species))) +
  scale_color_manual(values = color_species) + 
  theme_bw()
# theme(axis.text.x = element_text(angle = -20, hjust = -.05))
plot_SMR_position_significant

ggsave(filename = paste0(results_folder, "06_SCMR_position_significant.pdf"), plot = plot_SMR_position_significant, width = 6, height = 5)

## both together

plot_SMR_combined <- plot_SMR_comparison + plot_SMR_position_significant +
  plot_layout(guides = "collect", axes = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(filename = paste0(results_folder, "06_SCMR_combined.pdf"), plot = plot_SMR_combined, width = 10, height = 6)

## both all together

plot_SMR_combined_all <- plot_SMR_comparison_all + plot_SMR_position_all +
  plot_layout(guides = "collect", axes = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(filename = paste0(results_folder, "06_SCMR_combined_all.pdf"), plot = plot_SMR_combined_all, width = 10, height = 6)

#### Plot LOSO results ####

#### log-transformed
## load and rename data
factor_train <- factor("Training", levels = c("Training", "Testing"))
factor_test <- factor("Testing", levels = c("Training", "Testing"))
# AC
load("01_data/03_intermediate/007_data_comparison/no_AC_SVM_eps_val.Rdata")
AC_LOSO_plot_train_t <- SVM_eps_eval_chron_t$plot_train + theme(plot.title = element_blank())
AC_LOSO_plot_test_t <- SVM_eps_eval_chron_t$plot_test + theme(plot.title = element_blank())
AC_LOSO_AE <- rbind(data.frame(AE = SVM_eps_eval_chron_t$values_AE_test, type = factor_test),
                    data.frame(AE = SVM_eps_eval_chron_t$values_AE_train, type = factor_train))
# AS
load("01_data/03_intermediate/007_data_comparison/no_AS_SVM_eps_val.Rdata")
AS_LOSO_plot_train_t <- SVM_eps_eval_chron_t$plot_train + theme(plot.title = element_blank())
AS_LOSO_plot_test_t <- SVM_eps_eval_chron_t$plot_test + theme(plot.title = element_blank())
AS_LOSO_AE <- rbind(data.frame(AE = SVM_eps_eval_chron_t$values_AE_test, type = factor_test),
                    data.frame(AE = SVM_eps_eval_chron_t$values_AE_train, type = factor_train))
# EH
load("01_data/03_intermediate/007_data_comparison/no_EH_SVM_eps_val.Rdata")
EH_LOSO_plot_train_t <- SVM_eps_eval_chron_t$plot_train+ theme(plot.title = element_blank())
EH_LOSO_plot_test_t <- SVM_eps_eval_chron_t$plot_test + theme(plot.title = element_blank())
EH_LOSO_AE <- rbind(data.frame(AE = SVM_eps_eval_chron_t$values_AE_test, type = factor_test),
                    data.frame(AE = SVM_eps_eval_chron_t$values_AE_train, type = factor_train))
# ZF
load("01_data/03_intermediate/007_data_comparison/no_ZF_SVM_eps_val.Rdata")
ZF_LOSO_plot_train_t <- SVM_eps_eval_chron_t$plot_train + theme(plot.title = element_blank())
ZF_LOSO_plot_test_t <- SVM_eps_eval_chron_t$plot_test + theme(plot.title = element_blank())
ZF_LOSO_AE <- rbind(data.frame(AE = SVM_eps_eval_chron_t$values_AE_test, type = factor_test),
                    data.frame(AE = SVM_eps_eval_chron_t$values_AE_train, type = factor_train))

### plot AE for each LOSO

# AC
AC_plot_AE_SVM_t <- ggplot(AC_LOSO_AE, aes(x = type, y = AE, pattern = type)) +
  geom_boxplot_pattern(position = position_dodge(width = .9), outlier.shape = NA, pattern_fill = "transparent", pattern_color = "gray10") + 
  scale_pattern_manual(values = c("stripe", "circle")) + 
  labs(y = "Absolute error (years)", x = "Data", pattern = "Data") +
  theme_bw() + 
  theme(legend.position = "") +
  ylim(c(0,2.5)) 

# AS
AS_plot_AE_SVM_t <- ggplot(AS_LOSO_AE, aes(x = type, y = AE, pattern = type)) +
  geom_boxplot_pattern(position = position_dodge(width = .9), outlier.shape = NA, pattern_fill = "transparent", pattern_color = "gray10") + 
  scale_pattern_manual(values = c("stripe", "circle")) + 
  labs(y = "Absolute error (years)", x = "Data", pattern = "Data") +
  theme_bw() + 
  theme(legend.position = "") +
  ylim(c(0,4.5)) 

# EH
EH_plot_AE_SVM_t <- ggplot(EH_LOSO_AE, aes(x = type, y = AE, pattern = type)) +
  geom_boxplot_pattern(position = position_dodge(width = .9), outlier.shape = NA, pattern_fill = "transparent", pattern_color = "gray10") + 
  scale_pattern_manual(values = c("stripe", "circle")) + 
  labs(y = "Absolute error (years)", x = "Data", pattern = "Data") +
  theme_bw() + 
  theme(legend.position = "") +
  ylim(c(0,25)) 

# ZF
ZF_plot_AE_SVM_t <- ggplot(ZF_LOSO_AE, aes(x = type, y = AE, pattern = type)) +
  geom_boxplot_pattern(position = position_dodge(width = .9), outlier.shape = NA, pattern_fill = "transparent", pattern_color = "gray10") + 
  scale_pattern_manual(values = c("stripe", "circle")) + 
  labs(y = "Absolute error (years)", x = "Data", pattern = "Data") +
  theme_bw() + 
  theme(legend.position = "") +
  ylim(c(0,2)) 

plot_LOSO_all_SVM_log <- AC_LOSO_plot_train_t + AC_LOSO_plot_test_t + AC_plot_AE_SVM_t +
  AS_LOSO_plot_train_t + AS_LOSO_plot_test_t + AS_plot_AE_SVM_t +
  EH_LOSO_plot_train_t + EH_LOSO_plot_test_t + EH_plot_AE_SVM_t +
  ZF_LOSO_plot_train_t + ZF_LOSO_plot_test_t + ZF_plot_AE_SVM_t +
  plot_layout(guides = "collect", axes = "collect", ncol = 3, widths = c(3,3,1)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))
plot_LOSO_all_SVM_log

ggsave(filename = paste0(results_folder, "06_LOSO_SVM_all.pdf"), plot = plot_LOSO_all_SVM_log, width = 9, height = 11.5)


#### plotting SMR methylation ####

load("01_data/03_intermediate/004_methyl_values/all_meth_values_long.RData")
load("01_data/03_intermediate/005_correlation_data/all_mix_cor_CpG_common.RData")

all_meth_values_long_sel <- all_meth_values_long[all_meth_values_long$Site %in% all_mix_cor_CpG_common$Site,]

ggplot(AS_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[1], high = colpal_CB_c[5], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values AS (human rgenome)", color = "Age")

ggplot(AC_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[2], high = colpal_CB_c[6], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values AC (human rgenome)", color = "Age")

ggplot(EH_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[3], high = colpal_CB_c[7], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values EH (human rgenome)", color = "Age")

ggplot(ZF_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_sina(aes(color = age), alpha = 0.7) +
  geom_boxplot(aes(group = Site_f), alpha = 0.5) +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_color_gradient(low = colpal_CB_c[4], high = colpal_CB_c[8], guide = "legend") +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values ZF (human rgenome)", color = "Age")

## all boxplot
ggplot(all_meth_values_long, aes(x = Site, y = Methylation_Value)) +
  geom_boxplot(aes(group = Site_f, fill = species), color = NA, alpha = 0.9, outlier.shape = "") +
  facet_wrap(~SCMR, scale = "free_x") +
  scale_fill_manual(values = color_species) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (AS), European Hake (EH), Zebrafish (ZF) (human rgenome)")

ggplot(all_meth_values_long_sel, aes(x = Site, y = Methylation_Value)) +
  geom_boxplot(aes(group = Site_f, fill = species), outlier.shape = "") +
  facet_wrap(~SCMR, scale = "free_x", nrow = 5) +
  scale_fill_manual(values = color_species) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Methylation values Atlantic Cod (AC), Australasian Snapper (AS), European Hake (EH), Zebrafish (ZF) (human rgenome)")
