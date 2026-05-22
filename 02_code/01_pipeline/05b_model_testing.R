# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: Epigenetic model testing
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026-02

#### Settings ####
data_folder <- "01_data/"
intermediate_folder <- "01_data/03_intermediate/07_data_comparison/"
dir.create(intermediate_folder, recursive = TRUE, showWarnings = FALSE)

load(paste0(data_folder, "04_metadata/color_palettes.RData"))

#### Preparation ####
library(tidyverse)
library(tidymodels)
library(ggtext)
library(ggpattern)
library(patchwork)

#### Loading data ####
load("01_data/03_intermediate/06_model_creation/all_data.Rdata")
load("01_data/03_intermediate/06_model_creation/models.Rdata")

#### Helper functions ####
message("Helper functions...")
compute.metrics <- function(actual, predicted) {
  # Standardize to vectors in case they come as 1-column matrices
  actual <- as.vector(actual)
  predicted <- as.vector(predicted)
  
  data.frame(
    R   = round(cor(predicted, actual, method = "pearson"), 2),
    MSE = round(mean((predicted - actual)^2), 3),
    MAE = round(mean(abs(predicted - actual)), 3),
    N   = length(actual)
  )
}
plot.prediction.results <- function(df, metrics, title, colpal, x_lim, y_lim, x_title, y_title) {
  anno_x <- x_lim[1]
  anno_y <- y_lim[2] 
  
  ggplot(df, aes(x = age, y = age_predicted, color = species)) +
    geom_point(size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
    geom_smooth(method = "lm", se = TRUE, color = "black", fill = "grey85") +
    scale_color_manual(values = colpal) +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    labs(
      title = title,
      y = y_title,
      x = x_title,
      color = "Species"
    ) +
    annotate("richtext", 
             x = anno_x, y = anno_y, 
             label = paste0("**R** = ", metrics$R, "<br>**MAE** = ", metrics$MAE), 
             size = 3.5, hjust = 0, vjust = 1, fill = NA, label.color = NA) 
}

plot.error.comparison <- function(ae_train, ae_test, y_lim = c(0, 0.125), y_title) {
  df_plot <- rbind(
    data.frame(AE = ae_train$AE, type = "Training"),
    data.frame(AE = ae_test$AE, type = "Testing")
  )
  df_plot$type <- factor(df_plot$type, levels = c("Training", "Testing"))

  ggplot(df_plot, aes(x = type, y = AE, fill = type, pattern = type)) +
    geom_boxplot_pattern(
      position = position_dodge(width = .9),
      outlier.shape = NA,
      pattern_fill = "transparent",
      pattern_color = "gray10"
    ) + 
    scale_fill_manual(values = c("Training" = "grey95", "Testing" = "grey95")) +
    scale_pattern_manual(values = c("stripe", "circle")) + 
    coord_cartesian(ylim = y_lim) +
    labs(y = y_title, x = "Data", pattern = "Data", fill = "Data") + 
    theme(legend.position = "none") 
}

evaluate.model <- function(model, X_train, Y_train, X_test, Y_test, species_train,
                           species_test, transform = "none", col = color_species, 
                           plot_title = "Model evaluation:", 
                           type = "relative",
                           y_lim = c(0,.3), x_lim = c(0,.3), 
                           ae_lim = c(0, 0.125), 
                           CpGs = "not defined", s = NA) {
  
  # Define labs based on type
  if (type == "chronological") {
    mod_x_title <- "Chronological age (years)"
    mod_y_title <- "Estimated age (years)"
    ae_y_title  <- "Absolute Error (years)"
  } else {
    mod_x_title <- "Relative age"
    mod_y_title <- "Estimated age"
    ae_y_title  <- "Absolute Error"
  }
  
  # transform actual values to original scale
  if (transform == "log") {
    Y_train <- exp(Y_train)
    Y_test  <- exp(Y_test)
  } else if (transform == "double_log" || transform == TRUE) {
    Y_train <- exp(-exp(-Y_train))
    Y_test  <- exp(-exp(-Y_test))
  }
  
  get_preds <- function(X) {
    p <- if (!is.na(s)) predict(model, X, s = s) else predict(model, X)
    p <- as.vector(p) 
    
    if (transform == "log") {
      p <- exp(p)
    } else if (transform == "double_log" || transform == TRUE) {
      p <- exp(-exp(-p))
    }
    return(p)
  }
  
  preds_train <- get_preds(X_train)
  preds_test  <- get_preds(X_test)
  
  metrics_train <- compute.metrics(Y_train, preds_train)
  metrics_test  <- compute.metrics(Y_test, preds_test)
  
  ae_train <- data.frame(AE = round(abs(preds_train - Y_train), 4))
  ae_test  <- data.frame(AE = round(abs(preds_test - Y_test), 4))
  
  metrics_ttest <- t.test(ae_train$AE, ae_test$AE)
  
  # clean up
  df_train <- data.frame(age_predicted = preds_train, age = Y_train, species = species_train)
  df_test  <- data.frame(age_predicted = preds_test, age = Y_test, species = species_test)
  
  # Plot
  plot_train <- plot.prediction.results(
    df_train, metrics_train, paste(plot_title, "(Training)"), 
    col, x_lim, y_lim, 
    x_title = mod_x_title, y_title = mod_y_title
  )
  
  plot_test <- plot.prediction.results(
    df_test, metrics_test, paste(plot_title, "(Testing)"), 
    col, x_lim, y_lim,
    x_title = mod_x_title, y_title = mod_y_title
  )
  
  plot_ae_comparison <- plot.error.comparison(
    ae_train, ae_test, 
    y_lim = ae_lim, 
    y_title = ae_y_title
  )
  
  return(list(
    metrics_train = metrics_train, 
    metrics_test  = metrics_test, 
    plot_train    = plot_train, 
    plot_test     = plot_test, 
    plot_ae       = plot_ae_comparison,
    values_AE_train = ae_train, 
    values_AE_test  = ae_test, 
    t_test        = metrics_ttest
  ))
}

# setting ggplot theme 
theme_custom <- theme_minimal() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(color = "black", linewidth = 1),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 12), 
            axis.ticks = element_line(color = "black"),
            axis.ticks.length = unit(0.1, "cm"), 
            axis.title = element_text(size = 10), 
            legend.ticks = element_blank(), 
            legend.title = element_text(face = "bold"))

theme_set(theme_custom)

#### Model evaluation ####
message("Model evaluation...")
# evaluation
mlm_chron_eval <- evaluate.model(
  mlm_chron, data_train, metadata_train$age, data_test, metadata_test$age, 
  metadata_train$species, metadata_test$species, 
  col = color_species, plot_title = "Chronological Age", 
  y_lim = c(0, 12), x_lim = c(0, 12), ae_lim = c(0, 3.5), CpGs = "not defined", s = NA, type = "chronological")

mlm_rel_eval <- evaluate.model(
  mlm_rel, data_train, metadata_train$rel_age, data_test, metadata_test$rel_age, 
  metadata_train$species, metadata_test$species, 
  col = color_species, plot_title = "Relative Age", 
  y_lim = c(0, .3), x_lim = c(0, .3), ae_lim = c(0, .15), CpGs = "not defined", s = NA, type = "relative")

# comparison plot
mlm_comparison_plot <- mlm_chron_eval$plot_train + mlm_chron_eval$plot_test + mlm_chron_eval$plot_ae +
  mlm_rel_eval$plot_train + mlm_rel_eval$plot_test + mlm_rel_eval$plot_ae +
  plot_layout(nrow = 2, guides = "collect", width = c(3,3,1)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(mlm_comparison_plot, file = "03_results/01_figures/05_mlm_comparison.png", width = 10, height = 7)
ggsave(mlm_comparison_plot, file = "03_results/01_figures/05_mlm_comparison.pdf", width = 10, height = 7)

# Log transformed models
# chronological age
mlm_chron_log_eval <- evaluate.model(
  mlm_chron_log, data_train, log(metadata_train$age), data_test, log(metadata_test$age), 
  metadata_train$species, metadata_test$species, 
  transform = "log", # <--- ADD THIS
  col = color_species, plot_title = "log-transformed", 
  y_lim = c(0, 12), x_lim = c(0, 12), ae_lim = c(0, 3.5), CpGs = "not defined", s = NA, type = "chronological")

# relative age
mlm_rel_log_eval <- evaluate.model(
  mlm_rel_log, data_train, -log(-log(metadata_train$rel_age)), data_test, -log(-log(metadata_test$rel_age)), 
  metadata_train$species, metadata_test$species, 
  transform = "double_log", # <--- ADD THIS
  col = color_species, plot_title = "log-transformed", 
  y_lim = c(0, .3), x_lim = c(0, .3), ae_lim = c(0, .15), CpGs = "not defined", s = NA, type = "relative")

# comparison plot
mlm_comparison_plot_log <- mlm_chron_log_eval$plot_train + mlm_chron_log_eval$plot_test + mlm_chron_log_eval$plot_ae +
  mlm_rel_log_eval$plot_train + mlm_rel_log_eval$plot_test + mlm_rel_log_eval$plot_ae +
  plot_layout(nrow = 2, guides = "collect", width = c(3,3,1)) + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(mlm_comparison_plot_log, file = "03_results/01_figures/05_mlm_comparison_log.png", width = 10, height = 7)
ggsave(mlm_comparison_plot_log, file = "03_results/01_figures/05_mlm_comparison_log.pdf", width = 10, height = 7)