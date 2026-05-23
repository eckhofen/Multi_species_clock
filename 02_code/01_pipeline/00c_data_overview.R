# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: Overview of data
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026-02

# Settings ----------------------------------------------------------------
library(tidyverse)
library(ggforce)
library(patchwork)

data_folder <- "01_data/"
results_folder <- "03_results/01_figures/"
dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)

# Loading data ------------------------------------------------------------

# Loading colors 
load(paste0(data_folder, "04_metadata/color_palettes.RData"))

# Loading metadata
load(paste0(data_folder, "04_metadata/metadata.RData"))

# plot style
source("02_code/02_helpers/plot_style.R")

#### Analysis #### 

# Counts per species 
n_per_species <- meta_data %>% 
  group_by(species, common_name) %>% 
  summarise(n = n())

#### Plots ####
# Settings for overview plot
plot_overview <- function(p) {
  # Getting y position for number of species
  y_var_name <- rlang::as_name(rlang::get_expr(p$mapping$y))
  y_values <- p$data[[y_var_name]]
  max_y_value <- max(y_values, na.rm = TRUE)
  
  p +
    geom_jitter(,
                inherit.aes = TRUE, 
                cex = 1.5,
                position = position_jitter(height = .01, width = .1),
                alpha = .8,
                show.legend = c(color = FALSE, alpha = TRUE)) +
    geom_text(data = n_per_species, aes(x = species, y = 0, label = n),
              inherit.aes = TRUE,
              size = 4.5,
              fontface = "bold", 
              nudge_y = -max_y_value*.07,
              show.legend = FALSE) +
    geom_violin(color = NA,
                alpha = .2,
                show.legend = (color = FALSE)) +

    # Adding mean
    stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 2.5,
                 color = "black",
                 alpha = .9,
                 show.legend = (color = FALSE)) +
    # Adding number of samples
    labs(x = "Species",
         color = "Species",
         fill = "Species") +
    scale_color_manual(values = color_species_common) +
    scale_fill_manual(values = color_species_common)
}


### Creating plots --------------------------------------------------------
# Overview for age 
fig_overview_age <- meta_data %>% 
  ggplot(aes(x = species, y = age, color = common_name, fill = common_name)) %>% 
  plot_overview() + 
  theme_custom + 
  labs(y = "Age (years)", title = "Sample age distribution")

# Saving plot
save_plot(filename = paste0(results_folder, "00_fig_overview_age"), plot = fig_overview_age, width = 7, height = 7)

# Overview for relative age
fig_overview_age_rel <- meta_data %>% 
  ggplot(aes(x = species, y = age_rel, color = common_name, fill = common_name)) %>% 
  plot_overview() + 
  theme_custom + 
  labs(y = "Age (relative to max. lifespan)", title = "Sample age (relative) distribution")

# Saving plot
save_plot(filename = paste0(results_folder, "00_fig_overview_age_rel"), plot = fig_overview_age_rel, width = 7, height = 7)

### Combining plots ---------------------------------------------------------
# Combined overview figure
fig_overview_combined <- 
  fig_overview_age + fig_overview_age_rel +
  plot_layout(nrow = 1, guides = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"), legend.position = "right", plot.title = element_blank())

# Saving plot
save_plot(filename = paste0(results_folder, "00_fig_overview_combined"), plot = fig_overview_combined, width = 9, height = 4)

