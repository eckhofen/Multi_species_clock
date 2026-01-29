
# Metadata ----------------------------------------------------------------
# Project: Piscine Multispecies Epigenetic Clock
# Description: Overview of data
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2025 06

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


# Analysis ----------------------------------------------------------------
## General analysis  ------------------------------------------------------

# Counts per species 
n_per_species <- meta_data %>% 
  group_by(species, common_name) %>% 
  summarise(n = n())

## Plots ------------------------------------------------------------------
### Plot settings -----------------------------------------------------------
# Custom theme settings
theme_custom <- function() {
  theme(panel.background = element_blank(), 
        plot.title = element_text(size = 18, face = "bold", hjust = .5), 
        axis.title = element_text(size = 14),
        panel.grid = element_blank(), 
        axis.line = element_line(),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10))
}

# Settings for overview plot
plot_overview <- function(p) {
  # Getting y position for number of species
  y_var_name <- rlang::as_name(rlang::get_expr(p$mapping$y))
  y_values <- p$data[[y_var_name]]
  max_y_value <- max(y_values, na.rm = TRUE)
  
  p +
    geom_jitter(aes(shape = sex),
                inherit.aes = TRUE, 
                cex = 2.5,
                alpha = .8,
                position = position_jitter(seed = 1999,
                                           width = .3,
                                           height = .005)) +
    geom_text(data = n_per_species, aes(x = species, y = 0, label = n),
              inherit.aes = TRUE,
              size = 4.5,
              fontface = "bold", 
              nudge_y = -max_y_value*.07,
              show.legend = FALSE) +
    geom_violin(color = NA,
                alpha = .2,
                show.legend = FALSE) +

    # Adding mean
    stat_summary(fun = mean,
                 geom = "point",
                 shape = 18,
                 size = 5,
                 color = "black",
                 alpha = .9,
                 show.legend = FALSE) +
    # Adding number of samples
    labs(x = "Species",
         color = "Species",
         fill = "Species",
         shape = "Sex") +
    scale_color_manual(values = color_species_common) +
    scale_fill_manual(values = color_species_common)
}


### Creating plots --------------------------------------------------------
# Overview for age 
fig_overview_age <- meta_data %>% 
  ggplot(aes(x = species, y = age, color = common_name, fill = common_name)) %>% 
  plot_overview() + 
  theme_custom() + 
  labs(y = "Age (years)", title = "Sample age distribution")

# Saving plot
ggsave(filename = paste0(results_folder, "01_fig_overview_age.PDF"), plot = fig_overview_age, width = 7, height = 7)

# Overview for relative age
fig_overview_age_rel <- meta_data %>% 
  ggplot(aes(x = species, y = age_rel, color = common_name, fill = common_name)) %>% 
  plot_overview() + 
  theme_custom() + 
  labs(y = "Age (relative to max. lifespan)", title = "Sample age (relative) distribution")

# Saving plot
ggsave(filename = paste0(results_folder, "01_fig_overview_age_rel.PDF"), plot = fig_overview_age_rel, width = 7, height = 7)

### Combining plots ---------------------------------------------------------
# Combined overview figure
fig_overview_combined <- 
  fig_overview_age + fig_overview_age_rel +
  plot_layout(nrow = 1, guides = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

# Saving plot
ggsave(filename = paste0(results_folder, "01_fig_overview_combined.PDF"), plot = fig_overview_combined, width = 12, height = 7)
