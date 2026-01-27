# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: June 2024

#### Overview ####
# retrieving age data on all samples for all species and plotting the distribution via boxplots

#### Settings ####
setwd("/Users/macether/Documents/2 - Studium/1 - Master/ZZ - Thesis/Repo_Multispecies_clock/Multi_species_clock/")
extension <- ".pdf"

#### Preparation ####
library(ggplot2)
library(ggforce)
library(dplyr)
library(patchwork)

#### Loading data ####
load("000_data/000_metadata/metadata_age.RData")

df_all <- rbind(AC_df, AS_df, EH_df, ZF_df)

#### Plotting ####

## setting up colors
# selected color which are colorblind-friendly
colpal_CB_c <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")

color_species_df <- data.frame(species = as.factor(c("AC","AS","EH","JM","ZF")), color = colpal_CB_c[c(1, 5, 3, 7, 8)])
color_species <- setNames(color_species_df$color, color_species_df$species)

## species counts
# calculate counts for each species
species_counts <- df_all %>%
  group_by(species) %>%
  summarise(count = n(), .groups = 'drop')

### plots
## relative age
plot_age <- ggplot(df_all, aes(y = rel_age, color = species)) +
  geom_sina(aes(x = species, shape = sex), cex = 3,  alpha = .7) +
  geom_boxplot(aes(x = species), fill = NA) +
  scale_color_manual(values = color_species) +
  theme_classic() +
  geom_text(data = species_counts, aes(x = species, y = -.015, label = count, color = species), size = 3.5) +
  labs(x = "Species", y = "Relative age", color = "Species", shape = "Sex")

ggsave(filename = paste0("002_plots/000_rel_age_distribution", extension), plot_age, width = 7, height = 7)

## chronological age
plot_direct_age <- ggplot(df_all, aes(y = age, color = species)) +
  geom_sina(aes(x = species, shape = sex), cex = 3,  alpha = .7) +
  geom_boxplot(aes(x = species), fill = NA) +
  scale_color_manual(values = color_species) +
  theme_classic() +
  geom_text(data = species_counts, aes(x = species, y = -.5, label = count, color = species), size = 3.5) +
  labs(x = "Species", y = "Chronological age (years)", color = "Species", shape = "Sex")

ggsave(filename = paste0("002_plots/000_age_distribution", extension), plot_direct_age, width = 7, height = 7)

## combined plot
plot_relative_direct_age <- plot_direct_age + plot_age +
  plot_layout(nrow = 1, guides = "collect") + 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

ggsave(filename = paste0("002_plots/000_age_distribution_both", extension), plot_relative_direct_age, width = 10, height = 5)
