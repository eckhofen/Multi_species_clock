# Metadata ----------------------------------------------------------------
# Project: Shared co-mehtylation regions (SCMRs)
# Description: SCMR gene analysis gene analysis
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2025 06

#### Libraries ####
library(tidyverse)
library(patchwork)
library(ggrepel)
library(ggtext)
library(ggnewscale)

#### Setting paths ####
figures_folder <- "03_results/01_figures/"
save_folder <- "01_data/03_intermediate/08_annotation/"

#### Loading data ####
# annotation
load("01_data/03_intermediate/08_annotation/annotation.RData")

# annotated CpGs counts
load("01_data/03_intermediate/08_annotation/SCMR_counts.RData")
load("01_data/03_intermediate/08_annotation/Selected_counts.RData")

# genes
load("01_data/03_intermediate/08_annotation/SCMR_genes.RData")
load("01_data/03_intermediate/08_annotation/Selected_genes.RData")

# colors
load("01_data/04_metadata/color_palettes.RData")

#### Plotting #### 
# color palette for annotation contexts
color_contexts <- c("gene" = "#1f77b4", "lncrna" = "#ff7f0e", "cpg" = "#2ca02c", "other" = "#d62728")

# color palette for annotation type
type_names <- c("islands", "shores", "shelves", "open sea", "1to5kb", "promoters", "5UTRs", "exons", "introns", "3UTRs", "lncrna")
type_names <- factor(type_names, levels = c("islands", "shores", "shelves", "open sea", "1to5kb", "promoters", "5UTRs", "exons", "introns", "3UTRs", "lncrna", "other"))
type_colors <- c('#d94701','#fd8d3c', '#fdbe85','#feedde', 
                 "#006d2c", "#31a354", "#74c476", "#a1d99b", "#c7e9c0", "#edf8e9", 
                 "#756bb1")
color_types <- setNames(type_colors, type_names)

# comine counts
scmr_counts$source <- factor("SCMR", levels = c("SCMR", "Selected"), ordered = TRUE)
scmr_counts$annot.type <- factor(scmr_counts$annot.type, levels = type_names)

selected_counts$source <- factor("Selected", levels = c("SCMR", "Selected"), ordered = TRUE)
selected_counts$annot.type <- factor(selected_counts$annot.type, levels = type_names)


total_counts <- rbind(scmr_counts, selected_counts) %>%
  dplyr::select(-species) %>%
  group_by(annot.type, source, context) %>%
  summarise(n = sum(n), .groups = "drop")

# ggplot theme settings
ggtheme <- theme_classic() +
  theme(strip.background = element_rect(fill = "grey95", color = NA), strip.text = element_text(color = "black", face = "bold"))

## SCMR characteristics

# total
p_scmr_counts <- ggplot(total_counts, aes(x = annot.type, y = n, fill = annot.type, group = context)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.1), show.legend = FALSE) +
  geom_label(aes(label = n), hjust = .5, vjust = -.3, show.legend = FALSE, border.color = NA) +
  facet_grid(source ~ context, scales = "free", space = "free_x") +
  ggtheme +
  scale_fill_manual(values = color_types) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(y = "CpGs per species", x = "") +
  theme(strip.text.x = element_blank(), strip.background.x = element_blank(), strip.background.y = element_rect(fill = "grey95", color = NA))

ggsave(file = paste0(figures_folder, "07_scmr_counts.pdf"), plot = p_scmr_counts, width = 8, height = 6)

# per species
p_scmr_counts_species <- ggplot(scmr_counts, aes(x = annot.type, y = n, fill = annot.type, group = context)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.1), show.legend = FALSE) +
  geom_label(aes(label = n), hjust = .5, vjust = -.3, show.legend = FALSE, color = "black", border.color = NA) +
  facet_grid(species ~ context, scales = "free", space = "free_x") +
  ggtheme +
  scale_fill_manual(values = color_types) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(y = "CpGs per species", x = "") +
  theme(strip.text.x = element_blank(), strip.background.x = element_blank(), strip.background.y = element_rect(fill = "grey95", color = NA))

ggsave(file = paste0(figures_folder, "07_scmr_counts_species.pdf"), plot = p_scmr_counts_species, width = 8, height = 6)


# Pie chart

p_scmr_pie <- ggplot(total_counts, aes(x = 3, y = n, fill = annot.type)) +
  geom_bar(stat = "identity", aes(alpha = source), show.legend = c(alpha = FALSE)) +
  scale_x_continuous(limits = c(1, NA)) +
  facet_wrap(source ~ context, scales = "free", nrow = 2, strip.position = "bottom") +
  geom_label(aes(x = 3, y = n, label = n, fill = annot.type), inherit.aes = FALSE, position = position_stack(vjust = .5), show.legend = FALSE, size = 3) +
  coord_polar(theta = "y", start = 0) +
  labs(x = "", y = "") +
  scale_fill_manual(values = color_types) +
  scale_alpha_manual(values = c("SCMR" = .5, "Selected" = 1)) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(.5, "cm"))

ggsave(file = paste0(figures_folder, "07_scmr_pie.pdf"), plot = p_scmr_pie, width = 8, height = 6)


## Selected CpGs

# per species
p_selected_counts_species <- ggplot(selected_counts) +
  geom_col(aes(x = annot.type, y = n, fill = annot.type, group = context), 
           position = position_dodge2(preserve = "single", padding = 0.1), show.legend = FALSE) +
  geom_label(aes(x = annot.type, y = n, label = n, fill = annot.type, group = context),
             hjust = .5, vjust = -.3, show.legend = FALSE, border.color = NA) +
  facet_grid(species ~ context, scales = "free", space = "free_x") +
  ggtheme +
  scale_fill_manual(values = color_types) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(y = "Selected CpGs for multispecies clock", x = "") +
  theme(strip.text.x = element_blank(), 
        strip.background.x = element_blank(), 
        strip.background.y = element_rect(fill = "grey95", color = NA))

ggsave(file = paste0(figures_folder, "07_selected_counts_species.pdf"), plot = p_selected_counts_species, width = 8, height = 6)

# pie chart
p_selected_pie <- ggplot(selected_counts, aes(x = 3, y = n, fill = annot.type)) +
  geom_bar(stat = "identity", aes(alpha = source), show.legend = c(alpha = FALSE)) +
  scale_x_continuous(limits = c(1, NA)) +
  facet_grid(context ~ species, scales = "free") +
  geom_label(aes(x = 3, y = n, label = n, fill = annot.type), inherit.aes = FALSE, position = position_stack(vjust = .5), show.legend = FALSE, size = 3) +
  coord_polar(theta = "y", start = 0) +
  labs(x = "", y = "") +
  scale_fill_manual(values = color_types) +
  scale_alpha_manual(values = c("SCMR" = .5, "Selected" = 1)) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(color = "black", face = "bold", size = 12),
        legend.key.size = unit(.5, "cm"))

ggsave(file = paste0(figures_folder, "07_selected_pie.pdf"), plot = p_selected_pie, width = 8, height = 4)
