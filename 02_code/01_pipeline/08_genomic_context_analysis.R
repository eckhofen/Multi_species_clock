# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: SMR gene analysis gene analysis
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
# correlation data
load("01_data/03_intermediate/05_correlation_data/cor_all.RData")

# annotated CpGs
load("01_data/03_intermediate/08_annotation/SMR_annotated.RData")
load("01_data/03_intermediate/08_annotation/Selected_annotated.RData")

# annotated CpGs counts
load("01_data/03_intermediate/08_annotation/SMR_counts.RData")
load("01_data/03_intermediate/08_annotation/Selected_counts.RData")

# genes
load("01_data/03_intermediate/08_annotation/SMR_genes.RData")
load("01_data/03_intermediate/08_annotation/Selected_genes.RData")

# colors
load("01_data/04_metadata/color_palettes.RData")

#### Plotting #### 
# ggplot theme settings
source("02_code/02_helpers/plot_style.R")

# define annotation type levels
type_names <- c("islands", "shores", "shelves", "open sea", "1to5kb", "promoters", "5UTRs", "exons", "introns", "3UTRs", "lncrna")
type_names <- factor(type_names, levels = c("islands", "shores", "shelves", "open sea", "1to5kb", "promoters", "5UTRs", "exons", "introns", "3UTRs", "lncrna", "other"))

# combine counts
SMR_counts$source <- factor("SMR", levels = c("SMR", "Selected"), ordered = TRUE)
SMR_counts$annot.type <- factor(SMR_counts$annot.type, levels = type_names)

selected_counts$source <- factor("Selected", levels = c("SMR", "Selected"), ordered = TRUE)
selected_counts$annot.type <- factor(selected_counts$annot.type, levels = type_names)

total_counts <- rbind(SMR_counts, selected_counts) %>%
  dplyr::select(-species) %>%
  group_by(annot.type, source, context) %>%
  summarise(n = sum(n), .groups = "drop")

# Summarise stats
cpg_context_wider <- total_counts %>%
  filter(context == "cpg") %>%
  pivot_wider(names_from = source, values_from = n) %>%
  mutate(
    pct_SMR = SMR / sum(SMR) * 100,
    pct_Selected = Selected / sum(Selected) * 100
  )

genic_context_wider <- total_counts %>%
  filter(context == "gene") %>%
  pivot_wider(names_from = source, values_from = n) %>%
  mutate(
    pct_SMR = SMR / sum(SMR) * 100,
    pct_Selected = Selected / sum(Selected) * 100
  )

# Test if different between selected and SMR
n_fisher_replicates <- 100000
# CpG context
chisq.test(cpg_context_wider %>% select(SMR, Selected))
# X-squared = 3.984, df = 3, p-value = 0.2632

set.seed(1999)
fisher.test(as.matrix(cpg_context_wider %>% select(SMR, Selected)), simulate.p.value = TRUE, B = n_fisher_replicates)
# p-value = 0.2812

chisq.test(genic_context_wider %>% select(SMR, Selected)) # Numbers are too small in some categories
# X-squared = 6.1334, df = 6, p-value = 0.4084

set.seed(1999)
fisher.test(as.matrix(genic_context_wider %>% select(SMR, Selected)), simulate.p.value = TRUE, B = n_fisher_replicates)
# p-value = 0.2878


# Comparison bar plot

total_counts_rev <- total_counts
total_counts_rev$annot.type <- factor(total_counts_rev$annot.type, levels = rev(type_names))
total_counts_rev$source <- factor(total_counts_rev$source, levels = rev(c("SMR", "Selected")))

p_comp <- total_counts_rev %>%
  group_by(source, context) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(y = source, x = pct, fill = annot.type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ context, ncol = 1, axes = "all") +
  labs(x = "Total CpGs (%)", fill = "Genomic context", y = "") +
  scale_fill_manual(values = color_types) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_discrete(expand = expansion(mult = c(0.25, 0.25))) +
  theme(strip.text = element_blank())

## SMR characteristics

# total
p_SMR_counts <- ggplot(total_counts, aes(x = annot.type, y = n, fill = annot.type, group = context)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.1), show.legend = FALSE) +
  geom_label(aes(label = n, fill = NULL), hjust = .5, vjust = -.3, show.legend = FALSE, border.color = NA) +
  facet_grid(source~., scales = "free", space = "free_x") +
  scale_fill_manual(values = color_types) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(y = "CpGs", x = "") + 
  theme(strip.text.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(paste0(figures_folder, "08_SMR_counts"), plot = p_SMR_counts, width = 8, height = 6)

# per species
p_SMR_counts_species <- ggplot(SMR_counts, aes(x = annot.type, y = n, fill = annot.type, group = context)) +
  geom_col(position = position_dodge2(preserve = "single", padding = 0.1), show.legend = FALSE) +
  geom_label(aes(label = n, fill = NULL), hjust = .5, vjust = -.3, show.legend = FALSE, color = "black", border.color = NA) +
  facet_grid(species ~ context, scales = "free", space = "free_x") +
  scale_fill_manual(values = color_types) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(y = "CpGs per species", x = "")

save_plot(paste0(figures_folder, "08_SMR_counts_species"), plot = p_SMR_counts_species, width = 8, height = 6)

# Pie chart

p_SMR_pie <- ggplot(total_counts, aes(x = 3, y = n, fill = annot.type)) +
  geom_bar(stat = "identity", aes(alpha = source), show.legend = c(alpha = FALSE)) +
  scale_x_continuous(limits = c(1, NA)) +
  facet_wrap(source ~ context, scales = "free", nrow = 2, strip.position = "bottom") +
  geom_label(aes(x = 3, y = n, label = n, fill = annot.type), inherit.aes = FALSE, position = position_stack(vjust = .5), show.legend = FALSE, size = 3) +
  coord_polar(theta = "y", start = 0) +
  labs(x = "", y = "") +
  scale_fill_manual(values = color_types) +
  scale_alpha_manual(values = c("SMR" = .5, "Selected" = 1)) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(.5, "cm"))

save_plot(paste0(figures_folder, "08_SMR_pie"), plot = p_SMR_pie, width = 8, height = 6)


## Selected CpGs

# per species
p_selected_counts_species <- ggplot(selected_counts) +
  geom_col(aes(x = annot.type, y = n, fill = annot.type, group = context), 
           position = position_dodge2(preserve = "single", padding = 0.1), show.legend = FALSE) +
  geom_label(aes(x = annot.type, y = n, label = n, fill = NULL, group = context),
             hjust = .5, vjust = -.3, show.legend = FALSE, border.color = NA) +
  facet_grid(species ~ context, scales = "free", space = "free_x") +
  scale_fill_manual(values = color_types) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(y = "Selected CpGs for multispecies clock", x = "")

save_plot(paste0(figures_folder, "08_selected_counts_species"), plot = p_selected_counts_species, width = 8, height = 6)

# pie chart
p_selected_pie <- ggplot(selected_counts, aes(x = 3, y = n, fill = annot.type)) +
  geom_bar(stat = "identity", aes(alpha = source), show.legend = c(alpha = FALSE)) +
  scale_x_continuous(limits = c(1, NA)) +
  facet_grid(context ~ species, scales = "free") +
  geom_label(aes(x = 3, y = n, label = n, fill = NULL), inherit.aes = FALSE, position = position_stack(vjust = .5), show.legend = FALSE, size = 3) +
  coord_polar(theta = "y", start = 0) +
  labs(x = "", y = "") +
  scale_fill_manual(values = color_types) +
  scale_alpha_manual(values = c("SMR" = .5, "Selected" = 1)) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        strip.text.y = element_blank(),
        strip.text = element_text(color = "black", face = "bold", size = 12),
        legend.key.size = unit(.5, "cm"))

save_plot(paste0(figures_folder, "08_selected_pie"), plot = p_selected_pie, width = 8, height = 4)
  
#### Genomic context specific analysis ####

# SMRs
cpgs_SMR_annotated <- cpgs_SMR_annotated %>% 
as.data.frame() %>%
# clean up names
mutate(annot.type = gsub("hg38_", "", annot.type)) %>%
# add context based on annot.type names
mutate(context = 
case_when(
    grepl("gene", annot.type) ~ "gene",
    grepl("lncrna", annot.type) ~ "gene",
    grepl("cpg", annot.type) ~ "cpg",
    TRUE ~ "other"
)) %>%
mutate(annot.type = gsub("genes_", "", annot.type)) %>%
mutate(annot.type = gsub("cpg_", "", annot.type)) %>%
mutate(annot.type = gsub("_gencode", "", annot.type)) %>%
mutate(annot.type = 
    case_when(
    annot.type == "inter" ~ "open sea",
    TRUE ~ annot.type)) %>%
mutate(annot.type = factor(annot.type, levels = type_names))

# Remove duplicated annotations (e.g. one position aligns to two introns)
cpgs_SMR_annotated_unique <- cpgs_SMR_annotated %>% 
  distinct(Site, annot.type, .keep_all = TRUE) %>%
  group_by(annot.type) %>%
  mutate(selected_ratio = sum(selected) / n()) %>%
  ungroup()

# Saving 
save(cpgs_SMR_annotated_unique, file = "01_data/03_intermediate/08_annotation/cpgs_SMR_annotated_unique.RData")

ggplot(cpgs_SMR_annotated_unique %>% filter(P_value < 0.05, selected == TRUE), aes(x = annot.type, y = abs(Correlation), fill = annot.type)) +
  geom_boxplot(show.legend = FALSE) +
  scale_fill_manual(values = color_types) +
  labs(y = "Correlation", x = "Annotation type")

p_gen_context_per_SMR_selected <- ggplot(cpgs_SMR_annotated_unique %>% filter(selected == TRUE), aes(x = annot.type, fill = annot.type)) +
  geom_bar() +
  facet_wrap(~SMR, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = color_types) +
  labs(y = "CpGs", x = "", fill = "") +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.position = "right",
        legend.key.size = unit(.5, "cm"))

save_plot(paste0(figures_folder, "08_gen_context_per_SMR_selected"), plot = p_gen_context_per_SMR_selected, width = 6, height = 9)

# Check selection ratio per group
selection_ratio <- cpgs_SMR_annotated_unique %>%
  group_by(annot.type, species) %>%
  summarise(selected_ratio = sum(selected) / n()) %>%
  mutate(selected_ratio_normalized = selected_ratio / max(selected_ratio))

ggplot(selection_ratio, aes(x = annot.type, y = selected_ratio, fill = annot.type)) +
  geom_col(show.legend = FALSE) +
  facet_grid(species ~ ., scales = "free_y") +
  scale_fill_manual(values = color_types) +
  labs(y = "Selection ratio", x = "Annotation type")

# Combined plots 
p_combined_count <- p_SMR_counts + p_comp + 
  plot_layout(ncol = 1, guides = "collect", axes = "collect_x", heights = c(3, 1)) + 
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 18, face = "bold"), legend.position = "none")

save_plot(paste0(figures_folder, "08_combined_count"), plot = p_combined_count, width = 6, height = 7)
save(p_combined_count, file = paste0(save_folder, "08_p_combined_count.RData"))

# Check genomic context per corelation direction
load("01_data/03_intermediate/08_annotation/SMR_annotated.RData")
load("01_data/03_intermediate/08_annotation/Selected_annotated.RData")

cpgs_SMR_annotated
cpgs_selected_annotated

cpgs_per_cor <- cpgs_selected_annotated %>% as_tibble() %>%
  mutate(cor_dir = if_else(Correlation > 0, "positive", "negative")) %>%
  count(cor_dir, annot.type, significant) %>%
  # clean up names
        mutate(annot.type = gsub("hg38_", "", annot.type)) %>%
    # add context based on annot.type names
        mutate(context = 
        case_when(
            grepl("gene", annot.type) ~ "gene",
            grepl("lncrna", annot.type) ~ "gene",
            grepl("cpg", annot.type) ~ "cpg",
            TRUE ~ "other"
        )) %>%
        mutate(annot.type = gsub("genes_", "", annot.type)) %>%
        mutate(annot.type = gsub("cpg_", "", annot.type)) %>%
        mutate(annot.type = gsub("_gencode", "", annot.type)) %>%
        mutate(annot.type = 
           case_when(
            annot.type == "inter" ~ "open sea",
            TRUE ~ annot.type)) %>%
        group_by(cor_dir, context, significant) %>%
        mutate(pct = n / sum(n)) %>%
        ungroup()

cpgs_per_cor$annot.type <- factor(cpgs_per_cor$annot.type, levels = rev(type_names))


p_cpg_cor <- ggplot(cpgs_per_cor %>% filter(significant == TRUE), aes(x = cor_dir, y = pct, fill = annot.type)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ context, ncol = 1, axes = "all") +
  labs(y = "Total CpGs (%)", fill = "Genomic context", x = "Age-correlation direction") +
  scale_fill_manual(values = color_types) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  scale_x_discrete(expand = expansion(mult = c(.5, .5))) +
  theme(strip.text = element_blank())

save_plot(paste0(figures_folder, "08_cpg_cor"), plot = p_cpg_cor, width = 4, height = 7)

# Test if distributions are different between correlation directions
cor_cpgs <- cpgs_per_cor %>% filter(significant) %>% select(-pct) %>% pivot_wider(names_from = cor_dir, values_from = n, values_fill = 0)

cor_cpgs_genic <- cor_cpgs %>% filter(context == "gene") %>% select(negative, positive)
cor_cpgs_cpg <- cor_cpgs %>% filter(context == "cpg") %>% select(negative, positive)

# Genic context
set.seed(1999)
fisher.test(as.matrix(cor_cpgs_genic), simulate.p.value = TRUE, B = n_fisher_replicates)
# p-value = 0.00494

# cpg context
set.seed(1999)
fisher.test(as.matrix(cor_cpgs_cpg), simulate.p.value = TRUE, B = n_fisher_replicates)
# p-value = 0.2783

message("Pipeline complete.")
