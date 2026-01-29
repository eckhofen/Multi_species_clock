# Author: Gabriel Ecker-Eckhofen (gabriel.eckhofen@imbrsea.eu)
# Date: November 2024

# DISCLAIMER:
# none

#### Overview ####
## gene analysis of SCMR data on the human reference genome 

#### Settings ####
data_folder <- "01_data/"
intermediate_folder <- paste0(data_folder, "03_intermediate/")
results_folder <- "03_results/01_figures/"
save_folder <- paste0(intermediate_folder, "009_gene_analysis/")

dir.create(results_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)

palette_path <- paste0(data_folder, "04_metadata/color_palettes.RData")
if (file.exists(palette_path)) {
  load(palette_path)
}

if (!exists("color_compare")) {
  color_compare <- c("#005AB5", "#DC3220")
}

if (!exists("color_species")) {
  colpal_CB_c <- c("#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", "#AA4499", "#882255")
  color_species_df <- data.frame(
    species = as.factor(c("AC", "AS", "EH", "ZF")),
    color = colpal_CB_c[c(1, 5, 3, 8)]
  )
  color_species <- setNames(color_species_df$color, color_species_df$species)
}

#### Preparation ####
library(tidyverse)
library(tidymodels)
library(tibble)
library(dplyr)
library(ggplot2)
library(patchwork)
library(annotatr)
library(rtracklayer)
library(AnnotationHub)


#### Loading data ####
### Annotation
## (A) loading annotation file (downloaded from NCBI)
annotation_folder_new <- paste0(intermediate_folder, "008_annotation/")
annotation_folder <- annotation_folder_new

annotation_gff_path <- paste0(annotation_folder, "GRCh38_p14_annotation.gff")
annotation_human <- import.gff(annotation_gff_path, genome = "GRCh38")
str(annotation_human)


## (B) getting annotation with AnnotationHub (alternative) 
hub <- AnnotationHub()
source_hub <- query(hub, c("UCSC","CpG Islands", "hg38")) # does not exist or work
source_hub <- query(hub, c("CpG Islands"))
source_hub <- query(hub, c("hg38"))
grep("Cp", source_hub$title, value = TRUE)

cpg_data <- hub[["AH5086"]]
sort(unique(mcols(cpg_data)))

## (C) annotation from UCSC for CpG Islands
CpG_isl <- read.table(paste0(annotation_folder, "UCSC_CpG_islands.gff"), sep = "\t", header = TRUE)

# transforming into GRanges object
anno_CpG <- GRanges(seqnames = CpG_isl$chrom, IRanges(start = CpG_isl$chromStart, end = CpG_isl$chromEnd))

# adding metadata
mcols(anno_CpG) <- CpG_isl[,-c(1:4)]

## changing name to seqinfo style (chr1 to NC000***) for the first 24 chromosomes
seqnames(anno_CpG)

# getting NC names for first 24 chromosomes
chr_names_NC <- grep("^NC", unique(seqnames(annotation_human)), value = TRUE) %>% .[1:24]

# getting chromosome name vector
chr_names_chr <-  sort(unique(seqnames(anno_CpG))) %>%
  .[1:24]

# adding levels to the GRanges object
seqlevels(anno_CpG) <- union(seqlevels(anno_CpG), chr_names_NC)

# making name map for NC format
name_map <- setNames(chr_names_NC, chr_names_chr)

## applying name map
current_seqnames <- as.character(seqnames(anno_CpG))
new_names <- ifelse(current_seqnames %in% as.character(chr_names_chr),
                             as.character(name_map[current_seqnames]),
                             current_seqnames)

# creating new granges object with new names
anno_CpG_NC <- GRanges(seqnames = new_names, ranges = ranges(anno_CpG), strand = "*", mcols = mcols(anno_CpG))

### loading CpGs
path_methylsites_new <- paste0(annotation_folder_new, "methylsites_all.RData")
load(path_methylsites_new)
methyl_sites_combined

# transforming df into GRanges 
CpGs <- GRanges(seqnames = methyl_sites_combined$chr_align, ranges = IRanges(methyl_sites_combined$pos_align, methyl_sites_combined$pos_align, 1))

# load correlation data for all CpGs
path_cor_all_new <- paste0(intermediate_folder, "005_correlation_data/cor_all.RData")
load(path_cor_all_new)
cor_all

# adding metadata and correlation test results as meta column 
mcols(CpGs) <- cbind(methyl_sites_combined[,-c(1:4)], cor_all)
CpGs

## selecting final CpGs
# loading selected CpGs
path_sel_cpg_new <- paste0(intermediate_folder, "005_correlation_data/all_mix_cor_CpG_common.RData")
load(path_sel_cpg_new)
all_mix_cor_CpG_common

# selecting final CpGs
CpGs_selected <- CpGs[CpGs$Site %in% all_mix_cor_CpG_common$Site]

#### Data annotation ####

### (A) UCSC annotation

## ALL CpGs
# finding overlapping annotations
overlaps <- findOverlaps(CpGs, annotation_human, type = "within")

# getting all annotations found for queried CpGs 
anno_index <- subjectHits(overlaps)
anno_hits <- annotation_human$type[anno_index]

# annotations only for genes
gene_hits <- annotation_human$gene[anno_index[which(annotation_human$type[anno_index] == "gene")]] %>% table()

lnc_RNA_hits <- annotation_human$product[anno_index[which(annotation_human$type[anno_index] == "lnc_RNA")]]

# exporting genes
write.csv(gene_hits, file = paste0(save_folder, "gene_list.csv"))

# counting the regions in the selected SMRs
CpG_regions <- split(anno_hits, queryHits(overlaps))

# only counting single matches
CpG_regions_uni <- lapply(CpG_regions, function(x) unique(x))

# adding correlation direction of CpGs (positive or negative) as dataframe
CpG_regions_pos <- CpG_regions_uni[CpGs$Correlation > 0] %>% 
  unlist() %>% table() %>% as.data.frame()
CpG_regions_neg <- CpG_regions_uni[CpGs$Correlation < 0] %>% 
  unlist() %>% table() %>% as.data.frame()

# saving into dataframe and discarding absence
CpG_regions_df <- data_frame(regions = CpG_regions_neg$.,
                                 pos_cor = CpG_regions_pos$Freq,
                                 neg_cor = CpG_regions_neg$Freq) %>% 
  mutate(check = c(pos_cor + neg_cor)) %>% 
  .[.$check > 0, ] %>% .[!.$regions %in% "mRNA",] 

# saving as dataframe which sums up both correlation groups
CpG_regions_all_df <- CpG_regions_df %>%
  mutate(CpGs = pos_cor + neg_cor) %>%
  dplyr::select(regions, CpGs)

# changing to long format for ggplot 
CpG_regions_df_long <- CpG_regions_df %>% 
  pivot_longer(cols = c(pos_cor, neg_cor), 
               names_to = "correlation", 
               values_to = "occurrences")

# plotting CpGs in regions
p_CpGs <- ggplot(CpG_regions_df_long, aes(x = regions, y = occurrences, fill = correlation)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = occurrences), position = position_dodge(width = 0.9), hjust = -.1, angle = 90, size = 2.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Regions", y = "CpGs", fill = "Correlation") + 
  scale_fill_manual(labels = c("pos_cor" = "Positive", "neg_cor" = "Negative"), 
                    values = c("pos_cor" = color_compare[1], "neg_cor" = color_compare[2])) +
  ylim(0, 1.1*max(CpG_regions_df_long$occurrences))
p_CpGs

ggsave(paste0(results_folder, "07_regions_all_cor_group.pdf"), p_CpGs, width = 8, height = 6)

# plotting combined 
p_CpGs_combined <- ggplot(CpG_regions_all_df, aes(x = regions, y = CpGs, fill = regions)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  geom_text(aes(label = CpGs), position = position_dodge(width = 0.9), hjust = -.1, angle = 90, size = 4.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12.5), legend.position = "none") +
  labs(x = "Regions", y = "CpGs") +
  ylim(0, 1.1*max(CpG_regions_all_df$CpGs))
p_CpGs_combined

ggsave(paste0(results_folder, "07_regions_all_CpGs_combined.pdf"), p_CpGs_combined, width = 8, height = 6)
 
# adding species of CpGs as dataframe
CpG_regions_AC <- CpG_regions_uni[CpGs$species == "AC"] %>% 
  unlist() %>% table() %>% as.data.frame()
CpG_regions_AS <- CpG_regions_uni[CpGs$species == "AS"] %>% 
  unlist() %>% table() %>% as.data.frame()
CpG_regions_EH <- CpG_regions_uni[CpGs$species == "EH"] %>% 
  unlist() %>% table() %>% as.data.frame()
CpG_regions_ZF <- CpG_regions_uni[CpGs$species == "ZF"] %>% 
  unlist() %>% table() %>% as.data.frame()

# saving into dataframe and discarding absence
CpG_regions_df_species <- data_frame(regions = CpG_regions_neg$.,
                                         AC = CpG_regions_AC$Freq,
                                         AS = CpG_regions_AS$Freq,
                                         EH = CpG_regions_EH$Freq,
                                         ZF = CpG_regions_ZF$Freq) %>% 
  .[!.$regions %in% "mRNA",] %>% mutate(check = c(AC + AS + EH + ZF)) %>% 
  .[.$check > 0, ]

# changing to long format for ggplot 
CpG_regions_df_species_long <- CpG_regions_df_species %>% 
  pivot_longer(cols = c(AC, AS, EH, ZF), 
               names_to = "correlation", 
               values_to = "occurrences")

# plotting CpGs in regions
p_CpGs_species <- ggplot(CpG_regions_df_species_long, aes(x = regions, y = occurrences, fill = correlation)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = occurrences), position = position_dodge(width = 0.9), 
            hjust = -0.1, size = 2.5, angle = 90) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Regions", y = "CpGs", fill = "Species") + 
  scale_fill_manual(values = color_species) +
  ylim(0,1.1*max(CpG_regions_df_species_long$occurrences))
p_CpGs_species

ggsave(paste0(results_folder, "07_regions_all_species.pdf"), p_CpGs_species, width = 8, height = 6)

## SELECTED CpGs
# finding overlapping annotations
overlaps_sel <- findOverlaps(CpGs_selected, annotation_human, type = "within")
overlaps_sel_pos <- findOverlaps(CpGs_selected[CpGs_selected$Correlation > 0], annotation_human, type = "within")
overlaps_sel_neg <- findOverlaps(CpGs_selected[CpGs_selected$Correlation < 0], annotation_human, type = "within")

# getting all annotations found for queried CpGs 
anno_index_sel <- subjectHits(overlaps_sel)
anno_hits_sel <- annotation_human$type[anno_index_sel]

# for pos. and neg. cor
anno_index_sel_pos <- subjectHits(overlaps_sel_pos)
anno_index_sel_neg <- subjectHits(overlaps_sel_neg)

### annotations for genes
gene_hits_sel <- annotation_human$gene[anno_index_sel[which(annotation_human$type[anno_index_sel] == "gene")]]

gene_hits_sel_pos <- annotation_human$gene[anno_index_sel_pos[which(annotation_human$type[anno_index_sel_pos] == "gene")]]

gene_hits_sel_neg <- annotation_human$gene[anno_index_sel_neg[which(annotation_human$type[anno_index_sel_neg] == "gene")]]

# count how many times each gene appears in each group of correlation direction
genes_sel_count_pos <- sapply(unique(gene_hits_sel), function(i) sum(gene_hits_sel_pos == i)) %>% as.data.frame()
genes_sel_count_neg <- sapply(unique(gene_hits_sel), function(i) sum(gene_hits_sel_neg == i)) %>% as.data.frame()

genes_sel_count <- data_frame(genes = rownames(genes_sel_count_pos),
                              pos_cor = genes_sel_count_pos$.,
                              neg_cor = genes_sel_count_neg$.)
# exporting genes
write.csv(table(gene_hits_sel), file = paste0(save_folder, "gene_list_sel.csv"))
write.csv(genes_sel_count, file = paste0(save_folder, "gene_list_count.csv"))

# compare with age related genes
GenAge_genes <- read_csv(paste0(annotation_folder, "genage_human.csv"))

genes_sel_count$genes %in% GenAge_genes$symbol 
names(gene_hits) %in% GenAge_genes$symbol 

## plotting affected genes

# long format
genes_sel_count_long <- pivot_longer(genes_sel_count, cols = c("pos_cor", "neg_cor"), names_to = "cor", values_to = "CpGs")

# plot
p_genes_sel <- ggplot(genes_sel_count_long, aes(y = genes, x = CpGs, fill = cor)) +
  geom_col(position = position_stack()) +
  scale_fill_manual(values = rev(color_compare), ) +
  labs(y = "Genes", x = "CpGs", fill = "Correlation") +
  scale_y_discrete(limits = rev(sort(genes_sel_count_long$genes))) +
  theme_minimal()
p_genes_sel

### regions
# counting the regions in the selected SMRs
CpG_sel_regions <- split(anno_hits_sel, queryHits(overlaps_sel))

# ridding of multiple exons/mRNAs etc.
CpG_sel_regions_uni <- lapply(CpG_sel_regions, function(x) unique(x))

# adding correlation direction of CpGs (positive or negative) as dataframe
CpG_sel_regions_pos <- CpG_sel_regions_uni[CpGs_selected$Correlation > 0] %>% 
  unlist() %>% table() %>% as.data.frame()
CpG_sel_regions_neg <- CpG_sel_regions_uni[CpGs_selected$Correlation < 0] %>% 
  unlist() %>% table() %>% as.data.frame()


# saving into dataframe and discarding absence
CpG_sel_regions_df <- data_frame(regions = CpG_sel_regions_neg$.,
                                 pos_cor = CpG_sel_regions_pos$Freq,
                                 neg_cor = CpG_sel_regions_neg$Freq) %>% 
  mutate(check = c(pos_cor + neg_cor)) %>% 
  .[.$check > 0, ] %>% .[!.$regions %in% "mRNA",] 

# saving both groups combined as dataframe 
CpG_sel_regions_comb_df <- CpG_sel_regions_df %>% 
  mutate(CpGs = pos_cor + neg_cor) %>% 
  dplyr::select(regions, CpGs)

# changing to long format for ggplot 
CpG_sel_regions_df_long <- CpG_sel_regions_df %>% 
  pivot_longer(cols = c(pos_cor, neg_cor), 
               names_to = "correlation", 
               values_to = "occurrences")

# plotting CpGs in regions
p_CpGs_sel <- ggplot(CpG_sel_regions_df_long, aes(x = regions, y = occurrences, fill = correlation)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = occurrences), position = position_dodge(width = 0.9), hjust = -.1, angle = 90, size = 2.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Regions", y = "CpGs", fill = "Correlation") + 
  scale_fill_manual(labels = c("pos_cor" = "Positive", "neg_cor" = "Negative"), 
                    values = c("pos_cor" = color_compare[1], "neg_cor" = color_compare[2])) +
  ylim(0,1.1*max(CpG_sel_regions_df_long$occurrences))
p_CpGs_sel

ggsave(paste0(results_folder, "07_regions_sel_cor_group.pdf"), p_CpGs_sel, width = 8, height = 6)


# plot of selected CpGs both groups combined
p_CpGs_sel_combined <- ggplot(CpG_sel_regions_comb_df, aes(x = regions, y = CpGs, fill = regions)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = CpGs), position = position_dodge(width = 0.9), hjust = -.1, angle = 90, size = 4.5) +
  theme_classic() +
  scale_fill_viridis_d(option = "cividis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12.5), legend.position = "none") +
  labs(x = "Regions", y = "CpGs") +
  ylim(0, 1.1*max(CpG_sel_regions_comb_df$CpGs))
p_CpGs_sel_combined

ggsave(paste0(results_folder, "07_regions_sel_CpGs_combined.pdf"), p_CpGs_sel_combined, width = 8, height = 6)

# adding species of CpGs as dataframe
CpG_sel_regions_AC <- CpG_sel_regions_uni[CpGs_selected$species == "AC"] %>% 
  unlist() %>% table() %>% as.data.frame()
CpG_sel_regions_AS <- CpG_sel_regions_uni[CpGs_selected$species == "AS"] %>% 
  unlist() %>% table() %>% as.data.frame()
CpG_sel_regions_EH <- CpG_sel_regions_uni[CpGs_selected$species == "EH"] %>% 
  unlist() %>% table() %>% as.data.frame()
CpG_sel_regions_ZF <- CpG_sel_regions_uni[CpGs_selected$species == "ZF"] %>% 
  unlist() %>% table() %>% as.data.frame()

# saving into dataframe and discarding absence
CpG_sel_regions_df_species <- data_frame(regions = CpG_sel_regions_neg$.,
                                 AC = CpG_sel_regions_AC$Freq,
                                 AS = CpG_sel_regions_AS$Freq,
                                 EH = CpG_sel_regions_EH$Freq,
                                 ZF = CpG_sel_regions_ZF$Freq) %>% 
  .[!.$regions %in% "mRNA",] %>% mutate(check = c(AC + AS + EH + ZF)) %>% 
  .[.$check > 0, ]

# changing to long format for ggplot 
CpG_sel_regions_df_species_long <- CpG_sel_regions_df_species %>% 
  pivot_longer(cols = c(AC, AS, EH, ZF), 
               names_to = "correlation", 
               values_to = "occurrences")

# plotting CpGs in regions
p_CpGs_sel_species <- ggplot(CpG_sel_regions_df_species_long, aes(x = regions, y = occurrences, fill = correlation)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = occurrences), position = position_dodge(width = 0.9), 
            hjust = -0.1, size = 2.5, angle = 90) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Regions", y = "CpGs", fill = "Species") + 
  scale_fill_manual(values = color_species) +
  ylim(0,1.1*max(CpG_sel_regions_df_species_long$occurrences))
p_CpGs_sel_species

ggsave(paste0(results_folder, "07_regions_sel_species.pdf"), p_CpGs_sel_species, width = 8, height = 6)

#### plots ####
p_all_CpGs <- p_CpGs + p_CpGs_species + p_CpGs_sel + p_CpGs_sel_species +
  plot_layout(guides = "collect") &
  theme(legend.position = "top") 
p_all_CpGs

ggsave(paste0(results_folder, "07_regions_all.pdf"), p_all_CpGs, width = 10, height = 7)

p_genes_sel <- p_genes_sel + theme_classic()

ggsave(paste0(results_folder, "07_genes_sel_cor.pdf"), p_genes_sel, width = 6, height = 8)

 ## (C)
overlaps_CGI <- findOverlaps(CpGs, anno_CpG_NC, type = "within")
overlaps_CGI_sel <- findOverlaps(CpGs_selected, anno_CpG_NC, type = "within")

anno_index_CGI <- subjectHits(overlaps_CGI)
anno_index_CGI_sel <- subjectHits(overlaps_CGI_sel)

sort(table(anno_CpG_NC$mcols.name[anno_index_CGI]))
sort(table(anno_CpG_NC$mcols.name[anno_index_CGI_sel]))







## renaming seqnames with chromosomes
# chr_names_NC <- grep("^NC", unique(seqnames(annotation_human)), value = TRUE) %>% .[1:24]
# chr_names_chr <-  sort(unique(seqnames(cpg_data))) %>% 
#   .[1:24]
# 
# original_seqnames <- unique(as.character(seqnames(cpg_data)))
# 
# seqlevels(cpg_data) <- union(seqlevels(cpg_data), chr_names_NC)
# 
# name_map <- setNames(chr_names_NC, chr_names_chr)
# 
# # Apply the mapping to seqnames in cpg_data
# current_seqnames <- as.character(seqnames(cpg_data))
# new_names <- ifelse(current_seqnames %in% as.character(chr_names_chr), 
#                              as.character(name_map[current_seqnames]), 
#                              current_seqnames)
# cpg_hub <- GRanges(seqnames = new_names, ranges = ranges(cpg_data), strand = "*", mcols = mcols(cpg_data))
# 
# # add seqlevels of the query (SMRs)
# seqlevels(cpg_hub) <- union(seqlevels(cpg_hub), seqlevels(CpGs))
# 
# findOverlaps(CpGs, cpg_hub)
