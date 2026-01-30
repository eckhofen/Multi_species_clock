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
annotation_folder_new <- paste0(intermediate_folder, "08_annotation/")
annotation_folder <- annotation_folder_new
dir.create(annotation_folder, recursive = TRUE, showWarnings = FALSE)

annotation_gff_path <- paste0(annotation_folder, "GRCh38_p14_annotation.gff")
annotation_ucsc_cpg_path <- paste0(annotation_folder, "UCSC_CpG_islands.gff")

required_inputs <- c(
  annotation_gff_path,
  annotation_ucsc_cpg_path,
  paste0(intermediate_folder, "05_correlation_data/cor_all.RData"),
  paste0(intermediate_folder, "05_correlation_data/all_mix_cor_CpG_common.RData")
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  message(
    "Missing required inputs for gene analysis.\n",
    "Missing files:\n- ",
    paste(missing_inputs, collapse = "\n- "),
    "\nSkipping gene analysis."
  )
  quit(status = 0)
}

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
CpG_isl <- read.table(annotation_ucsc_cpg_path, sep = "\t", header = TRUE)

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
if (file.exists(path_methylsites_new)) {
  load(path_methylsites_new)
} else {
  meth_sites_inputs <- c(
    paste0(intermediate_folder, "04_methyl_values/HS_AC_methyl_sites.Rdata"),
    paste0(intermediate_folder, "04_methyl_values/HS_AS_methyl_sites.Rdata"),
    paste0(intermediate_folder, "04_methyl_values/HS_EH_methyl_sites.Rdata"),
    paste0(intermediate_folder, "04_methyl_values/HS_ZF_methyl_sites.Rdata")
  )

  missing_meth_sites <- meth_sites_inputs[!file.exists(meth_sites_inputs)]
  if (length(missing_meth_sites) > 0) {
    message(
      "Missing methyl site inputs for gene analysis.\n",
      "Missing files:\n- ",
      paste(missing_meth_sites, collapse = "\n- "),
      "\nSkipping gene analysis."
    )
    quit(status = 0)
  }

  load(meth_sites_inputs[1])
  load(meth_sites_inputs[2])
  load(meth_sites_inputs[3])
  load(meth_sites_inputs[4])

  methyl_sites_combined <- dplyr::bind_rows(
    if (exists("AC_methyl_sites")) dplyr::mutate(AC_methyl_sites, species = "AC") else NULL,
    if (exists("AS_methyl_sites")) dplyr::mutate(AS_methyl_sites, species = "AS") else NULL,
    if (exists("EH_methyl_sites")) dplyr::mutate(EH_methyl_sites, species = "EH") else NULL,
    if (exists("ZF_methyl_sites")) dplyr::mutate(ZF_methyl_sites, species = "ZF") else NULL
  )
}

 if (!exists("methyl_sites_combined")) {
   message("methyl_sites_combined is not available after loading inputs. Skipping gene analysis.")
   quit(status = 0)
 }

if (("SMR" %in% colnames(methyl_sites_combined)) && !("SCMR" %in% colnames(methyl_sites_combined))) {
  methyl_sites_combined <- dplyr::rename(methyl_sites_combined, SCMR = SMR)
}

 if (!("Site" %in% colnames(methyl_sites_combined))) {
   if (!("pos_rgenome" %in% colnames(methyl_sites_combined))) {
     message("methyl_sites_combined is missing required column pos_rgenome to construct Site. Skipping gene analysis.")
     quit(status = 0)
   }

   if (!("Chr" %in% colnames(methyl_sites_combined)) || !("species" %in% colnames(methyl_sites_combined))) {
     message("methyl_sites_combined is missing required columns Chr and/or species to construct Site. Skipping gene analysis.")
     quit(status = 0)
   }

   methyl_sites_combined <- methyl_sites_combined %>%
     dplyr::mutate(
       Site = dplyr::case_when(
         species == "AC" & grepl("_[0-9]+$", Chr) ~ paste0("X", gsub("^.*_", "", Chr), ".", pos_rgenome),
         species == "AS" ~ paste0(Chr, "-", pos_rgenome),
         species == "EH" ~ paste0(gsub("^EH_", "", Chr), ".", pos_rgenome),
         species == "ZF" ~ paste0(gsub("^ZF_", "", Chr), ":", pos_rgenome),
         TRUE ~ NA_character_
       )
     )
 }

if (!("chr_align" %in% colnames(methyl_sites_combined)) || !("pos_align" %in% colnames(methyl_sites_combined))) {
  message("methyl_sites_combined is missing required columns chr_align and/or pos_align. Skipping gene analysis.")
  quit(status = 0)
}

# transforming df into GRanges 
CpGs <- GRanges(seqnames = methyl_sites_combined$chr_align, ranges = IRanges(methyl_sites_combined$pos_align, methyl_sites_combined$pos_align, 1))

# load correlation data for all CpGs
path_cor_all_new <- paste0(intermediate_folder, "05_correlation_data/cor_all.RData")
load(path_cor_all_new)
cor_all

# adding metadata and correlation test results as meta column
# (join by keys instead of cbind to avoid row-count mismatches)
join_keys_candidates <- list(
  c("Site", "species"),
  c("Site")
)

join_keys <- character(0)
for (k in join_keys_candidates) {
  if (all(k %in% colnames(methyl_sites_combined)) && all(k %in% colnames(cor_all))) {
    join_keys <- k
    break
  }
}

if (length(join_keys) == 0) {
  message(
    "Could not find common join keys between methyl_sites_combined and cor_all.\n",
    "Available columns in methyl_sites_combined: ", paste(colnames(methyl_sites_combined), collapse = ", "), "\n",
    "Available columns in cor_all: ", paste(colnames(cor_all), collapse = ", "), "\n",
    "Skipping gene analysis."
  )
  quit(status = 0)
}

methyl_sites_joined <- dplyr::left_join(
  methyl_sites_combined,
  cor_all,
  by = join_keys,
  suffix = c("", "_cor")
)

if (nrow(methyl_sites_joined) != length(CpGs)) {
  message(
    "Joined methylation/correlation table row count does not match CpGs length (",
    nrow(methyl_sites_joined), " vs ", length(CpGs), "). Skipping gene analysis."
  )
  quit(status = 0)
}

mcols(CpGs) <- S4Vectors::DataFrame(
  methyl_sites_joined %>% dplyr::select(-chr_align, -pos_align),
  check.names = FALSE
)
CpGs

## selecting final CpGs
# loading selected CpGs
path_sel_cpg_new <- paste0(intermediate_folder, "05_correlation_data/all_mix_cor_CpG_common.RData")
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
region_counts <- function(x) {
  tab <- table(unlist(x))
  if (length(tab) == 0) {
    return(tibble::tibble(regions = character(), Freq = integer()))
  }
  tibble::tibble(regions = names(tab), Freq = as.integer(tab))
}

CpG_regions_pos <- CpG_regions_uni[which(CpGs$Correlation > 0)] %>% region_counts()
CpG_regions_neg <- CpG_regions_uni[which(CpGs$Correlation < 0)] %>% region_counts()

CpG_regions_df <- dplyr::full_join(
  CpG_regions_pos %>% dplyr::rename(pos_cor = Freq),
  CpG_regions_neg %>% dplyr::rename(neg_cor = Freq),
  by = "regions"
) %>%
  dplyr::mutate(
    pos_cor = tidyr::replace_na(pos_cor, 0L),
    neg_cor = tidyr::replace_na(neg_cor, 0L),
    check = pos_cor + neg_cor
  ) %>%
  .[.$check > 0, ] %>%
  .[!.$regions %in% "mRNA", ]

if (nrow(CpG_regions_df) == 0) {
  message("No CpG regions could be summarised from overlaps/correlation direction with the provided data. Skipping downstream gene analysis plots.")
  quit(status = 0)
}

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
  ylim(0, 1.1*max(CpG_regions_df_long$occurrences, na.rm = TRUE))
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
  ylim(0, 1.1*max(CpG_regions_all_df$CpGs, na.rm = TRUE))
p_CpGs_combined

ggsave(paste0(results_folder, "07_regions_all_CpGs_combined.pdf"), p_CpGs_combined, width = 8, height = 6)
 
# adding species of CpGs as dataframe
species_region_counts <- function(sel_species) {
  region_counts(CpG_regions_uni[which(CpGs$species == sel_species)]) %>%
    dplyr::rename(!!sel_species := Freq)
}

CpG_regions_df_species <- dplyr::full_join(species_region_counts("AC"), species_region_counts("AS"), by = "regions") %>%
  dplyr::full_join(species_region_counts("EH"), by = "regions") %>%
  dplyr::full_join(species_region_counts("ZF"), by = "regions") %>%
  dplyr::mutate(
    AC = tidyr::replace_na(AC, 0L),
    AS = tidyr::replace_na(AS, 0L),
    EH = tidyr::replace_na(EH, 0L),
    ZF = tidyr::replace_na(ZF, 0L)
  ) %>%
  .[!.$regions %in% "mRNA", ] %>%
  dplyr::mutate(check = AC + AS + EH + ZF) %>%
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
  ylim(0, ifelse(all(is.na(CpG_regions_df_species_long$occurrences)), 1, 1.1*max(CpG_regions_df_species_long$occurrences, na.rm = TRUE)))
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
gene_levels <- sort(unique(gene_hits_sel))
genes_sel_count <- tibble::tibble(
  genes = gene_levels,
  pos_cor = as.integer(table(factor(gene_hits_sel_pos, levels = gene_levels))),
  neg_cor = as.integer(table(factor(gene_hits_sel_neg, levels = gene_levels)))
)
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
  scale_fill_manual(values = rev(color_compare)) +
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
CpG_sel_regions_pos <- CpG_sel_regions_uni[which(CpGs_selected$Correlation > 0)] %>% region_counts()
CpG_sel_regions_neg <- CpG_sel_regions_uni[which(CpGs_selected$Correlation < 0)] %>% region_counts()

CpG_sel_regions_df <- dplyr::full_join(
  CpG_sel_regions_pos %>% dplyr::rename(pos_cor = Freq),
  CpG_sel_regions_neg %>% dplyr::rename(neg_cor = Freq),
  by = "regions"
) %>%
  dplyr::mutate(
    pos_cor = tidyr::replace_na(pos_cor, 0L),
    neg_cor = tidyr::replace_na(neg_cor, 0L),
    check = pos_cor + neg_cor
  ) %>%
  .[.$check > 0, ] %>%
  .[!.$regions %in% "mRNA", ]

if (nrow(CpG_sel_regions_df) == 0) {
  message("No selected-CpG regions could be summarised with the provided data. Skipping selected-CpG region plots.")
  quit(status = 0)
}

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
  ylim(0, 1.1*max(CpG_sel_regions_df_long$occurrences, na.rm = TRUE))
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
  ylim(0, 1.1*max(CpG_sel_regions_comb_df$CpGs, na.rm = TRUE))
p_CpGs_sel_combined

ggsave(paste0(results_folder, "07_regions_sel_CpGs_combined.pdf"), p_CpGs_sel_combined, width = 8, height = 6)

# adding species of CpGs as dataframe
species_region_counts_sel <- function(sel_species) {
  region_counts(CpG_sel_regions_uni[which(CpGs_selected$species == sel_species)]) %>%
    dplyr::rename(!!sel_species := Freq)
}

CpG_sel_regions_df_species <- dplyr::full_join(species_region_counts_sel("AC"), species_region_counts_sel("AS"), by = "regions") %>%
  dplyr::full_join(species_region_counts_sel("EH"), by = "regions") %>%
  dplyr::full_join(species_region_counts_sel("ZF"), by = "regions") %>%
  dplyr::mutate(
    AC = tidyr::replace_na(AC, 0L),
    AS = tidyr::replace_na(AS, 0L),
    EH = tidyr::replace_na(EH, 0L),
    ZF = tidyr::replace_na(ZF, 0L)
  ) %>%
  .[!.$regions %in% "mRNA", ] %>%
  dplyr::mutate(check = AC + AS + EH + ZF) %>%
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
  ylim(0, 1.1*max(CpG_sel_regions_df_species_long$occurrences, na.rm = TRUE))
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
