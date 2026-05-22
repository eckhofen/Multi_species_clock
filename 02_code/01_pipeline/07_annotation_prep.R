# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: SMR gene analysis data preparation
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2025 06

#### Preparation ####
library(tidyverse)
library(AnnotationHub)
library(GenomicFeatures)
library(GenomeInfoDb)
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

#### Loading data ####
load("01_data/03_intermediate/05_correlation_data/cor_all.RData")

#### Preparing data ####
# select CpGs
cpgs_SMR <- GRanges(seqnames = cor_all$chr_align, 
                     ranges = IRanges(cor_all$pos_align, cor_all$pos_align, 1))

mcols(cpgs_SMR) <- cor_all

cpgs_selected <- cpgs_SMR[cpgs_SMR$selected == TRUE,]

# changing chromosome style to NCBI
chrom_info_human <- getChromInfoFromNCBI("GRCh38")
rename_vector <- chrom_info_human$UCSCStyleName
names(rename_vector) <- chrom_info_human$RefSeqAccn

# rename chromosomes
cpgs_selected <- renameSeqlevels(cpgs_selected, rename_vector)
cpgs_SMR <- renameSeqlevels(cpgs_SMR, rename_vector)

# importing & building annotations
annots <- c("hg38_cpgs", "hg38_basicgenes", "hg38_lncrna_gencode")
annotations <- build_annotations(genome = "hg38", annotations = annots)

# setting chromosome style to UCSC to match with annotations
seqlevelsStyle(cpgs_SMR) <- "UCSC"
seqlevelsStyle(cpgs_selected) <- "UCSC"

# annotate CpGs
cpgs_SMR_annotated <- annotate_regions(cpgs_SMR, annotations, ignore.strand = TRUE)
cpgs_selected_annotated <- annotate_regions(cpgs_selected, annotations, ignore.strand = TRUE)

# save
save(cpgs_SMR_annotated, file = "01_data/03_intermediate/08_annotation/SMR_annotated.RData")
save(cpgs_selected_annotated, file = "01_data/03_intermediate/08_annotation/Selected_annotated.RData")

# check annotations
summarize_annotations(cpgs_SMR_annotated)
summarize_annotations(cpgs_selected_annotated)

# get annot. types
annot_types <- unique(cpgs_SMR_annotated$annot.type)
annot_names <- unique(cpgs_SMR_annotated$annot.name)

# helper function to get counts per species
get_counts <- function(cpgs_annotated) {
    # group by species
    species_list <- split(cpgs_annotated, cpgs_annotated$species)

    # count per annotation type
    counts_list <- lapply(species_list, function(x) {
      summarize_annotations(x)
    })
    # add species column
    for (species in names(counts_list)) {
        counts_list[[species]]$species <- species
    }
    
    # combine into tibble
    counts <- do.call(rbind, counts_list) %>%
        dplyr::select(species, everything()) %>%
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
            TRUE ~ annot.type))
    return(counts)
}

# select per species
SMR_counts <- get_counts(cpgs_SMR_annotated)
selected_counts <- get_counts(cpgs_selected_annotated)

# save
save(SMR_counts, file = "01_data/03_intermediate/08_annotation/SMR_counts.RData")
save(selected_counts, file = "01_data/03_intermediate/08_annotation/Selected_counts.RData")

# get affected genes and lncRNAs helper function
get_genes <- function(annotated_regions) {
    df <- as.data.frame(annotated_regions) %>%
      group_by(species)

    # get genes
    genes <- df %>%
      dplyr::select(Site, species, annot.gene_id, annot.symbol, SMR, seq_number, annot.type) %>%
      mutate(type = ifelse(grepl("gene", annot.type), "gene", "lncrna")) %>%
      dplyr::select(-annot.type) %>%
      distinct() %>%
      na.omit()

    return(genes)
}

# get genes
genes_SMR <- get_genes(cpgs_SMR_annotated)
genes_selected <- get_genes(cpgs_selected_annotated)

# save genes
save(genes_SMR, file = "01_data/03_intermediate/08_annotation/SMR_genes.RData")
save(genes_selected, file = "01_data/03_intermediate/08_annotation/Selected_genes.RData")

cpgs_meta_data <- cor_all %>% dplyr::select(Site, Correlation, P_value)

gene_list_SMR <- genes_SMR %>% 
  left_join(cpgs_meta_data, by = "Site") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

gene_list_selected <- genes_selected %>% 
  left_join(cpgs_meta_data, by = "Site") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

write.csv(gene_list_SMR, file = "03_results/02_tables/gene_list_SMR.csv")
write.csv(gene_list_selected, file = "03_results/02_tables/gene_list_selected.csv")

