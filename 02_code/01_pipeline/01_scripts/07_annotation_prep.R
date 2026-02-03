# Metadata ----------------------------------------------------------------
# Project: Piscine Multispecies Epigenetic Clock
# Description: SCMR gene analysis data preparation
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
load("01_data/03_intermediate/08_annotation/annotation.RData")
load("01_data/03_intermediate/08_annotation/SCMRs_and_CpGs.RData")

# importing & building annotations
annots <- c("hg38_cpgs", "hg38_basicgenes", "hg38_lncrna_gencode")
annotations <- build_annotations(genome = "hg38", annotations = annots)

# setting chromosome style to UCSC to match with annotations
seqlevelsStyle(cpgs_scmr) <- "UCSC"
seqlevelsStyle(cpgs_selected) <- "UCSC"

# removing non-standard chromosomes
cpgs_scmr <- keepStandardChromosomes(cpgs_scmr, species = "Homo_sapiens", pruning.mode = "coarse")
cpgs_selected <- keepStandardChromosomes(cpgs_selected, species = "Homo_sapiens", pruning.mode = "coarse")

# annotate CpGs
cpgs_scmr_annotated <- annotate_regions(cpgs_scmr, annotations, ignore.strand = TRUE)
cpgs_selected_annotated <- annotate_regions(cpgs_selected, annotations, ignore.strand = TRUE)

# check annotations
summarize_annotations(cpgs_scmr_annotated)
summarize_annotations(cpgs_selected_annotated)

# get annot. types
annot_types <- unique(cpgs_scmr_annotated$annot.type)
annot_names <- unique(cpgs_scmr_annotated$annot.name)

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
scmr_counts <- get_counts(cpgs_scmr_annotated)
selected_counts <- get_counts(cpgs_selected_annotated)

# save
save(scmr_counts, file = "01_data/03_intermediate/08_annotation/SCMR_counts.RData")
save(selected_counts, file = "01_data/03_intermediate/08_annotation/Selected_counts.RData")

# get affected genes and lncRNAs helper function
get_genes <- function(annotated_regions) {
    df <- as.data.frame(annotated_regions) %>%
      group_by(species)

    # get genes
    genes <- df %>%
      dplyr::select(species, annot.gene_id, annot.symbol, SCMR, seq_number, annot.type) %>%
      mutate(type = ifelse(grepl("gene", annot.type), "gene", "lncrna")) %>%
      dplyr::select(-annot.type) %>%
      distinct() %>%
      na.omit()

    return(genes)
}

# get genes
genes_scmr <- get_genes(cpgs_scmr_annotated)
genes_selected <- get_genes(cpgs_selected_annotated)

# save genes
save(genes_scmr, file = "01_data/03_intermediate/08_annotation/SCMR_genes.RData")
save(genes_selected, file = "01_data/03_intermediate/08_annotation/Selected_genes.RData")

write.csv(genes_scmr, file = "01_data/03_intermediate/08_annotation/gene_list_scmr.csv")
write.csv(genes_selected, file = "01_data/03_intermediate/08_annotation/gene_list_selected.csv")

