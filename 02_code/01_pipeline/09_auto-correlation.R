# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: SMR co-correlation analysis
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026 02

#### Libraries ####
library(patchwork)
library(ggrepel)
library(ggnewscale)
library(vegan)
library(reshape2)
library(Hmisc)
library(data.table)
library(tidyverse)

#### Setting paths ####
figures_folder <- "03_results/01_figures/"
save_folder <- "01_data/03_intermediate/09_co-methylation/"
dir.create(save_folder, recursive = TRUE, showWarnings = FALSE)

#### Loading data ####
# annotated CpGs
load("01_data/03_intermediate/08_annotation/cpgs_SMR_annotated_unique.RData")

# corrlation data & metadata
load("01_data/03_intermediate/05_correlation_data/cor_all.RData")

# methylation values
load("01_data/03_intermediate/06_model_creation/all_meth_values_selected.RData")
load("01_data/03_intermediate/04_methyl_values/all_meth_values.Rdata")

# colors
load("01_data/04_metadata/color_palettes.RData")

# ggplot theme settings
source("02_code/02_helpers/plot_style.R")

# loading whole methyl data
xyz <- load("01_data/01_raw_private/AC_methyldata.Rdata")
assign("AC_meth_data", get(xyz))

xyz <- load("01_data/01_raw_private/AS_methyldata.Rdata")
assign("AS_meth_data", get(xyz))

xyz <- load("01_data/01_raw_private/EH_methyldata.Rdata")
assign("EH_meth_data", get(xyz))

xyz <- load("01_data/01_raw_private/ZF_methyldata.Rdata")
assign("ZF_meth_data", get(xyz))

#### Helper functions ####
# for metadata extractions 

get_meta_AC <- function(meth_data) {
  cpgs <- colnames(meth_data)
  clean_cpgs <- sub("^X", "", cpgs)
  
  data.table(
    Site = cpgs,
    chr = sub("\\..*$", "", clean_cpgs),
    pos = as.numeric(sub("^.*\\.", "", clean_cpgs))
  )
}

get_meta_AS <- function(meth_data) {
  cpgs <- colnames(meth_data)
  
  suppressWarnings(data.table(
    Site = cpgs,
    chr = sub("-[0-9]+$", "", cpgs),
    pos = as.numeric(sub("^.*-", "", cpgs))
  ))
}

get_meta_EH <- function(meth_data) {
  cpgs <- colnames(meth_data)
  
  data.table(
    Site = cpgs,
    chr = sub("\\..*$", "", cpgs),
    pos = as.numeric(sub("^.*\\.", "", cpgs))
  )
}

get_meta_ZF <- function(meth_data) {
  cpgs <- colnames(meth_data)
  
  data.table(
    Site = cpgs,
    chr = sub(":.*$", "", cpgs),
    pos = as.numeric(sub("^.*:", "", cpgs)) 
  )
}

# Autocorrelation calculation

calc_distance_autocorrelation <- function(meth_data, 
                                          meta_dt, 
                                          group_var = "chr", 
                                          max_bp_dist = 500, 
                                          window_size = 10, 
                                          mem_threshold = 200000) {
  
  cat("Cleaning metadata and finding valid pairs...\n")
  
  # Standardize Metadata
  meta_dt <- as.data.table(meta_dt) %>%
    rename_with(tolower)
  
  group_var <- str_to_lower(group_var)
  setnames(meta_dt, group_var, "grp")
  
  # Filter and set keys for non-equi join
  meta_dt <- meta_dt[!is.na(get("grp")) & !is.na(pos)]
  setorder(meta_dt, grp, pos)
  
  # Find pairs within max_bp_dist within the same group (chr or SMR)
  meta_dt[, pos_max := pos + max_bp_dist]
  
  pairs_dt <- meta_dt[meta_dt, 
                      on = .(grp, pos > pos, pos <= pos_max), 
                      nomatch = NULL, 
                      .(Site1 = x.site, Site2 = i.site, dist = i.pos - x.pos)]
  
  pairs_dt[, dist_bin := floor(dist / window_size) * window_size]
  
  cat(sprintf("Found %s pairs. Preparing methylation data...\n", format(nrow(pairs_dt), big.mark=",")))
  
  # Prepare Methylation Matrix
  meth_matrix <- as.matrix(meth_data)
  suppressWarnings(mode(meth_matrix) <- "numeric")
  
  sample_names <- rownames(meth_matrix) %||% as.character(1:nrow(meth_matrix))
  cpg_names <- colnames(meth_matrix)
  
  # Filter pairs to only those present in the methylation matrix
  pairs_dt <- pairs_dt[Site1 %in% cpg_names & Site2 %in% cpg_names]
  
  # Determine Mode based on number of positions
  n_positions <- nrow(meta_dt)
  
  if (n_positions > mem_threshold) {
    cat(sprintf("Positions (%d) > threshold. Using LOW MEMORY mode...\n", n_positions))
    
    n_samples <- length(sample_names)
    res_list <- vector("list", n_samples)
    
    for (i in 1:n_samples) {
      if (i %% 5 == 0 || i == 1) cat(sprintf("   -> Sample %d of %d\n", i, n_samples))
      
      meth_vec <- meth_matrix[i, ]
      valid_idx <- !is.na(meth_vec)
      samp_dt <- data.table(Site = cpg_names[valid_idx], Meth = meth_vec[valid_idx])
      
      tmp <- pairs_dt[samp_dt, on = .(Site1 = Site), nomatch = NULL, .(Site2, dist, dist_bin, Meth1 = i.Meth)]
      tmp <- tmp[samp_dt, on = .(Site2 = Site), nomatch = NULL, .(dist, dist_bin, Meth1, Meth2 = i.Meth)]
      
      res_list[[i]] <- tmp[, .(
        Sample = sample_names[i],
        mean_actual_dist = mean(dist),
        n_pairs = .N,
        pearson_r = if(.N >= 3) cor(Meth1, Meth2, method = "pearson") else NA_real_
      ), by = dist_bin]
    }
    final_res <- rbindlist(res_list)
    
  } else {
    cat(sprintf("Positions (%d) <= threshold. Using FAST mode...\n", n_positions))
    
    # Melt the whole matrix to long format
    meth_dt <- as.data.table(meth_matrix, keep.rownames = "Sample")
    meth_long <- melt(meth_dt, id.vars = "Sample", variable.name = "Site", value.name = "Meth", na.rm = TRUE)
    setkey(meth_long, Site)
    
    # Join Site 1 values
    meth_pairs <- merge(pairs_dt, meth_long, by.x = "Site1", by.y = "Site", allow.cartesian = TRUE)
    setnames(meth_pairs, "Meth", "Meth1")
    
    # Join Site 2 values (matching on Site AND Sample)
    meth_pairs <- merge(meth_pairs, meth_long, by.x = c("Site2", "Sample"), by.y = c("Site", "Sample"))
    setnames(meth_pairs, "Meth", "Meth2")
    
    # Calculate all correlations at once
    final_res <- meth_pairs[, .(
      mean_actual_dist = mean(dist),
      n_pairs = .N,
      pearson_r = if(.N >= 3) cor(Meth1, Meth2, method = "pearson") else NA_real_
    ), by = .(Sample, dist_bin)]
  }
  
  cat("Finalizing results...\n")
  setorder(final_res, Sample, dist_bin)
  return(as_tibble(final_res))
}

# Helper function for cor plot 
plot_cpg_correlation <- function(data, 
                                 y_col, 
                                 y_label = "Pearson Correlation (R)", 
                                 y_limits = c(0, 1), 
                                 colors = color_species, 
                                 show_samples = FALSE, 
                                 facet_species = FALSE, 
                                 loess = TRUE) {
  y_breaks = seq(min(y_limits), max(y_limits), by = 0.2)

  p <- ggplot(data, aes(x = abs(mean_actual_dist), y = .data[[y_col]], group = Sample)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "CpG Distance (bp)", y = y_label, color = "Species") +
    scale_x_continuous(breaks = seq(0, 300, by = 50), limits = c(0, 300)) +
    scale_y_continuous(breaks = y_breaks) +
    scale_color_manual(values = colors)

  if (show_samples) {
    p <- p + geom_line(aes(color = species), alpha = 0.05, show.legend = FALSE)
  }
  if (loess) {
  p <- p + geom_smooth(aes(group = species, color = species), 
                       method = "loess", 
                       se = FALSE, 
                       span = 0.2)
  }
  if (facet_species) {
    p <- p + facet_wrap(~ species, ncol = 1)
  }

  p <- p + coord_cartesian(ylim = y_limits)

  return(p)
}

#### Analysis ####
# SMR data
# Prepare meta_data for grouping
smr_meta_data <- cor_all %>% select(Site, species, SMR, pos_align) %>% rename(pos = pos_align)

# Calculate autocorrelation
AC_smr_cor <- calc_distance_autocorrelation(AC_meth_values, smr_meta_data, group_var = "SMR")
AS_smr_cor <- calc_distance_autocorrelation(AS_meth_values, smr_meta_data, group_var = "SMR")
EH_smr_cor <- calc_distance_autocorrelation(EH_meth_values, smr_meta_data, group_var = "SMR")
ZF_smr_cor <- calc_distance_autocorrelation(ZF_meth_values, smr_meta_data, group_var = "SMR")

AC_smr_cor$species <- "AC"
AS_smr_cor$species <- "AS"
EH_smr_cor$species <- "EH"
ZF_smr_cor$species <- "ZF"

combined_smr_cor <- rbind(AC_smr_cor, AS_smr_cor, EH_smr_cor, ZF_smr_cor)

# All data 
AC_meta    <- get_meta_AC(AC_meth_data)
AC_all_cor <- calc_distance_autocorrelation(AC_meth_data, AC_meta)
AC_all_cor$species <- "AC"

AS_meta    <- get_meta_AS(AS_meth_data)
AS_all_cor <- calc_distance_autocorrelation(AS_meth_data, AS_meta)
AS_all_cor$species <- "AS"

EH_meta    <- get_meta_EH(EH_meth_data)
EH_all_cor <- calc_distance_autocorrelation(EH_meth_data, EH_meta)
EH_all_cor$species <- "EH"

ZF_meta    <- get_meta_ZF(ZF_meth_data)
ZF_all_cor <- calc_distance_autocorrelation(ZF_meth_data, ZF_meta)
ZF_all_cor$species <- "ZF"

combined_all_cor <- rbind(AC_all_cor, AS_all_cor, EH_all_cor, ZF_all_cor)

# Normalize r values 
combined_all_cor_normed <- combined_all_cor %>%
  group_by(species) %>%
  mutate(pearson_r_norm = pearson_r / max(abs(pearson_r), na.rm = TRUE))

p_all_sample <- plot_cpg_correlation(combined_all_cor, "pearson_r", "Pearson Correlation (r)", show_samples = TRUE, facet_species = TRUE)
p_all_loess <- plot_cpg_correlation(combined_all_cor, "pearson_r", "Pearson Correlation (r)")
p_all_loess_normed <- plot_cpg_correlation(combined_all_cor_normed, "pearson_r_norm", "Pearson Correlation normalized (r)")

save_plot("03_results/01_figures/09_all_sample_auto_correlation", p_all_sample, width = 6, height = 8)
save_plot("03_results/01_figures/09_all_loess_auto_correlation", p_all_loess, width = 6, height = 4)
save_plot("03_results/01_figures/09_all_loess_auto_correlation_normed", p_all_loess_normed, width = 6, height = 4)

# SMR results
# Normalize as well
combined_smr_cor_normed <- combined_smr_cor %>%
  group_by(species) %>%
  mutate(pearson_r_norm = pearson_r / max(abs(pearson_r), na.rm = TRUE))

p_smr_sample <- plot_cpg_correlation(combined_smr_cor_normed, "pearson_r", "Pearson Correlation (r)", show_samples = TRUE, facet_species = TRUE, y_limits = c(-0.4, 1))
p_smr_loess <- plot_cpg_correlation(combined_smr_cor_normed, "pearson_r", "Pearson Correlation (r)", y_limits = c(-0.4, 1))
p_smr_loess_sample <- plot_cpg_correlation(combined_smr_cor_normed, "pearson_r", "Pearson Correlation (r)", y_limits = c(-0.4, 1), show_samples = TRUE)
p_smr_loess_normed <- plot_cpg_correlation(combined_smr_cor_normed, "pearson_r_norm", "Pearson Correlation normalized (r)", y_limits = c(-0.4, 1))

save_plot("03_results/01_figures/09_smr_sample_auto_correlation", p_smr_sample, width = 6, height = 8)
save_plot("03_results/01_figures/09_smr_loess_auto_correlation", p_smr_loess, width = 6, height = 4)
save_plot("03_results/01_figures/09_smr_loess_sample_auto_correlation", p_smr_loess_sample, width = 6, height = 4)
save_plot("03_results/01_figures/09_smr_loess_auto_correlation_normed", p_smr_loess_normed, width = 6, height = 4)

# Combined plots in one
combined_smr_cor$dataset <- "SMR"
combined_all_cor$dataset <- "All"

combined_all <- rbind(combined_smr_cor, combined_all_cor)

p_combined_all <- plot_cpg_correlation(combined_all, "pearson_r", "Pearson Correlation (r)", y_limits = c(-0.4, 1), loess = FALSE) + 
  geom_smooth(aes(group = interaction(species, dataset), color = species, linetype = dataset), method = "loess", se = FALSE, span = 0.2) +
  labs(linetype = "CpGs")

save_plot("03_results/01_figures/09_combined_all.png", p_combined_all, width = 6, height = 6)

# Combined panel plots
p_comparison_smr_all <- (p_smr_loess + theme(legend.position = c(0.89, 0.80))) + (p_all_loess + theme(legend.position = "none")) + 
  plot_layout(nrow = 2, axes = "collect_x") + 
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 18, face = "bold")) &
  labs(color = "")

save_plot("03_results/01_figures/09_comparison_smr_all.png", p_comparison_smr_all, width = 6, height = 6)

# Normalized
p_comparison_smr_all_normed <- p_smr_loess_normed + p_all_loess_normed + 
  plot_layout(nrow = 2, guides = "collect", axes = "collect_x") + 
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

save_plot("03_results/01_figures/09_comparison_smr_all_normed.png", p_comparison_smr_all_normed, width = 6, height = 6)

# Combined plot for figure 2/3
# loading plot saved by 08_genomic_context_analysis.R
load("01_data/03_intermediate/08_annotation/08_p_combined_count.RData")

p_combined_cor_context <- (p_comparison_smr_all | p_combined_count) +
  plot_layout(axes = "collect_x") + 
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 18, face = "bold"))

save_plot("03_results/01_figures/09_combined_cor_context.png", p_combined_cor_context, width = 10, height = 8)
