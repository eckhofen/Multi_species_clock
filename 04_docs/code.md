# code inventory

## Repository-wide objects (files / datasets)

- **`000_data/000_metadata/color_palettes.RData`**
  - **Objects**: `color_compare`, `color_species`, `color_sex`
  - **Created by**: `001_scripts/00_data_preparation.R`

- **`000_data/000_metadata/metadata.RData`**
  - **Objects**: `meta_data`
  - **Created by**: `001_scripts/00_data_preparation.R`

- **`000_data/004_methyl_values/HS_*_meth_values(.Rdata/.csv)`**
  - **Objects**: `AC_meth_values`, `AS_meth_values`, `EH_meth_values`, `ZF_meth_values` / `ZF_meth_values_imputed`
  - **Created by**: `001_scripts/003_methylation_extraction.R`

- **`000_data/004_methyl_values/HS_*_methyl_sites(.Rdata/.csv)`**
  - **Objects**: `AC_methyl_sites`, `AS_methyl_sites`, `EH_methyl_sites`, `ZF_methyl_sites`
  - **Created by**: `001_scripts/003_methylation_extraction.R`

- **`000_data/004_methyl_values/HS_all_age.Rdata`**
  - **Objects**: `AC_age`, `AS_age`, `EH_age`, `EH_sex`, `ZF_age`, `meth_sites_names_tmp_AC`
  - **Created by**: `001_scripts/003_methylation_extraction.R`

- **`000_data/004_methyl_values/all_meth_values_long.RData`**
  - **Objects**: `all_meth_values_long`, `AC_meth_values_long`, `AS_meth_values_long`, `EH_meth_values_long`, `ZF_meth_values_long`
  - **Created by**: `001_scripts/003_methylation_extraction.R`

- **`000_data/005_correlation_data/cor_all.RData`**
  - **Objects**: `cor_all`
  - **Created by**: `001_scripts/004_correlation_testing.R`

- **`000_data/005_correlation_data/all_mix_cor_CpG_common.RData`**
  - **Objects**: `all_mix_cor_CpG_common`
  - **Created by**: `001_scripts/004_correlation_testing.R`

- **`000_data/006_model_creation/all_meth_values_selected.RData`**
  - **Objects**: `AC_meth_values_selected`, `AS_meth_values_selected`, `EH_meth_values_selected`, `ZF_meth_values_selected`, `all_meth_values_selected`
  - **Created by**: `001_scripts/004_correlation_testing.R`

- **`000_data/006_model_creation/all_meth_values_all_SMR.RData`**
  - **Objects**: `AC_meth_values_selected_all`, `AS_meth_values_selected_all`, `EH_meth_values_selected_all`, `ZF_meth_values_selected_all`, `all_meth_values_all_SMR`
  - **Created by**: `001_scripts/004_correlation_testing.R`

## Key functions by script

- **`001_scripts/001_CpGs_to_sequences.R`**
  - **`fix.seq(seq, rgenome, seq_width)`**: clamps extended GRanges to valid reference genome boundaries.
  - **`create.MethylPos(seqs, seqs_GR, methylsites, name = "CpGs")`**: returns per-sequence CpG positions and counts based on overlaps.

- **`001_scripts/002_rgHS_conserved_sequences.R`**
  - **`find.Overlap(...)`**: iteratively finds overlapping aligned reads across multiple `GAlignments`/`GRanges` inputs and returns overlap subsets per input.

- **`001_scripts/003_methylation_extraction.R`**
  - **`get.methyl.sites(seqs_aligned, species = "undefined", SMRs = "undefined")`**: maps original CpG positions into aligned coordinates using CIGAR parsing; assigns CpGs to SMRs.

- **`001_scripts/01_data_overview.R`**
  - **`theme_custom()`**: project ggplot theme.
  - **`plot_overview(p)`**: adds jitter/violin/mean + sample counts to a base ggplot.

- **`001_scripts/004_correlation_testing.R`**
  - **`cor.test.age(methyl_values, age, SMR = "not_defined", species = "undefined", method = "pearson")`**: computes correlation for each CpG with age.
  - **`cor.test.age.filter(input, p_value = 0.05)`**: annotates a correlation table with `significant`.
  - **`select.max.cor(cor_tibble, filter_significant = FALSE)`**: selects the strongest (absolute) correlation per SMR.

- **`001_scripts/005a_statistical_models_REL.R`** / **`005b_statistical_models.R`**
  - **`evaluate.model(...)`**: unified evaluation (R/MSE/MAE + plots + AE distributions + t-test) for regression models.

## Key plots and where they are created

- **`figures/01_fig_overview_age.PDF`**: `01_data_overview.R`
- **`figures/01_fig_overview_age_rel.PDF`**: `01_data_overview.R`
- **`figures/01_fig_overview_combined.PDF`**: `01_data_overview.R`

- **`002_plots/003_SMR_position_hist_all.pdf`**: `003_methylation_extraction.R`
- **`002_plots/003_SMR_position_lines_all.pdf`**: `003_methylation_extraction.R`
- **`002_plots/003_SMR_position.pdf`**: `003_methylation_extraction.R`
- **`002_plots/003_PCA_all.pdf`**: `003_methylation_extraction.R`

- **`002_plots/005_*`** model evaluation figures: `005a_statistical_models_REL.R` and `005b_statistical_models.R`
- **`002_plots/006_*`** comparison / LOSO plots: `006_data_comparison.R`
- **`002_plots/007_*`** region/gene annotation plots: `007_gene_analysis.R`
