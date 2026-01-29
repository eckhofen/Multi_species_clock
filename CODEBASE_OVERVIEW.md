# Codebase overview

## High-level purpose

This repository implements a comparative epigenetics workflow to identify **Shared Conserved Methylation Regions (SCMRs)** across fish species and use them to build a **multispecies piscine epigenetic clock**.

Conceptual pipeline (as implemented by the scripts):

1. **Sample metadata curation** (age, sex, species, max lifespan; relative age).
2. **(Upstream / not fully reproducible from this repo)** Expand CpG coordinates to sequences and align them to a reference genome (human used in legacy scripts).
3. Identify overlaps between aligned sequences to define **SCMRs**.
4. Extract CpG methylation values belonging to those SCMRs.
5. Test CpG–age correlations; select best CpG per SCMR per species and keep SCMRs with consistent cross-species direction.
6. Build regression models to predict age (relative or chronological).
7. Compare models and interpret results via genomic annotation / gene context.

## Repository structure

- **`01_data/`**
  - Inputs and derived artifacts (including external server outputs and intermediate data).
- **`02_code/`**
  - Refactored analysis code (pipeline scripts, helpers).
- **`03_results/`**
  - Result outputs (figures/tables/logs). Not tracked in git.
- **`04_docs/`**
  - Documentation.

- **Legacy code**
  - Older scripts and deprecated structures are kept under `02_code/99_archive/` (gitignored) and are not part of the active pipeline.

## External inputs (server-only / not reproducible from this repo)

The following steps depend on large private datasets and/or a server alignment workflow and are treated as **external inputs** to the locally reproducible pipeline:

- **Sequence expansion + metadata**: files copied into `01_data/02_external_server_outputs/01_sequences/` (e.g. `*_metadata_*1000bp.csv`).
- **Conserved/aligned overlap objects**: files copied into `01_data/02_external_server_outputs/02_conserved_seq/` (e.g. `HS_AC_AS_EH_ZF_overlaps.Rdata`).

Downstream scripts should start from these external inputs and write derived artifacts to `01_data/03_intermediate/` and figures to `03_results/`.

## Data artifacts (authoritative in-repo I/O)

### Metadata

- **Inputs**: `01_data/01_raw_private/*_age.csv` (AC/AS/EH/ZF)
- **Outputs**:
  - `01_data/04_metadata/color_palettes.RData`
  - `01_data/04_metadata/metadata.RData`

### SMR and methylation extraction

- **Outputs (key files in repo)**: `01_data/03_intermediate/`
  - `004_methyl_values/HS_*_meth_values(.Rdata/.csv)`
  - `004_methyl_values/HS_*_methyl_sites(.Rdata/.csv)`
  - `004_methyl_values/HS_all_age.Rdata`
  - `004_methyl_values/all_meth_values_long.RData`
  - `005_correlation_data/all_mix_cor_CpG_common.RData`

### Correlation + model creation

- `01_data/03_intermediate/005_correlation_data/cor_all.RData`
- `01_data/03_intermediate/005_correlation_data/all_mix_cor_CpG_common.RData`
- `01_data/03_intermediate/006_model_creation/all_meth_values_selected.RData`
- `01_data/03_intermediate/006_model_creation/all_meth_values_all_SMR.RData`

## Script-by-script overview

### `02_code/01_pipeline/01_scripts/00_data_preparation.R`

- **Goal**: define global constants (species max lifespan), build unified metadata, define project color palettes.
- **Reads**:
  - `01_data/01_raw_private/AC_age.csv`
  - `01_data/01_raw_private/AS_age.csv`
  - `01_data/01_raw_private/EH_age.csv`
  - `01_data/01_raw_private/ZF_age.csv`
- **Creates (objects)**:
  - `age_max_AC`, `age_max_AS`, `age_max_EH`, `age_max_ZF`
  - `color_species`, `color_species_common`, `color_species_sci`, `color_compare`, `color_sex`
  - `AC`, `AS`, `EH`, `ZF` (raw metadata tibbles)
  - `meta_data` (combined metadata with `age_max` and `age_rel`)
- **Writes**:
  - `01_data/04_metadata/color_palettes.RData`
  - `01_data/04_metadata/metadata.RData`

### `02_code/01_pipeline/01_scripts/01_data_overview.R`

- **Goal**: visualize sample age distributions (absolute and relative).
- **Reads**:
  - `01_data/04_metadata/color_palettes.RData`
  - `01_data/04_metadata/metadata.RData`
- **Creates (objects)**:
  - `n_per_species`
  - `theme_custom()`
  - `plot_overview(p)`
  - `fig_overview_age`, `fig_overview_age_rel`, `fig_overview_combined`
- **Writes (plots)**:
  - `03_results/01_figures/01_fig_overview_age.PDF`
  - `03_results/01_figures/01_fig_overview_age_rel.PDF`
  - `03_results/01_figures/01_fig_overview_combined.PDF`

### `02_code/99_archive/legacy_scripts_from_pipeline/000_age_metadata.R` (legacy)

- **Goal**: older standalone age distribution plots.
- **Notes**:
  - Uses `setwd("/Users/.../Repo_Multispecies_clock/..." )` (hard-coded).
  - Reads `000_data/000_metadata/metadata_age.RData` and expects objects like `AC_df` etc.
- **Writes (plots)**:
  - `002_plots/000_rel_age_distribution.pdf` (via `extension`)
  - `002_plots/000_age_distribution.pdf`
  - `002_plots/000_age_distribution_both.pdf`

### `02_code/99_archive/legacy_scripts_from_pipeline/001_CpGs_to_sequences.R` (upstream / external paths)

- **Goal**: expand CpG loci to fixed-length sequences (default 1000bp), write FASTA + metadata.
- **Reads**:
  - Multiple methyl position tables and reference genomes from external locations (`/powerplant/...`, `raw-data/...`).
- **Core objects/functions**:
  - `bp_ext`, `data_folder`, `save_folder`, `file_ext`
  - `fix.seq(seq, rgenome, seq_width)`
  - `create.MethylPos(seqs, seqs_GR, methylsites, name = "CpGs")`
  - Species GRanges: `EH`, `AC`, `AS`, `ZF` and methyl ranges `*_methyl`
  - Extracted sequences: `*_seq` (`DNAStringSet`)
  - Per-sequence CpG metadata: `*_metadata`
- **Writes**:
  - `000_data/001_sequences/*_CpG_1000bp.fasta`
  - `000_data/001_sequences/*_metadata_1000bp.csv` (naming varies per species)

### `02_code/99_archive/legacy_scripts_from_pipeline/002_rgHS_conserved_sequences.R` (upstream / external paths)

- **Goal**: determine conserved/overlapping aligned sequences across species based on BAM alignments.
- **Reads**:
  - BAM files from `results-data/bowtie2/*.bam`
  - metadata csv from `000_data/001_sequences/*_metadata_*.csv`
- **Creates**:
  - `HS_*_1000_bt2` (`GAlignments`)
  - `find.Overlap(...)`
  - `HS_overlap_seqs`
- **Writes**:
  - `000_data/002_conserved_seq/HS_AC_AS_EH_ZF_overlaps.Rdata`

### `02_code/01_pipeline/01_scripts/03_methylation_extraction.R`

- **Goal**: construct SMRs from overlaps, map CpG positions onto aligned coordinates, extract methylation matrices for CpGs within SMRs, and do initial QC/visualizations.
- **Reads**:
  - `01_data/02_external_server_outputs/02_conserved_seq/HS_AC_AS_EH_ZF_overlaps.Rdata`
  - Private methylation inputs under `01_data/01_raw_private/` (not tracked).
- **Creates (key objects)**:
  - SMR ranges: `HS_SMR_b` (reduced genomic ranges of overlaps)
  - `get.methyl.sites(seqs_aligned, species, SMRs)`
  - `*_methyl_sites` for AC/AS/EH/ZF
  - `methyl_sites_combined` and normalized version `methyl_sites_combined_nor`
  - Methyl matrices: `AC_meth_values`, `AS_meth_values`, `EH_meth_values`, `ZF_meth_values`
  - Imputed ZF matrix: `ZF_meth_values_imputed`
  - PCA objects: `PCA_AC`, `PCA_AS`, `PCA_EH`, `PCA_ZF` (+ plotting data frames)
  - Long-format methylation: `*_meth_values_long`, `all_meth_values_long`
- **Writes (data)**:
  - `01_data/03_intermediate/03_SCMR/HS_SCMR.RData`
  - `01_data/03_intermediate/008_annotation/methylsites_all.RData`
  - `01_data/03_intermediate/007_data_comparison/methyl_sites_combined_nor.Rdata`
  - `01_data/03_intermediate/004_methyl_values/HS_*_meth_values(.Rdata/.csv)`
  - `01_data/03_intermediate/004_methyl_values/HS_*_methyl_sites(.Rdata/.csv)`
  - `01_data/03_intermediate/004_methyl_values/HS_all_age.Rdata`
  - `01_data/03_intermediate/004_methyl_values/HS_ZF_meth_values_imputed(.Rdata/.csv)`
  - `01_data/03_intermediate/004_methyl_values/all_meth_values_long.RData`
- **Writes (plots)**:
  - `03_results/01_figures/03_SCMR_position_hist_all.pdf`
  - `03_results/01_figures/03_SCMR_position_lines_all.pdf`
  - `03_results/01_figures/03_SCMR_position.pdf`
  - `03_results/01_figures/03_PCA_all.pdf`

### `02_code/01_pipeline/01_scripts/04_correlation_testing.R`

- **Goal**: compute per-CpG correlation with age, filter by significance, select max-|r| CpG per SMR per species, and export selected feature matrices.
- **Reads**:
  - `01_data/03_intermediate/004_methyl_values/HS_*_meth_values.Rdata`
  - `01_data/03_intermediate/004_methyl_values/HS_*_methyl_sites.Rdata`
  - `01_data/03_intermediate/004_methyl_values/HS_all_age.Rdata`
- **Creates (objects/functions)**:
  - `cor.test.age(...)`
  - `cor.test.age.filter(...)`
  - `select.max.cor(...)`
  - `*_cor_age_filtered_pearson` and combined `cor_all`
  - Selected CpG lists: `all_pos_cor_CpG_common`, `all_neg_cor_CpG_common`, `all_mix_cor_CpG_common`
  - Feature matrices per species:
    - `*_meth_values_selected` (one CpG per SMR; SMRs with 4 species; mixed rule)
    - `*_meth_values_selected_all` (all pos+neg SMRs, with `_pos/_neg` labels)
  - Combined matrices:
    - `all_meth_values_selected`
    - `all_meth_values_all_SMR`
- **Writes**:
  - `01_data/03_intermediate/005_correlation_data/cor_all.RData`
  - `01_data/03_intermediate/005_correlation_data/all_mix_cor_CpG_common.RData`
  - `01_data/03_intermediate/006_model_creation/all_meth_values_selected.RData`
  - `01_data/03_intermediate/006_model_creation/all_meth_values_all_SMR.RData`

### `02_code/01_pipeline/01_scripts/05a_statistical_models_REL.R`

- **Goal**: build and evaluate predictive models for **relative age**.
- **Reads**:
  - `01_data/03_intermediate/006_model_creation/all_meth_values_selected.RData`
- **Key created objects**:
  - Train/test splits per species (`*_split`) and combined `meth_train`, `meth_test`
  - Evaluation helper: `evaluate.model(...)`
  - Models:
    - Elastic Net regression via `cv.glmnet` (normal and transformed response)
    - Multivariate linear regression (`lm`) (normal and transformed response; optional “sig-only” version)
    - Random forest (`randomForest`) (normal and transformed)
    - SVM regression (`e1071::svm`) (eps and nu)
  - AE (absolute error) aggregation tables (used for later comparison)
- **Writes (plots)**: a variety of `03_results/01_figures/05_*` model evaluation PDFs.

### `02_code/01_pipeline/01_scripts/05b_statistical_models.R`

- **Goal**: build and evaluate predictive models for **chronological age** (same modeling scaffold as `005a`, different target transformation and plot axes).
- **Reads**:
  - `01_data/03_intermediate/006_model_creation/all_meth_values_selected.RData`
- **Writes (plots)**: a variety of `03_results/01_figures/05_*` model evaluation PDFs.

### `02_code/01_pipeline/01_scripts/06_data_comparison.R`

- **Goal**: compare model coefficients / SMR importance between relative vs chronological models; compile additional comparison plots (including LOSO blocks).
- **Reads**:
  - `01_data/03_intermediate/007_data_comparison/mlm_age_summary.Rdata`
  - `01_data/03_intermediate/007_data_comparison/mlm_rel_age_summary.Rdata`
  - `01_data/03_intermediate/007_data_comparison/methyl_sites_combined_nor.Rdata`
  - `01_data/03_intermediate/005_correlation_data/all_mix_cor_CpG_common.RData`
  - Various LOSO result `.Rdata` (expected to exist from prior runs)
- **Creates**:
  - `df_SMR_comparison` / `df_SMR_comparison_all` (regression coefficient summaries)
- **Writes (plots)**: `03_results/01_figures/06_*` (comparison / LOSO plots)

### `02_code/01_pipeline/01_scripts/07_gene_analysis.R`

- **Goal**: annotate CpGs/SMRs on the human reference genome; summarize impacted genomic regions and genes; compare with GenAge gene lists.
- **Reads**:
  - `01_data/03_intermediate/008_annotation/GRCh38_p14_annotation.gff`
  - `01_data/03_intermediate/008_annotation/UCSC_CpG_islands.gff`
  - `01_data/03_intermediate/008_annotation/methylsites_all.RData`
  - `01_data/03_intermediate/005_correlation_data/cor_all.RData`
  - `01_data/03_intermediate/005_correlation_data/all_mix_cor_CpG_common.RData`
  - `01_data/03_intermediate/008_annotation/genage_human.csv`
- **Creates**:
  - `annotation_human` (GRanges)
  - `anno_CpG_NC` (CpG islands GRanges with normalized seqnames)
  - `CpGs` (all CpGs GRanges with metadata)
  - `CpGs_selected` (selected CpGs GRanges)
  - Region tallies and gene lists for all vs selected CpGs, stratified by correlation direction and species
- **Writes**:
  - `01_data/03_intermediate/009_gene_analysis/gene_list.csv`
  - `01_data/03_intermediate/009_gene_analysis/gene_list_sel.csv`
  - `01_data/03_intermediate/009_gene_analysis/gene_list_count.csv`
  - `03_results/01_figures/07_regions_*.pdf`
  - `03_results/01_figures/07_genes_sel_cor.pdf`

## Important implementation notes / assumptions

- Many scripts include **hard-coded `setwd()` paths** pointing to user-specific folders and older repository locations. Reproducibility currently depends on locally mirroring those paths or editing scripts.
- “Upstream” parts (`001`, `002`, parts of `003`) refer to external raw data and alignment outputs that are not present in this repo snapshot.
- The downstream pipeline (from `004` onward) is consistent with the README statement that scripts `>= 003` rely on non-public inputs, while `004+` can reproduce paper results when required intermediate `.RData` are present.
