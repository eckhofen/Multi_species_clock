# Suggestions and improvements

## Reproducibility and portability (highest impact)

- **Remove hard-coded `setwd()` paths**
  - Several scripts hardcode absolute paths (e.g., `/Users/...` or `/powerplant/...`).
  - Replace with a project-root based approach (e.g., `here::here()`), or an explicit `config.R`.

- **Define a single source of truth for folder layout**
  - Standardize output folders: e.g. `figures/`, `002_plots/`, `000_data/...`.
  - Right now some scripts save to `figures/` while others save to `002_plots/`.

- **Add dependency capture**
  - Use `renv` to snapshot package versions.
  - Document system dependencies for Bioconductor packages.

- **Make scripts runnable from a clean session**
  - Many scripts assume objects exist from previous runs.
  - Add explicit `load()` statements and guard rails (e.g., stop with clear error if expected `.RData` is missing).

## Code quality and maintainability

- **Convert repeated blocks into shared functions**
  - Examples repeated across scripts:
    - color palette definitions
    - train/test splitting
    - plotting templates
    - the `evaluate.model()` function (duplicated between `005a` and `005b`)

- **Use consistent naming conventions**
  - Prefer snake_case and consistent suffixes, e.g. `*_meth_values`, `*_meth_values_imputed`, `*_methyl_sites`.
  - Avoid ambiguous names like `xx`.

- **Reduce duplication in analysis scripts**
  - `005a_statistical_models_REL.R` and `005b_statistical_models.R` are mostly identical.
  - Factor out:
    - the model training/eval code
    - only parameterize the target variable and transformations.

- **Make plotting deterministic**
  - Some uses of jitter are seeded (good), but many steps rely on default randomness.
  - Standardize `set.seed()` per script for plot reproducibility.

## Data handling improvements

- **Avoid using `rbind()` on tibbles with different schemas**
  - Prefer `dplyr::bind_rows()` and ensure consistent columns.

- **Explicitly store sample identifiers**
  - Some datasets likely depend on row order matching ages.
  - Add an explicit `sample_id` column and join by ID rather than relying on order.

- **Imputation strategy**
  - `zoo::na.aggregate()` is a simple aggregation imputer.
  - Consider reporting imputation rate per site/sample and evaluating sensitivity:
    - site filtering by missingness threshold
    - alternative imputation (kNN / missForest) depending on scale

## Statistical / modeling workflow

- **Use a unified modeling framework**
  - Consider `tidymodels` end-to-end (recipes + parsnip + workflows) instead of mixing caret/base/glmnet.

- **Cross-validation design**
  - The current split is stratified per species and then merged.
  - Consider nested CV or grouped CV by species if the goal is cross-species generalization.

- **Report uncertainty and calibration**
  - Add confidence intervals / prediction intervals where possible.
  - Add residual diagnostics and calibration plots.

## Outputs and documentation

- **Add a single “run order” document**
  - Current ordering is implicit in script numbering.
  - Provide a minimal “how to reproduce” section: which scripts can run with public data vs private.

- **Track generated outputs**
  - Add a lightweight manifest (e.g., `outputs_manifest.csv`) listing plots/data files and which script generates them.

## Quick wins (can be done immediately)

- **Create a project-level `config.R`**
  - Provide paths for `data_dir`, `out_dir`, `plots_dir`.

- **Replace hard-coded constants with a single parameter file**
  - Max lifespans, species mappings, naming maps.

- **Add input validation**
  - Check required columns exist in input dataframes before running major steps.
