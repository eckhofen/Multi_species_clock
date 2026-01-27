# Refactoring roadmap

This roadmap focuses on improving maintainability, reproducibility, and clarity while preserving current scientific behavior.

## Phase 0: Baseline and safety (1–2 sessions)

- **Freeze current behavior**
  - Identify which scripts are "source-of-truth" for current paper figures.
  - Run them once (locally) and record the expected output files.

- **Introduce `renv` (optional but recommended)**
  - Capture package versions used for the current results.

- **Add minimal runtime checks**
  - At the start of each script, assert that required input files exist.

## Phase 1: Path/config cleanup (high priority)

- **Remove all hard-coded `setwd()` calls**
  - Introduce `here::here()` (or equivalent) and require running from repo root.
  - Add a single `R/config.R` (or similar) that defines:
    - `paths$metadata_dir`
    - `paths$data_dir`
    - `paths$plots_dir`
    - `paths$figures_dir`

- **Standardize output locations**
  - Decide: either always write to `002_plots/` or `figures/`.
  - Update scripts to match.

## Phase 2: Extract shared utilities (medium priority)

- **Create a utilities module** (e.g. `R/utils_plotting.R`, `R/utils_modeling.R`, `R/utils_metadata.R`)
  - Move duplicated pieces:
    - palettes + species naming maps
    - `theme_custom()` and other plot scaffolds
    - `evaluate.model()`
    - correlation helpers (`cor.test.age*`, `select.max.cor`)

- **Add consistent naming and object schemas**
  - Standardize required columns across datasets (age, rel_age, species, sample_id).

## Phase 3: Consolidate modeling scripts (high value)

- **Unify `005a` and `005b`**
  - Create a single modeling driver that takes:
    - target (`rel_age` vs chronological)
    - transformation (none / log / -log-log)
    - y-axis/x-axis limits
  - Keep old scripts as thin wrappers for backward compatibility.

- **Move to `tidymodels` workflows**
  - Use `recipes` for preprocessing.
  - Use `parsnip` for models; avoid mixing caret/base.

## Phase 4: Pipeline orchestration (optional, but powerful)

- **Turn the numbered scripts into a pipeline**
  - Use `targets` or `drake` so outputs are cached and dependency-tracked.
  - Benefits:
    - re-run only what changed
    - automatic parallelism
    - clear DAG of data dependencies

- **Add an output manifest**
  - Each target registers produced plots and datasets.

## Phase 5: Testing and validation (recommended)

- **Add lightweight tests**
  - Use `testthat` for:
    - correlation selection rules
    - SMR filtering rules (direction agreement)
    - basic invariants (dimensions, missingness thresholds)

- **Add a smoke-test dataset**
  - Provide a tiny synthetic dataset to verify pipeline runs end-to-end without private data.

## Phase 6: Packaging / API surface (long-term)

- **Convert core logic into an R package**
  - Expose:
    - `build_metadata()`
    - `run_correlation_screen()`
    - `select_smr_features()`
    - `fit_clock_models()`
    - `annotate_cpgs()`
  - Keep analyses/figures as vignettes or scripts.
