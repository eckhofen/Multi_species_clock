# changes

## Unreleased

- Updated `02_code/01_pipeline/01_scripts/03_methylation_extraction.R` to create required intermediate output folders before `save()` and to exit cleanly when private raw methylation inputs under `01_data/01_raw_private/` are not present.

- Updated `02_code/01_pipeline/01_scripts/04_correlation_testing.R` to check for required `01_data/03_intermediate/004_methyl_values/*.Rdata` inputs and exit cleanly with a message when they are missing.

- Updated `02_code/01_pipeline/01_scripts/05a_statistical_models_REL.R` to check for required `01_data/03_intermediate/006_model_creation/all_meth_values_selected.RData` input and exit cleanly with a message when it is missing.

- Updated `02_code/01_pipeline/01_scripts/05b_statistical_models.R` to check for required `01_data/03_intermediate/006_model_creation/all_meth_values_selected.RData` input and exit cleanly with a message when it is missing.

- Updated `02_code/01_pipeline/01_scripts/06_data_comparison.R` to check for required comparison/correlation intermediate inputs and exit cleanly with a message when they are missing.

- Updated `02_code/01_pipeline/01_scripts/07_gene_analysis.R` to check for required annotation and correlation intermediate inputs and exit cleanly with a message when they are missing.

- Add documentation: codebase overview, suggested improvements, and refactoring roadmap.
- Add `changes.md` and `code.md` for change tracking and object/function inventory.
- Update `02_code/01_pipeline/01_scripts/00_data_preparation.R` to exclude JM and use `01_data/04_metadata/` for palettes and metadata.
- Move server-only sequence metadata CSVs to `01_data/02_external_server_outputs/01_sequences/` and archive JM sequence metadata under `01_data/99_archive/JM/`.
- Refactor pipeline scripts `03_methylation_extraction.R`, `04_correlation_testing.R`, `05a_statistical_models_REL.R`, `05b_statistical_models.R`, and `06_data_comparison.R` to use the `01_data/` + `03_results/` folder structure, remove JM, and standardize SCMR naming (with backward-compatible handling of older `SMR` columns where needed).
- Run parse-only smoke tests (`Rscript -e "parse(file=...)"`) for the refactored scripts to confirm no syntax/path string issues.
- Update `README.md` to use SCMR terminology and describe the refactored folder structure and server-only vs locally reproducible pipeline split.
- Document the expected external server-only inputs in `README.md` (`01_data/02_external_server_outputs/*`) and which downstream scripts consume them.
- Initialize a git repository and update `.gitignore` to exclude result outputs (`03_results/`) and archive folders.
- Update `CODEBASE_OVERVIEW.md` (and `04_docs/` copy) to SCMR terminology and add explicit documentation of server-only external inputs vs locally reproducible downstream steps.
- Update `code.md` and `04_docs/code.md` to label external inputs vs locally reproducible intermediates and to point key plot outputs to `03_results/01_figures/`.
- Update `04_docs/README.md` to describe the tidy folder structure and external inputs in the same way as the top-level `README.md`.
- Refactor `02_code/01_pipeline/01_scripts/07_gene_analysis.R` to repo-relative paths, remove JM palette usage, and route outputs to `01_data/03_intermediate/009_gene_analysis/` and `03_results/01_figures/`.
- Run parse-only smoke test for `007_gene_analysis.R`.
- Archive deprecated legacy root folders (e.g. `001_scripts/`, `000_data/`, `002_plots/`, `code/`, `data/`) under `02_code/99_archive/deprecated_root/` (gitignored) and remove them from git tracking.
- Promote the active pipeline scripts folder to `02_code/01_pipeline/01_scripts/`.
- Enforce two-digit indexing for active pipeline script filenames (`03_...` instead of `003_...`) and move legacy/non-refactored scripts out of the active pipeline into a gitignored archive under `02_code/99_archive/`.
- Track `01_data/02_external_server_outputs/` and `01_data/04_metadata/` in git, while keeping derived intermediates (`01_data/03_intermediate/`), private raw inputs (`01_data/01_raw_private/`), and archives (`01_data/99_archive/`) excluded via `.gitignore`.
- Tweak point sizes in `02_code/01_pipeline/01_scripts/01_data_overview.R` and refresh `01_data/04_metadata/*.(RData)` artifacts.
- Restore server-only upstream scripts into a tracked folder `02_code/01_pipeline/00_server_only/01_scripts/` (two-digit names) for completeness of the pipeline documentation.
- Document the post-data-tracking tweaks: adjusted plotting point sizes in 01_data_overview.R and refreshed metadata/palette RData artifacts.
- Recover missing external overlap object `HS_AC_AS_EH_ZF_overlaps.Rdata` into `01_data/02_external_server_outputs/02_conserved_seq/` from the archived legacy `000_data/` snapshot to enable running `03_methylation_extraction.R`.
