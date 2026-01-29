# changes

## Unreleased

- Add documentation: codebase overview, suggested improvements, and refactoring roadmap.
- Add `changes.md` and `code.md` for change tracking and object/function inventory.
- Update `02_code/01_pipeline/01_scripts/00_data_preparation.R` to exclude JM and use `01_data/04_metadata/` for palettes and metadata.
- Move server-only sequence metadata CSVs to `01_data/02_external_server_outputs/01_sequences/` and archive JM sequence metadata under `01_data/99_archive/JM/`.
- Refactor pipeline scripts `003_methylation_extraction.R`, `004_correlation_testing.R`, `005a_statistical_models_REL.R`, `005b_statistical_models.R`, and `006_data_comparison.R` to use the `01_data/` + `03_results/` folder structure, remove JM, and standardize SCMR naming (with backward-compatible handling of older `SMR` columns where needed).
- Run parse-only smoke tests (`Rscript -e "parse(file=...)"`) for the refactored scripts to confirm no syntax/path string issues.
- Update `README.md` to use SCMR terminology and describe the refactored folder structure and server-only vs locally reproducible pipeline split.
- Document the expected external server-only inputs in `README.md` (`01_data/02_external_server_outputs/*`) and which downstream scripts consume them.
- Initialize a git repository and update `.gitignore` to exclude result outputs (`03_results/`) and archive folders.
- Update `CODEBASE_OVERVIEW.md` (and `04_docs/` copy) to SCMR terminology and add explicit documentation of server-only external inputs vs locally reproducible downstream steps.
- Update `code.md` and `04_docs/code.md` to label external inputs vs locally reproducible intermediates and to point key plot outputs to `03_results/01_figures/`.
- Update `04_docs/README.md` to describe the tidy folder structure and external inputs in the same way as the top-level `README.md`.
- Refactor `02_code/01_pipeline/01_scripts/007_gene_analysis.R` to repo-relative paths, remove JM palette usage, and route outputs to `01_data/03_intermediate/009_gene_analysis/` and `03_results/01_figures/` (with legacy fallbacks for required annotation inputs where needed).
- Run parse-only smoke test for `007_gene_analysis.R`.
- Archive deprecated legacy root folders (e.g. `001_scripts/`, `000_data/`, `002_plots/`, `code/`, `data/`) under `02_code/99_archive/deprecated_root/` (gitignored) and remove them from git tracking.
- Promote the active pipeline scripts folder to `02_code/01_pipeline/01_scripts/`.
