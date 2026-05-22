# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: Create the complete folder structure used by the pipeline.
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es

# Each entry below is created if missing
folders <- c(
  # Pre-processing outputs (per-species subfolders are created by the
  # 02_code/00_pre-processing/<SPECIES>/ scripts themselves).
  "01_data/00_pre-processing",

  # Raw / private inputs the user must populate (see README).
  "01_data/01_raw_private",
  "01_data/01_raw_private/genomes",   # per-species reference FASTAs
  "01_data/01_raw_private/rgenome",   # human bowtie2 index used as alignment target

  # External / server-produced inputs copied into the repo.
  "01_data/02_external_server_outputs/01_sequences",
  "01_data/02_external_server_outputs/02_conserved_seq",

  # Pipeline-generated intermediate files.
  "01_data/03_intermediate/03_SMR",
  "01_data/03_intermediate/04_methyl_values",
  "01_data/03_intermediate/05_correlation_data",
  "01_data/03_intermediate/06_model_creation",
  "01_data/03_intermediate/07_data_comparison",
  "01_data/03_intermediate/08_annotation",
  "01_data/03_intermediate/09_co-methylation",

  # Metadata.
  "01_data/04_metadata",

  # Results.
  "03_results/01_figures",
  "03_results/02_tables"
)

invisible(lapply(folders, dir.create, recursive = TRUE, showWarnings = FALSE))

message("Created (or confirmed) ", length(folders), " folders:")
for (f in folders) message("  ", f)
