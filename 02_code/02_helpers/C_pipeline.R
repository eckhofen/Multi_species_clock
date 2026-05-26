# Metadata ----------------------------------------------------------------
# Project: Shared methylation regions
# Description: Orchestrator function for workflow "C"
# Author: Gabriel Ecker-Eckhofen
# eckhofen@icm.csic.es
# Date: 2026 05

# This script allows to reproduce the results from the paper. 
# Please check the README.md file for more information.

message("Running pipeline C...\n")

tryCatch({
    # 1/7 Data preparation
    source("02_code/01_pipeline/00a_folder_structure.R")
    source("02_code/01_pipeline/00b_data_preparation.R")
    source("02_code/01_pipeline/00c_data_overview.R")
    rm(list = ls())

    # 2/7 Data analysis
    source("02_code/01_pipeline/03_SMR_plotting.R")
    source("02_code/01_pipeline/04_correlation_testing.R")
    rm(list = ls())

    # 3/7 Model creation and testing
    source("02_code/01_pipeline/05a_model_creation.R")
    source("02_code/01_pipeline/05b_model_testing.R")
    rm(list = ls())

    # 4/7 SMR analysis
    source("02_code/01_pipeline/06_SMR_analysis.R")
    rm(list = ls())

    # 5/7 Annotation preparation
    source("02_code/01_pipeline/07_annotation_prep.R")
    rm(list = ls())

    # 6/7 Annotation
    source("02_code/01_pipeline/08_genomic_context_analysis.R")
    rm(list = ls())

    # Optional (files required)
    try(
        source("02_code/01_pipeline/09_auto-correlation.R"),
        rm(list = ls())
    )
    
    # 7/7 Additional data comparisons
    source("02_code/01_pipeline/10_data_comparison.R")
    rm(list = ls())

}, error = function(e) {
    message("Error in pipeline C: ", e$message)
    stop(e)
})

message("Pipeline C completed successfully!\n")