---
title: "Shared conserved methylation regions and multispecies piscine epigenetic clock"
author: "Gabriel Ecker-Eckhofen"
date: "June 2024"
---


### Overview
This repository has been made for a masters thesis project: **Shared Conserved Methylation Regions (SCMRs): a New Method in Comparative Epigenetics Used to Create a Multispecies Piscine Epigenetic Clock**, which has not been made publicly available yet. 


### Introduction
In this project we introduced a new workflow for identifying comparable methylation data across species. This workflow involves 1) expanding CpG methylation sites from methylation data to sequences (in our case 1kb long) using the species genome, 2) aligning these sequences onto a common reference genome for all species and 3) identifying overlapping aligned sequences which we named **shared conserved methylation regions (SCMRs)**.

These SCMRs now included CpG sites from all species which 4) we then used for age correlation testing. *Note, that this could be done with any variable you can test correlation for*. 5) SCMRs were grouped into two classes. One class contained CpGs which were positively correlating with age and the other was containing negatively correlating ones. 6) Now, we selected only the highest correlating CpG for each species and CpG in every SCMR group. This leaves one CpG (for each species) in each SCMR. 7) If there was no agreement in correlation direction within an SCMR (not all species had CpGs correlating in the same direction), the SCMR was excluded. If there were SCMRs with positively correlating and negatively correlating CpGs, only the negatively ones were retained. 

This finally allowed us to 8) use the retained SCMRs as independent variables for creating a methylation-based age prediction model (also known as an "epigenetic clock"). 9) We tested various models such as multivariate linear regression and non-parametric random forest regression. We achieved notable accuracy in age prediction for four species using a single model ("multispecies epigegentic clock"), the results of which will be published in a paper separate from the thesis. 


### Using this repo
As of now, crucial data which has been used in this project is not published yet. It is to be expected that, during the course of this year, all data sets will be made publicly available in separate articles.

This repository is currently being refactored into a tidy structure:

- `01_data/`: raw/private inputs, external server outputs, intermediate (derived) files, metadata
- `02_code/`: analysis code (pipeline scripts, helpers)
- `03_results/`: figures/tables/logs
- `04_docs/`: documentation

The pipeline includes a clear split between:

- **Server-only / external preprocessing inputs**: alignment outputs and sequence-level metadata produced outside this repo and copied into `01_data/02_external_server_outputs/`.
- **Locally reproducible downstream analysis**: scripts that start from those external inputs and produce intermediate datasets in `01_data/03_intermediate/` and figures in `03_results/`.

Server-only / external inputs expected by the downstream pipeline:

- `01_data/02_external_server_outputs/01_sequences/*_metadata_*1000bp.csv`
  - Consumed by `02_code/01_pipeline/01_scripts/002_rgHS_conserved_sequences.R`.
- `01_data/02_external_server_outputs/02_conserved_seq/HS_AC_AS_EH_ZF_overlaps.Rdata`
  - Consumed by `02_code/01_pipeline/01_scripts/003_methylation_extraction.R`.

To reproduce results which have been shown in the master thesis, scripts can be looked at and (depending on data availability) run.


### Documentation
- `CODEBASE_OVERVIEW.md`: Script-by-script overview (data inputs/outputs, objects created, plots generated).
- `code.md`: Inventory of key functions, objects, and stored datasets.
- `SUGGESTIONS_AND_IMPROVEMENTS.md`: Suggested improvements and best-practice recommendations.
- `REFACTORING_ROADMAP.md`: A phased roadmap for refactoring and maintainability.
- `changes.md`: Changelog for repository changes.


#### Disclaimer 
- Sequencing data is not published yet and therefore, only scripts involved in the downstream processes can be used as of now.
- This repository is under development and is therefore subject to changes.



For any questions, recommendations or remarks, reach out to me via eckhofen@icm.csic.es
