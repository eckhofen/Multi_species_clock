---
title: "A framework for comparative epigenomics using heterogeneous datasets validated in multispecies clocks"
author: "Gabriel Ecker-Eckhofen"
date: "May 2026"
---

# Overview
This repository comes along with the project: **A framework for comparative epigenomics using heterogeneous datasets validated in multispecies clocks**. 

# Introduction
In this project we introduced a new workflow for identifying comparable methylation data across species. This workflow involves 1) expanding CpG methylation sites from methylation data to sequences (in our case 1kb long) using the species genome, 2) aligning these sequences onto a common reference genome for all species and 3) identifying overlapping aligned sequences which we named **shared methylation regions (SMRs)**.

These SMRs now included CpG sites from all species which 4) we then used for age correlation testing. *Note, that this could be done with any variable you can test correlation for*. 5) SMRs were grouped into two classes. One class contained CpGs which were positively correlating with age and the other was containing negatively correlating ones. 6) Now, we selected only the highest correlating CpG for each species and CpG in every SMR group. This leaves one CpG (for each species) in each SMR. 7) If there was no agreement in correlation direction within an SMR (not all species had CpGs correlating in the same direction), the SMR was excluded. If there were SMRs with positively correlating and negatively correlating CpGs, only the negatively ones were retained. 

This finally allowed us to 8) use the retained SMRs as independent variables for creating a methylation-based age prediction model (also known as an "epigenetic clock"). 9) We tested various models such as multivariate linear regression and non-parametric random forest regression. We achieved notable accuracy in age prediction for four species using a single model ("multispecies epigegentic clock"), the results of which will be published in a paper separate from the thesis. 

# Installation guide

Close the github repo:
```bash
git clone https://github.com/eckhofen/Multi_species_clock.git
```

To run this pipeline make sure that the [requirements](#requirements) are met, based on either workflow (A), (B), or (C). See below for more information.

## Instructions for use

Based on the available data, follow the appropriate workflow:
- **(A)** If raw methylation reads are available, then follow the `02_code/00_pre-processing/` workflow. 
- **(B)** If methylation positions for each species are available, then proceed to `02_code/01_pipeline.R` and follow the scripts via indices starting with `00a_data_preparation.R`. 
- **(C)** If no data is available, then follow the `02_code/01_pipeline/03_SMR_plotting.R` and continue chronologically from there, or use the convenient orchestrator script: `02_code/02_helpers/C_pipeline.R`. This allows to reproduce all results from the paper apart from the autocorrelaiton of the whole methylation data (`02_code/01_pipeline/09_auto-correlation.R`). Please request these files if reproduction of those results is required.

## Demo / Replicating results

See `replicate_results.md` for detailed instructions on reproducing all numbers and figures from the paper.

For reviewers, please read the `for_reviewers.md` file for additional information and access to required data.

# Requirements

### Hardware requirements
This pipeline has been tested on **macOS 15.5** and **Red Hat Enterprise Linux 9.7**. Hardware requirements are minimal (8 GB RAM recommended), apart from the pre-processing and alignment which largely benefit from having sufficient RAM and storage space. No special hardware needed. 

### Software requirements

#### Command-line tools (pre-processing / alignment) (A, B)
Used by the scripts in `02_code/00_pre-processing/`:

- [trim_galore](https://github.com/FelixKrueger/TrimGalore) v0.4.1
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (any recent version)
- [MultiQC](https://multiqc.info/) v1.11
- [Bismark](https://felixkrueger.github.io/Bismark/) v0.23.0
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) (used by Bismark for alignment / indexing)
- `samtools` (used by Bismark internally)
- A SLURM-style scheduler (`sbatch` + the local `abatch` wrapper) on the cluster — scripts can also be adapted to run locally.

#### R (downstream analysis) (A, B, C)
Tested with `R 4.0`, `R 4.3.3` and `R 4.5.2`

CRAN packages:

In R, use the helper script to install all required packages:
```r
source("02_code/02_helpers/install_requirements.R")
```

Or install manually:
```r
install.packages(c(
  "tidyverse", "tidymodels", "data.table", "reshape2",
  "patchwork", "ggrepel", "ggtext", "ggnewscale", "ggpattern",
  "ggforce", "ggfortify", "RColorBrewer",
  "caret", "Hmisc", "vegan", "zoo"
))
```

Bioconductor packages (install via `BiocManager::install(...)`):
```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "Biostrings", "GenomicRanges", "GenomicAlignments", "GenomicFeatures",
  "GenomeInfoDb", "AnnotationHub", "annotatr",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ggbio"
))
```

The hg38 genomic-context annotations used by `07_annotation_prep.R`
(`hg38_cpgs`, `hg38_basicgenes`, `hg38_lncrna_gencode`) are pulled from
AnnotationHub by `annotatr::build_annotations()` the first time the script
is run and are cached locally afterwards, meaning **no separate manual download
is required** beyond the reference genomes and the GenAge file listed below. 
The installation time may vary depending on the internet connection, but typically takes about 10-20 mins.

### Methylation pre-processing (required for (A))
The RRBS / BisRADseq data has to be preprocessed according to the methodology described in the respective publication. If raw reads are not available, then proceed with `02_code/01_pipeline.R` scripts. The zebrafish RRBS preprocessing + Bismark alignment workflow lives in `02_code/00_pre-processing/ZF/` (scripts `01_…08_`); outputs are written to `01_data/00_pre-processing/ZF/`.

The methylation data has to be requested for Australasian snapper and European hake, while the others can be downloaded here: 
- Atlantic cod: [methylation data](https://doi.org/10.5061/dryad.tmpg4f565) and related [paper](https://doi.org/10.1111/1755-0998.70109)
- Zebrafish: [methylation data](https://data.csiro.au/collection/csiro:46344v3) and related [paper](https://doi.org/10.18632/aging.202400)

### Reference genomes (required for (A, B))
The reference genomes used in this project are available from the respective sources:
- European hake: [fMerMel2.1_cnag1.scaffolds.fa](https://denovo.cnag.cat/filebrowser/download/1819)
- Atlantic cod: [GCF_902167405.1_gadMor3.0_genomic.fasta](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902167405.1/)
- Australasian snapper: Chrysophrys_auratus.v.1.0.all.male.map.fasta. Note, that this genome is managed by Aotearoa Genomic Data Repository (AGDR; [64](https://doi.org/10.1016/j.ygeno.2024.110929)) and can be accessed via application at https://data.agdr.org.nz/.
- Zebrafish: [GCF_000002035.6_GRCz11_genomic.fasta](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002035.6/)
- Human: [GRCh38.fa](https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/)

Place the per-species FASTAs into `01_data/01_raw_private/genomes/`.

The human reference (`GRCh38.fa`) is the alignment *target* used by `01b_index_rgenome.sh` and `01c_align_to_rgenome.sh`; place the FASTA at `01_data/01_raw_private/GRCh38.fa` and the resulting bowtie2 index will be written to `01_data/01_raw_private/rgenome/` automatically.

### Aging gene reference (required for (A, B, C))
- GenAge human: [genage_human.csv](https://genomics.senescence.info/genes/human_genes.zip). Place the unzipped csv into the `01_data/03_intermediate/08_annotation` folder.

# License
This project is licensed under the **GNU General Public License v3.0.**

## Contact

For any questions or inquiries, contact me at `eckhofen@icm.csic.es`.