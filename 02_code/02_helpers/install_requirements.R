# Helper script to install all R software requirements

# CRAN
install.packages(c(
  "tidyverse", "tidymodels", "data.table", "reshape2",
  "patchwork", "ggrepel", "ggtext", "ggnewscale", "ggpattern",
  "ggforce", "ggfortify", "RColorBrewer",
  "caret", "Hmisc", "vegan", "zoo"
))

# Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "Biostrings", "GenomicRanges", "GenomicAlignments", "GenomicFeatures",
  "GenomeInfoDb", "AnnotationHub", "annotatr",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ggbio"
))