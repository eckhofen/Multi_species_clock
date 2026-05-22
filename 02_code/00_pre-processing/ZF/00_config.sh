#!/bin/bash
# Shared configuration for the ZF RRBS preprocessing + Bismark alignment workflow.

# PROJECT is the root output directory (repo-relative; override via `export PROJECT=...`)
: "${PROJECT:=01_data/00_pre-processing/ZF}"

# Per-species pre-processing root inside the SMR repo

RAW_READS="${PROJECT}/raw-reads"
FASTQC_RAW="${PROJECT}/01_fastqc_raw"
TRIM_GALORE="${PROJECT}/02_trim_galore"
FASTQC_TRIMMED="${PROJECT}/03_fastqc_trimmed"
GENOME="01_data/01_raw_private/genomes"
ALIGNMENTS="${PROJECT}/04_alignments_bismark_local"
METH_EXTRACT="${PROJECT}/05_meth_extraction"
MULTIQC_OUT="${PROJECT}/06_multiqc"

# trim_galore v0.4.1
# bismark v0.23.0
# MultiQC v1.11
# Reference genome: GRCz11 (e.g. GCF_000002035.6_GRCz11_genomic.fasta)

export PROJECT RAW_READS FASTQC_RAW TRIM_GALORE FASTQC_TRIMMED \
       GENOME ALIGNMENTS METH_EXTRACT MULTIQC_OUT
