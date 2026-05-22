#!/bin/bash

# modules
echo "loading modules..."
module load bowtie2


# setting up variables
path_raw=01_data/01_raw_private
path_rgenome=01_data/01_raw_private/rgenome
path_results=01_data/02_external_server_outputs/02_conserved_seq
path_sequences=01_data/02_external_server_outputs/01_sequences

mkdir -p ${path_rgenome}

# indexxing reference genome
echo "building bowtie2 index..."
bowtie2-build ${path_raw}/GRCh38.fa ${path_rgenome}/GRCh38

echo "done!" 