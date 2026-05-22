#!/bin/bash
# Align sequences to reference genome
# Requires to run 02_code/01_pipeline/01_scripts/01_CpGs_to_sequences.R first

# modules
echo "loading modules..."
# module load bowtie2
# module load samtools

# setting up variables 
echo "setting up variables..."
path_raw=01_data/01_raw_private/
path_rgenome=01_data/01_raw_private/rgenome/GRCh38
path_results=01_data/02_external_server_outputs/02_conserved_seq/
path_sequences=01_data/02_external_server_outputs/01_sequences/

# filenames
seq_1000bp=("AC_CpG_1000bp" "AS_CpG_1000bp" "EH_CpG_1000bp" "ZF_CpG_1000bp")

# file ending
suffix=".fasta"

# arguments
bowtie2_args="--very-sensitive --local -p 8"

# alignment
echo "aligning to rgenome GRCh38..."
echo "using bowtie2 with arguments: $bowtie2_args"

for (( i=0; i<${#seq_1000bp[@]}; i++ )); 
do 
    SAMPLE_NAME=${seq_1000bp[$i]}
    SAM_FILE="${path_results}HS_${SAMPLE_NAME}.sam"
    BAM_FILE="${path_results}HS_${SAMPLE_NAME}.bam"

    echo "Processing $SAMPLE_NAME..."

    # 1. Aligning
    bowtie2 $bowtie2_args -x ${path_rgenome} \
    -f -U ${path_sequences}${SAMPLE_NAME}$suffix -S "$SAM_FILE" -N 1

    # 2. Convert SAM to BAM, Sort it, and Save
    samtools view -bS "$SAM_FILE" | samtools sort -o "$BAM_FILE"

    # 3. Create a BAM Index
    samtools index "$BAM_FILE"

    # 4. Clean up the large SAM file
    rm "$SAM_FILE"

    echo "Done with $SAMPLE_NAME. BAM created and indexed."
done

echo "done!"