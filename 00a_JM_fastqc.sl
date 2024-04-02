#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=2
#SBATCH --array=1-47%10
#SBATCH --output=/workspace/cfngle/raw-data/JM/fastqc_%A_%a.out

# Load the FastQC module
module load FastQC

# Define your project directory
PROJECT="/workspace/cfngle/raw-data/JM"

# Define input and output directories
IN="${PROJECT}/raw-reads"
OUT="${PROJECT}/001.fastqc_raw"

# Find all .fastq files and put them into an array
FILES=($(find ${IN} -type f -name "*.fastq"))

# Calculate index for this task
INDEX=$((SLURM_ARRAY_TASK_ID-1))

# Select file based on task array ID
FILE=${FILES[$INDEX]}

# Define command to run FastQC on the selected file
COMMAND="fastqc --nogroup -q -t ${SLURM_CPUS_PER_TASK} -o ${OUT} ${FILE}"

# Execute the command
echo "Executing: $COMMAND"
$COMMAND