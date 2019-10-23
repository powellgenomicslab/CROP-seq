#!/bin/bash

# Arguments - Input FASTQ files
# Paired end from MiSeq
INPUT_DIR=/share/ScratchGeneral/annsen/data/experimental_data/RAW/CROP-seq
OUTPUT_DIR=${INPUT_DIR}/processed
R1_FASTQ=${INPUT_DIR}/fastq/POOLED1_S7_L001_R1_001.fastq
R2_FASTQ=${INPUT_DIR}/fastq/POOLED1_S7_L001_R2_001.fastq

# Define outputs
mkdir -p $OUTPUT_DIR
R1_FILENAME=$( basename "$R1_FASTQ" )
R2_FILENAME=$( basename "$R2_FASTQ" )

R1_TAGGED=${OUTPUT_DIR}/${R1_FILENAME}
R2_TAGGED=${OUTPUT_DIR}/${R2_FILENAME}

R1_TRIMMED=${OUTPUT_DIR}/${R1_FILENAME}
R2_TRIMMED=${OUTPUT_DIR}/${R2_FILENAME}

# Activate conda environment
# Conda environment has umi_tools and cutadapt from bioconda
conda activate umi_tools

# BC pattern - C is cell barcode sequence, N is UMI sequence
# Get cell number from QC/Cell Ranger
umi_tools whitelist --stdin $R1_FASTQ --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --set-cell-number=10000 --log2stderr > whitelist.txt

# Extract barcodes and UMIs and add to read names
umi_tools extract -I $R1_FASTQ --extract-method=string --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN --read2-in=$R2_FASTQ --stdout=$R1_TAGGED --read2-out=$R2_TAGGED --filter-cell-barcode --whitelist whitelist.txt

# Trim polyA/polyT tail from reads
# We know polyA tail is 30 nucleotides in length from the documentation
#cutadapt -g "^T{30}" -e 0.25 -m 20 -o $R1_TRIMMED -p $R2_TRIMMED $R1_TAGGED $R2_TAGGED
