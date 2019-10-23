#!/bin/bash

# Arguments - Input FASTQ files
# Paired end from MiSeq
R1_FASTQ=$1
R2_FASTQ=$2
OUTPUT_DIR=$3


# Define outputs
mkdir -p $OUTPUT_DIR
R1_TRIMMED=${OUTPUT_DIR}/$( basename ${R1_FASTQ} )
R2_TRIMMED=${OUTPUT_DIR}/$( basename {R2_FASTQ} )

# Activate conda environment
# Conda environment has umi_tools and cutadapt from bioconda
conda activate umi_tools

# BC pattern - C is cell barcode sequence, N is UMI sequence
# Get cell number from QC/Cell Ranger
umi_tools whitelist --stdin $R1_FASTQ \
--bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
--set-cell-number=9900 \
--log2stderr > whitelist.txt || exit 1;

# Extract barcodes and UMIs and add to read names
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNNNN \
--stdin $R1_FASTQ \
--stdout $R1_EXTRACTED \
--read2-in $R2_FASTQ \
--read2-out=$R2_EXTRACTED \
--filter-cell-barcode \
--whitelist=whitelist.txt || exit 1;

# Trim polyA/polyT tail from reads
# We know polyA tail is 30 nucleotides in length from the documentation
cutadapt -g "^T{30}" -e 0.25 -m 20 -o $R1_TRIMMED -p $R2_TRIMMED $R1_FASTQ $R2_FASTQ



