#!/bin/bash
for FILE in *.fastqsanger; do
	# GET SAMPLE ID FROM FILENAME
	SAMPLE_ID=${FILE%.*}
	SAMPLE_ID=${SAMPLE_ID%%_*}
	SAMPLE_ID=$( cut -d "[" -f 2 <<<  "$SAMPLE_ID" )
	
	PROCESSING_FILE=${SAMPLE_ID}_Processed.fastq
	OUTPUT_FILE=${SAMPLE_ID}.fastq

	# Run fastp to label reads with cell barcodes and UMIs
	fastp -i ${FILE}  -o ${PROCESSING_FILE} -U --umi_loc read1 --umi_len 28  -w 4
	
	# Run cutadapt to strip polyT from beginning of read
	cutadapt -j 4 -g "T{200}" -o $OUTPUT_FILE $PROCESSING_FILE
done
