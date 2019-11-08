#!/usr/bin/env python
import argparse
import pysam
import pandas as pd

def quantifyGuides(x, bam_input = None):
    columns = ["gRNA", "cell_barcode", "umi", "reads"]
    molecule_df = pd.DataFrame(columns = columns)

    for read in bam_input.fetch(x + "_chrom"):
        # Parse read to get barcodes
        tag = (read.query_name).split(":")[-1]
        cell_barcode = tag[0:16]
        umi = tag[16:-1]
        
        row_df = pd.DataFrame([(x, cell_barcode, umi, 1)], columns = columns) 
        molecule_df = molecule_df.append(row_df, ignore_index = True)

    return(molecule_df)        

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type = str, help = "Enrichment alignment file.")
    parser.add_argument("-o", "--output", type = str, help = "Name of the output.")
    args = parser.parse_args()
    return((args.input, args.output))

if __name__ == "__main__":
    # Parse arguments
    input_dir,output_dir = parseArgs()

    # Read in guide RNA
    guide_filepath = "../example_data/CropSeq_sgRNA_2019.csv"
    guide_df = pd.read_csv(guide_filepath)

    # Read in BAM file
    bam_filepath = input_dir

    # Gather cell barcodes
    bam_input = pysam.AlignmentFile(bam_filepath, "rb")
    guide_df_list = [quantifyGuides(x, bam_input = bam_input) for x in guide_df["sgRNA"].tolist() ] 
    bam_input.close()

    molecule_df = pd.concat(guide_df_list)
    molecule_tally = (molecule_df.groupby(["gRNA", "cell_barcode", "umi"])["reads"].sum()).reset_index()
    molecule_tally.to_csv(output_dir, index = False, sep = "\t")
