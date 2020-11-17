#!/usr/bin/env python
import argparse
import pysam
import pandas as pd
import numpy as np

def quantifyGuides(x, bam_input = None):
    columns = ["gRNA", "cell_barcode", "umi", "reads"]
    molecule_df = pd.DataFrame(columns = columns)
    row_df = pd.DataFrame({"gRNA":np.chararray(10000),"cell_barcode":np.chararray(10000), \
        "umi":np.chararray(10000),"reads":np.zeros(10000)})
    i = 0

    for read in bam_input.fetch(x + "_chrom"):
        # Parse read to get barcodes
        tag = (read.query_name).split(":")[-1]
        cell_barcode = tag[0:16]
        umi = tag[16:-1]
        row_df.iloc[i % 10000,:] = [x, cell_barcode, umi, 1]
        i += 1

        if i % 10000 == 0:
            molecule_df = molecule_df.append(row_df, ignore_index = True)
            row_df = pd.DataFrame({"gRNA":np.chararray(10000),"cell_barcode":np.chararray(10000), \
                "umi":np.chararray(10000),"reads":np.zeros(10000)})

    return(molecule_df.iloc[:i-1,:])

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
