#!/usr/bin/env python
# Python 3.7 
import argparse
import pandas as pd

# Global variables
fasta_lines = []
gtf_lines = []

# Taken from crop-seq
fasta_header_template = ">{chrom}_chrom dna:chromosome chromosome:GRCh38:{chrom}_chrom:1:{length}:1 REF"

gtf_template = """{chrom}_chrom\thavana\tgene\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA";
{chrom}_chrom\thavana\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana";
{chrom}_chrom\thavana\texon\t1\t{length}\t.\t+\t.\tgene_id "{id}_gene"; transcript_id "{id}_transcript"; exon_number "1"; gene_name "{id}_gene"; gene_source "ensembl_havana"; gene_biotype "lincRNA"; transcript_name "{id}_transcript"; transcript_source "havana"; exon_id "{id}_exon";
"""

# Build sequence
def buildRef(row):
    # Get row with these columns from the CSV file of gRNAs
    sgRNA = row["sgRNA"]
    ref_seq = row["sequence"]
    
    # WRITE FASTA
    header = fasta_header_template.format(chrom=sgRNA, length=len(ref_seq))
    fasta_lines.append(header)
    fasta_lines.append(ref_seq)

    # WRITE GTF
    gtf_lines.append(gtf_template.format(chrom=sgRNA, id=sgRNA, length=len(ref_seq)))
    
# Parsers
def loadFiles(args):
    guide_df = pd.read_csv(args.input)
    return(guide_df)

# Argument parser
def parseArgs():
    parser = argparse.ArgumentParser(prog = "PrepRef")
    parser.add_argument("-i", "--input", type = str, help = "Input gRNA file")
    parser.add_argument("-o", "--output", type = str, help = "Basename of output files")
    args = parser.parse_args()
    return(args)

if __name__ ==  "__main__":
    # Parse arguments
    args = parseArgs()

    # Parse spreadsheets
    guide_df = loadFiles(args)

    guide_df.apply(buildRef, axis = 1)

    # Write to file
    output_fasta = args.output + ".fa"
    output_gtf = args.output + ".gtf"
    
    with open(output_fasta, "w") as fasta_handle:
        fasta_handle.writelines("\n".join(fasta_lines))
    with open(output_gtf, "w") as gtf_handle:
        gtf_handle.writelines(gtf_lines)
