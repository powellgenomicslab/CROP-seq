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
def getSeq(x, cropseq_df):
    seq = cropseq_df.loc[cropseq_df["component"] == x, "sequence"].tolist()[0]
    return(seq)

def buildCas9(cropseq_df):
    seqs = [getSeq(element, cropseq_df = cropseq_df) for element in ["cas9", "nls", "flag", "p2a", "blast", "space", "virus_ltr"]]
    cas9_seq = ("").join(seqs)
    return(cas9_seq)

def buildRef(row, cropseq_df = None):
    sgRNA = row['sgRNA_Name']
    ref_seq = buildSeq(row["Sequence"], cropseq_df = cropseq_df)
    
    # WRITE FASTA
    header = fasta_header_template.format(chrom=sgRNA, length=len(ref_seq))
    fasta_lines.append(header)
    fasta_lines.append(ref_seq)

    # WRITE GTF
    gtf_lines.append(gtf_template.format(chrom=sgRNA, id=sgRNA, length=len(ref_seq)))
    
def buildSeq(x, cropseq_df = None):
    u6_seq = getSeq("u6", cropseq_df)
    rest_seq = getSeq("rest", cropseq_df)
    seq = u6_seq + x + rest_seq
    return(seq)

# Parsers
def loadFiles(args):
    cropseq_df = pd.read_csv(args.cropseq)
    guide_df = pd.read_csv(args.input)
    return(cropseq_df, guide_df)

# Argument parser
def parseArgs():
    parser = argparse.ArgumentParser(prog = "PrepRef")
    parser.add_argument("-i", "--input", type = str, help = "Input gRNA file")
    parser.add_argument("-c", "--cropseq", type = str, help = "CROP-seq file")
    parser.add_argument("-o", "--output", type = str, help = "Basename of output files")
    args = parser.parse_args()
    return(args)

if __name__ ==  "__main__":
    # Parse arguments
    args = parseArgs()

    # Parse spreadsheets
    cropseq_df, guide_df = loadFiles(args)

    guide_df.apply(buildRef, cropseq_df = cropseq_df, axis = 1)

    seq_name = "Cas9_blast"
    sequence = buildCas9(cropseq_df)
    header = fasta_header_template.format(chrom=seq_name, length=len(sequence))

    # add cas9 entry
    fasta_lines.append(header)
    fasta_lines.append(sequence)
    gtf_lines.append(gtf_template.format(chrom=seq_name, id=seq_name, length=len(sequence)))

    # write to file
    output_fasta = args.output + ".fa"
    output_gtf = args.output + ".gtf"
    
    with open(output_fasta, "w") as fasta_handle:
        fasta_handle.writelines("\n".join(fasta_lines))
    with open(output_gtf, "w") as gtf_handle:
        gtf_handle.writelines(gtf_lines)
