# CROP-seq
Repository for the processing and analysis of CROP-seq data. This repository contains code to do the following steps for CROP-seq analysis:

1. Reference preparation
2. Guide quantification from scRNA-seq data
3. Guide quantificationf from enriched RNAseq data

## Reference preparation
The reference requires the following information:

- Primary assembly of the organism ([FASTA](ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz) and [GTF](ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.chr.gtf.gz))
- [Spreadsheet with gRNA sequence and target genes](example_data/GWAS_sgRNAs.csv)
- [Spreadsheet with CROP-seq vector elements](‎⁨example_data/CROP-seq.csv⁩)

1. Use the following script: [buildReference.py](py/buildReference.py) to build the guide reference.

```bash
python buildReference.py -i example_data/GWAS_sgRNAs.csv -c example_data/CROP\-seq.csv -o Test
```

2. Concatenate the generated files to your primary assembly files, to create a new set of reference files.

```bash
> cat Homo_sapiens.GRCh38.dna.primary_assembly.fa CROPseq.fa > Homo_sapiens.CROPseq.fa
> cat Homo_sapiens.GRCh38.98.chr.gtf CROPseq.gtf > Homo_sapiens.CROPseq.gtf
```

3. Convert the files to Cell Ranger reference files using *cellranger mkref*.

```bash
# Path to input files
GENOME_NAME=Homo_sapiens.CROPseq
INPUT_FASTA=Homo_sapiens.CROPseq.fa
INPUT_GTF=Homo_sapiens.CROPseq.gtf
INPUT_FILTERED_GTF=Homo_sapiens.CROPseq.Filtered.gtf

# Filter GTF
cellranger mkgtf ${INPUT_GTF} ${INPUT_FILTERED_GTF}  --attribute=gene_biotype:protein_coding --attribute=gene_biotype:lincRNA --attribute=gene_biotype:antisense || exit 1

# Create STAR reference
cellranger mkref --genome=${GENOME_NAME} --fasta=${INPUT_FASTA} --genes=${INPUT_FILTERED_GTF} || exit 1
```

## FASTQ Preparation
1. If not done already, split FASTQ into paired-end reads using BBMAP suite.

```bash
reformat.sh in=$INPUT_FASTQ out1=$INPUT_FASTQ1 out2=$INPUT_FASTQ2
```

2. Trim and label reads with cell barcodes and UMI barcodes using umi_tools.
3. Trim polyA from R2 sequence with cutadapt.

Both steps can be done with [prepareFastqs.sh](shell/prepareFastqs.sh). Note that you may need to adjust the nucleotide sequences for different chemistries of the Chromium system.

## Guide Mapping
Map guides to the reference generated in the first section using STAR aligner.

```bash
STAR --runThreadN 8 \
--genomeDir $REF_DIR \
--readFilesIn $R2_EXTRACTED \
--readFilesCommand zcat \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Fastx \
--quantMode GeneCounts || exit 1;
```


