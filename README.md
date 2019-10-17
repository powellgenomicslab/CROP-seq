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

Use the following script: [buildReference.py](py/buildReference.py) to build the guide reference.

```bash
python buildReference.py -i example_data/GWAS_sgRNAs.csv -c example_data/CROP\-seq.csv -o Test
```


