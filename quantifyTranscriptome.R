library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)

open_gzip <- function(x, type = c("r", "rb")){
  filestream <- gzfile(x, open = type)
  return(filestream)
}

input_dir <- "/Volumes/LACIE/CROP-seq/Transcriptome/AH_array-6/filtered_feature_bc_matrix/"
matrix <- readMM(open_gzip(paste0(input_dir, "matrix.mtx.gz"), type = "rb"))
barcodes <- read.csv(open_gzip(paste0(input_dir, "barcodes.tsv.gz"), type = "r"), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
genes <- read.csv(open_gzip(paste0(input_dir, "features.tsv.gz"), type = "r"), header = FALSE, stringsAsFactors = FALSE, sep = "\t")

colnames(barcodes) <- c("cell_barcode")
colnames(genes) <- c("ensembl_gene_id", "gene_name", "feature_type")
genes$gene_name <- make.unique(genes$gene_name)

colnames(matrix) <- barcodes$cell_barcode
rownames(matrix) <- genes$gene_name

# Grab list of gRNAs
grna_list <- grep("_gene", genes$gene_name, value = TRUE)
grna_matrix <- matrix[grna_list, ]

# Identify cells with at least one guide
postives <- as.data.frame(table(apply(grna_matrix, 2, function(x) any(x > 0))))
colnames(postives) <- c("gRNA", "nCells")
