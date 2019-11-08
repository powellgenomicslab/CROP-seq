#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)

open_gzip <- function(x, type = c("r", "rb")){
  filestream <- gzfile(x, open = type)
  return(filestream)
}

input_dir <- args[1]
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
positives <- as.data.frame(table(apply(grna_matrix, 2, function(x) any(x > 0))))
colnames(positives) <- c("gRNA", "nCells")

# Identify number of cells per guide
cells_per_guide <- as.data.frame(apply(grna_matrix, 1, function(x) length(which(x != 0))))
cells_per_guide$gRNA <- rownames(cells_per_guide)
rownames(cells_per_guide) <- NULL
cells_per_guide <- cells_per_guide[, 2:1]
colnames(cells_per_guide) <- c("gRNA", "Number of Cells")
write.csv(cells_per_guide, "Array2-CellsPerGuide.csv", row.names = FALSE)