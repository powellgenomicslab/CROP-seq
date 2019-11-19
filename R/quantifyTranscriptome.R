#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Seurat)
library(readr)

input_dir <- args[1]
sample_name <- args[2]

data <- Read10X(input_dir)
seurat_obj <- CreateSeuratObject(data, project = "CROP-seq", assay = "RNA")

# Grab list of gRNAs
grna_list <- grep("-gene", rownames(seurat_obj), value = TRUE)
counts <- seurat_obj[["RNA"]]@counts
grna_matrix <- counts[grna_list, ]

# For each cell, identify most abundant guide, discard if there are other molecules
transcriptome_assignments <- apply(grna_matrix, 2, function(x){
  if(any(x > 0)){
    max_val <- max(x)
    max_guide <- names(x)[which(x == max_val)]
    return(max_guide)
  } else{
    return(NA)
    }
})

transcriptome_df <- as.data.frame(unlist(transcriptome_assignments))
transcriptome_df$cell_barcode <- rownames(transcriptome_df)
rownames(transcriptome_df) <- NULL
colnames(transcriptome_df) <- c("gRNA", "cell_barcode")
transcriptome_df <- transcriptome_df[, c("cell_barcode", "gRNA")]
write_tsv(transcriptome_df, sprintf("%s_TranscriptAssignments.tsv", sample_name))