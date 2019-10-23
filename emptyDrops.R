library(Seurat)
library(ggplot2)

input_dir <- "/Volumes/LACIE/CROP-seq/CROP-seq/AH_pool-1/filtered_feature_bc_matrix/"
data <- Read10X(input_dir)
seurat_obj <- CreateSeuratObject(data, project = "CROP-seq")
metadata <- seurat_obj@meta.data

umi_content <- ggplot(metadata, aes(nCount_RNA)) + geom_histogram(binwidth = 1000) + geom_vline(xintercept = 77000, color = "red") + theme_classic()
feature_content <- ggplot(metadata, aes(nFeature_RNA)) + geom_histogram(binwidth = 100) + geom_vline(xintercept = 3000, color = "red") + geom_vline(xintercept = 7500, color = "red")

seurat_obj <- subset(seurat_obj, nCount_RNA < 77000 & nFeature_RNA > 3000 & nFeature_RNA < 7500)
cells <- Cells(seurat_obj)