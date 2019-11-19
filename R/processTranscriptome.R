library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)

sample_name <- "POOLED2"
input_dir <- "/Volumes/LACIE/CROP-seq/Transcriptome/AH_pool-2/filtered_feature_bc_matrix/"
assignment_df <- read_tsv(sprintf("~/Dropbox (Garvan)/Anne/CROP-seq/FinalCellAssignments/%s_FinalAssignments.tsv", sample_name))
data <- Read10X(input_dir)
seurat_obj<- CreateSeuratObject(data, project = "CROP-seq", assay = "RNA")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")

metadata <- seurat_obj@meta.data
metadata$cell_barcode <- rownames(metadata)
assignment_df <- as.data.frame(assignment_df)
rownames(assignment_df) <- assignment_df$cell_barcode
colnames(assignment_df) <- c("cell_barcode", "type", "gRNA", "target_gene")

metadata <- left_join(metadata, assignment_df, by = "cell_barcode")
#assignment_df <- assignment_df[, c("type", "gRNA",)]

# Add back to seurat obj
seurat_obj[["type"]] <- metadata$type
seurat_obj[["gRNA"]] <- metadata$gRNA
seurat_obj[["target_gene"]] <- metadata$target_gene

# QC plots
umi_hist <- ggplot(metadata, aes(nCount_RNA)) + geom_histogram(binwidth = 1000) + 
  geom_vline(xintercept = 9000, color = "red") + geom_vline(xintercept = 60000, color = "red") + theme_classic() + 
  ggtitle(sprintf("Total UMIs for %s", sample_name)) + xlab("Total UMIs") + ylab("Number of cells")
feature_hist <- ggplot(metadata, aes(nFeature_RNA)) + geom_histogram(binwidth = 100) + 
  geom_vline(xintercept = 3000, color = "red") + theme_classic() + 
  ggtitle(sprintf("Total Features for %s", sample_name)) + xlab("Total Features") + ylab("Number of cells")
mt_hist <- ggplot(metadata, aes(percent.mt)) + geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = 20, color = "red") + theme_classic() + 
  ggtitle(sprintf("Mitochondrial Expression for %s", sample_name)) + xlab("% Mt Expression") + ylab("Number of cells")
rb_hist <- ggplot(metadata, aes(percent.rb)) + geom_histogram(binwidth = 1) + 
  geom_vline(xintercept = 40, color = "red") + theme_classic() + 
  ggtitle(sprintf("Ribosomal Expression for %s", sample_name)) + xlab("% Rb Expression") + ylab("Number of cells")

pdf(sprintf("%s_QC.pdf", sample_name), width = 11, height = 8.5)
grid.arrange(umi_hist, feature_hist, mt_hist, rb_hist, ncol = 2)
dev.off()

# Filter data
seurat_obj <- subset(seurat_obj, nCount_RNA > 9000 & nCount_RNA < 60000 & nFeature_RNA > 3000 & 
                       percent.mt < 20 & percent.rb < 40)
saveRDS(seurat_obj, sprintf("/Volumes/LACIE/CROP-seq/%s_SeuratObject.rds", sample_name))


