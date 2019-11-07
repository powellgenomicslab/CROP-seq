library(Seurat)
library(tidyr)
library(dplyr)
library(ggplot2)

# Filtered in previous step
data <- Read10X("/Volumes/LACIE/CROP-seq/Transcriptome/AH_pool-1/filtered_feature_bc_matrix/")
seurat_obj <- readRDS("/Volumes/LACIE/CROP-seq/AH_pool-1_FilteredObject.rds")
metadata <- seurat_obj@meta.data

# Get transcript assignments
control_guides <- grep("NonTargeting", rownames(seurat_obj), value = TRUE)
guides <- grep("GUIDES", rownames(seurat_obj), value = TRUE)

# Add guide assignments
seurat_obj[["sgRNA"]] <- rep(NULL, ncol(seurat_obj))

# Normalisation via variance stabilisation
# Regress out mitochondrial and ribosomal expression
#seurat_obj <- SCTransform(seurat_obj, vars.to.regress = c("percent.mt", "percent.rb"), 
#                          return.only.var.genes = FALSE)

#saveRDS(seurat_obj, "/Volumes/LACIE/CROP-seq/AH_pool-1_NormalisedObject.rds")

# Count matrix
counts <- seurat_obj[["RNA"]]@counts
metadata[, "gRNA"] <- NULL
metadata$sgRNA <- "None"
metadata$sgRNA[which(rownames(metadata) %in% rownames(cell_assignments))] <- cell_assignments$gRNA
