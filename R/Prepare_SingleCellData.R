library(Seurat)
# Load aggregated array data from Cell Ranger
data <- Read10X("/Volumes/LACIE/CROP-seq/Transcriptome/ARRAYED_AGGR/outs/filtered_feature_bc_matrix/")

# Create Seurat object
seurat_obj <- CreateSeuratObject(data, project = "CROP-seq", assay = "RNA")

# Calculate percentage Mt and Rb expression
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")

# Read in list of post-filter guide assignments
assignment_file_list <- lapply(1:6, function(x) 
  read.csv(sprintf("/Users/anne/Dropbox (Garvan)/Anne/CROP-seq/FinalCellAssignments/ARRAY%d_PostFilterAssignments.tsv", x), 
           sep = "\t", stringsAsFactors = FALSE))

# Combine guide assignments into one data frame
names(assignment_file_list) <- as.character(1:6)
assignment_df <- bind_rows(assignment_file_list, .id = "ARRAY")

# Attatch array number as a suffix 
assignment_df$cell_barcode <- paste0(assignment_df$cell_barcode, "-", assignment_df$ARRAY)
assignment_df <- assignment_df %>% drop_na()

# Get list of cells to keep - this replaces the filtering step as we have already
# done this on individual samples
kept_cells <- intersect(Cells(seurat_obj), assignment_df$cell_barcode)

# Filter cells based on list of cells we want to keep
seurat_obj <- subset(seurat_obj, cells = kept_cells)

# Rearrange data frame so we can add this information to the Seurat object
assignment_df <- assignment_df[match(Cells(seurat_obj), assignment_df$cell_barcode), ]
seurat_obj[["ARRAY"]] <- assignment_df$ARRAY
seurat_obj[["TYPE"]] <- assignment_df$type
seurat_obj[["gRNA"]] <- assignment_df$gRNA
seurat_obj[["TARGET"]] <- assignment_df$target_gene

# Save to external HD
saveRDS(seurat_obj, "/Volumes/LACIE/CROP-seq/ARRAY_AGGR_OBJ.rds")

# Normalise using SCtransform, with regression for percent.mt and percent.rb,
# in addition to batch normalisation
seurat_obj <- SCTransform(seurat_obj, assay = "RNA", new.assay.name = "SCT",
                          vars.to.regress = c("percent.mt", "percent.rb"), 
                          batch_var = "ARRAY", conserve.memory = TRUE)
saveRDS(seurat_obj, "/Volumes/LACIE/CROP-seq/ARRAY_AGGR_OBJ.rds")

# Run general visualisation and clustering - this is just in case we want
# to take a look at this
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
saveRDS(seurat_obj, "/Volumes/LACIE/CROP-seq/ARRAY_AGGR_OBJ.rds")