library(Seurat)
library(readr)
library(tidyr)
library(dplyr)

sample_name <- "POOLED2"
input_path <- sprintf("/Volumes/LACIE/CROP-seq/%s_SeuratObject.rds", sample_name)
seurat_obj <- readRDS(input_path)
metadata <- seurat_obj@meta.data
metadata$cell_barcode <- rownames(metadata)
metadata <- metadata %>% mutate(sample = sample_name) %>% select(cell_barcode, sample, type, gRNA, target_gene) 
write_tsv(metadata, sprintf("/Users/anne/Dropbox (Garvan)/Anne/CROP-seq/FinalCellAssignments/%s_PostFilterAssignments.tsv", sample_name))

# REad all in
file_list <- Sys.glob("/Users/anne/Dropbox (Garvan)/Anne/CROP-seq/FinalCellAssignments/*_PostFilterAssignments.tsv")
df_list <- lapply(file_list, read_tsv)
assignment_df <- bind_rows(df_list)

assignment_df