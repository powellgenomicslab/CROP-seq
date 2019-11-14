args = commandArgs(trailingOnly=TRUE)
library(readr)
library(tidyr)
library(dplyr)
library(purrr)

sample_name <- args[1]
transcriptome_df <- read_tsv(sprintf("/Users/anne/Dropbox (Garvan)/Anne/CROP-seq/TranscriptAssignments/%s_TranscriptAssignments.tsv",sample_name))
read_df <- read_tsv(sprintf("/Users/anne/Dropbox (Garvan)/Anne/CROP-seq/ReadAssignments/%s_ReadAssignments.tsv", sample_name))
transcriptome_df <- separate(transcriptome_df, gRNA, "-", into = c("type", "gRNA", "gene", "feature"))
transcriptome_df[which(transcriptome_df$gRNA == "Human"), "gRNA"] <- paste0("nt",transcriptome_df$gene[which(transcriptome_df$gRNA == "Human")])
transcriptome_df[which(transcriptome_df$type == "NonTargeting"), "gene"] <- "NonTargeting"
transcriptome_df <- transcriptome_df[, c("cell_barcode", "type", "gRNA", "gene")]

# Combine assignments
combined_df <- left_join(transcriptome_df, read_df, by = "cell_barcode", suffix = c(".UMI", ".read"))

# Filter conflicts and all NA values
filtered_combined_df <- combined_df %>% filter_at(vars(gRNA.UMI, gRNA.read), any_vars(!is.na(.)))
filtered_combined_df <- filtered_combined_df %>% filter_at(vars(gene, target), any_vars(is.na(gene) | is.na(target) | (gene == target)))

mergeResult <- function(x){
  if (any(is.na(x))){
    if (all(is.na(x[2:4]))){
      x <- x %>% select(cell_barcode = cell_barcode, type = type.read, gRNA = gRNA.read, target = target)
    } else{
      x <- x %>% select(cell_barcode = cell_barcode, type = type.UMI, gRNA = gRNA.UMI, target = gene)
    }
  } else{
    x <- x %>% select(cell_barcode = cell_barcode, type = type.UMI, gRNA = gRNA.UMI, target = gene)
  }
  return(x)
}

filtered_rows <- lapply(1:nrow(filtered_combined_df), function(x) mergeResult(filtered_combined_df[x, ]))
final_assignments <- bind_rows(filtered_rows)
write_tsv(final_assignments, sprintf("/Users/anne/Dropbox (Garvan)/Anne/CROP-seq/FinalCellAssignments/%s_FinalAssignments.tsv",sample_name))
