#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
library(tidyr)
library(dplyr)
library(readr)
library(tibble)
library(Seurat)

# For pooled puro data data
input_dir <- args[1]
sample_name <- args[2]

grna_count_df <- read_tsv(input_dir)

summarized_reads <- grna_count_df %>% group_by(cell_barcode, gRNA) %>% summarize(total_reads = sum(reads))
summarized_reads <- summarized_reads %>% separate(gRNA, "_", into = c("type", "gRNA", "target"))
summarized_reads$gRNA[which(summarized_reads$type == "NonTargeting")] <- paste0("nt", summarized_reads$target[which(summarized_reads$type == "NonTargeting")])
summarized_reads$target[which(summarized_reads$type == "NonTargeting")] <- "NonTargeting"
summarized_reads <- ungroup(summarized_reads)

selectGuide <- function(x, y){
  if (nrow(x) > 1){
    # 1. If all of them are 1, we will not assign it something
    if (all(x["total_reads"] == 1)){
      output <- data.frame(cell_barcode = y, gRNA = NA)
    } 
    # 2. Otherwise, if all guides are equal, we will not assign it something
    else if (n_distinct(x["total_reads"]) == 1){
      output <- data.frame(cell_barcode = y, gRNA = NA)
    }
    # 3. If one is larger than the rest, and it is more than (sum(others) * 3) we will keep it
    else{
      # Find largest value and retrieve rows with largest value
      max_val <- max(x["total_reads"])
      max_row <- x[which(x["total_reads"] == max_val), ]
      
      # If there are multiple rows that are tied, check if they target the 
      # same gene
      # If they do, run standard check
      if (nrow(max_row) > 1){
        if (n_distinct(max_row["target"]) == 1){
          max_row <- max_row[1, ]
          other_rows <- x[which(!(x["gRNA"] %in% max_row["gRNA"])), ]
          if (max_row["total_reads"] > (sum(other_rows["total_reads"]))*3){
            output <- data.frame(cell_barcode = y, gRNA = max_row["gRNA"])
          } else{
            output <- data.frame(cell_barcode = y, gRNA = NA)
          }
        } else{
          output <- data.frame(cell_barcode = y, gRNA = NA)
        }        
      } else{
        # Find largest value and retrieve rows with largest value
        other_rows <- x[which(!(x["gRNA"] %in% max_row["gRNA"])), ]
        if (max_row["total_reads"] > (sum(other_rows["total_reads"]))*3){
          output <- data.frame(cell_barcode = y, gRNA = max_row["gRNA"])
        } else{
          output <- data.frame(cell_barcode = y, gRNA = NA)
        }
      }
  }} 
  else{
    output <- data.frame(cell_barcode = y, gRNA = x["gRNA"])
  }
  return(output)
}

# Remove null values
summarized_reads <- summarized_reads %>% group_by(cell_barcode)
assignment_list <- group_map(summarized_reads, selectGuide)
assignments_df <- bind_rows(assignment_list)
assignments_df <- left_join(assignments_df, summarized_reads[, c("cell_barcode", "type", "gRNA", "target")], by = c("cell_barcode", "gRNA"))
write_tsv(assignments_df, sprintf("%s_ReadAssignments.tsv", sample_name))
