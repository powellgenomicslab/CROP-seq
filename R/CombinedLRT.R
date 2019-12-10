#!/usr/bin/env Rscript
# Scripts for calculating combined likelihood ratio
library(Seurat)
library(tidyverse)
library(BiocParallel)

# Load environment
MinMax <- function(data, min, max) {
  data2 <- data
  data2[data2 > max] <- max
  data2[data2 < min] <- min
  return(data2)
}

# Secondary function to address values that don't fit in the model
bimodLikData_FixNoVar <- function(x, xmin = 0) {
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  xal <- MinMax(
    data = length(x = x2) / length(x = x),
    min = 1e-5,
    max = (1 - 1e-5)
  )
  
  likA <- length(x = x1) * log(x = 1 - xal)
  
  #likelihood of positivec cells 
  likB <- length(x = x2) * log(x = xal)
  
  return(likA + likB)
}

#internal function to run mcdavid et al. DE test
#
#' @importFrom stats sd dnorm
#
bimodLikData <- function(x, xmin = 0) {
  # x1 and x2 are 2 vectors representing 2 modes
  # x1 for 0 values -> on/off distribution model 
  # x2 for positive values -> normal, continuous distribution 
  x1 <- x[x <= xmin]
  x2 <- x[x > xmin]
  
  # estimate proportion of positive cells 
  # use 1e-5 as min and 1-1e-5 as max (i.e. if there is only 1 nonzero among 100K cells)
  xal <- MinMax(
    data = length(x = x2) / length(x = x),
    min = 1e-5,
    max = (1 - 1e-5)
  )
  
  # likelihood for observing x1, 1-xal is expected ratio of 0 values 
  likA <- length(x = x1) * log(x = 1 - xal)
  
  # calculate variance for x2, to be used in dnorm to calculate prob distribution
  if (length(x = x2) < 2) {
    mysd <- 1
  } else {
    mysd <- sd(x = x2)
  }
  
  # Likelihood for observing x2
  likB <- length(x = x2) *
    log(x = xal) +
    sum(dnorm(x = x2, mean = mean(x = x2), sd = mysd, log = TRUE))
  return(likA + likB)
}

# x = counts for test gene
# y = counts for control gene
DifferentialLRT <- function(x, y, xmin = 0) {
  lrtX <- bimodLikData(x = x)
  lrtY <- bimodLikData(x = y)
  lrtZ <- bimodLikData(x = c(x, y))
  lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  
  # Check to account for results that do not conform to expected model
  if (is.infinite(lrt_diff) || (lrt_diff < 0) || is.nan(lrt_diff) || is.na(lrt_diff)){
    lrtX <- bimodLikData_FixNoVar(x = x) 
    lrtY <- bimodLikData_FixNoVar(x = y)
    lrtZ <- bimodLikData_FixNoVar(x = c(x, y))    
    lrt_diff <- 2 * (lrtX + lrtY - lrtZ)
  }
  
  return(pchisq(q = lrt_diff, df = 3, lower.tail = F))
}

# Wrapper to all the fancy stuff
runCombinedLRT <- function(gene, test_matrix = NULL, control_matrix = NULL){
  # Retrieve counts
  test_counts <- as.numeric(unlist(test_matrix[gene, ]))
  control_counts <- as.numeric(unlist(control_matrix[gene, ]))
  
  # Perform Combined LRT
  de.result <- DifferentialLRT(test_counts, control_counts)    
  return(de.result)
}

worker <- function(gene_list, test_matrix = NULL, control_matrix = NULL){
  # Calculate p-values
  p_vals <- lapply(gene_list, runCombinedLRT, test_matrix = test_matrix, control_matrix)
  names(p_vals) <- gene_list
  p_vals <- unlist(p_vals)
  padj_vals <- stats::p.adjust(p_vals, method = "BH")
  
  # Calculate mean expression
  test_means <- Matrix::rowMeans(test_matrix[gene_list, ])
  control_means <- Matrix::rowMeans(control_matrix[gene_list, ])
  
  # FoldChange
  foldChange <- test_means/control_means
  log2FoldChange <- log2(foldChange)
  
  output_df <- data.frame(gene_id = gene_list,
                          test_mean = test_means, 
                          control_mean = control_means, 
                          pval = p_vals, padj = padj_vals, 
                          foldChange = foldChange, 
                          log2FoldChange = log2FoldChange)
  return(output_df)  
}


# MC DAVID TEST
# Discrete/continuous model for single cell expression data based on a mixture 
# a point mass at zero and log-normal distribution
# Use log normal counts
# LRT 

target <- "NonTargeting-Human-0004-gene"
seurat_obj <- readRDS("/Volumes/LACIE/CROP-seq/ARRAY_AGGR_OBJ.rds")

# Retrieve metadata from Seurat object
metadata <- FetchData(seurat_obj, vars = c("ARRAY", "TYPE", "gRNA", "TARGET"))

# Check if we are testing a target - if so, exclude it from the control pool
if (startsWith(target, "NonTargeting")){
  control_guides <- grep("^NonTargeting", metadata$TARGET, value = TRUE)
  control_guides <- unique(control_guides[which(control_guides != target)])
  control_info <- subset(metadata, metadata$TARGET %in% control_guides)
} else{
  control_info <- subset(metadata, metadata$TARGET == "NonTargeting")
}

# Subset out groups that we are going to test
test_info <- subset(metadata, metadata$TARGET == target)


# Use UMI-corrected count data and split into two matrices
normcounts <- seurat_obj[["SCT"]]@counts
test_counts <- normcounts[ , rownames(test_info)]
control_counts <- normcounts[, rownames(control_info)]

# Remove guide-associated genes from query list
gene_names <- rownames(normcounts)[which(!(rownames(normcounts) %in% grep("-gene$", rownames(normcounts), value = TRUE)))]

chunked_genes <- split(gene_names, seq(1,4))

test_results <- lapply(chunked_genes, worker, test_matrix = test_counts, control_matrix = control_counts)
de_results <- do.call("rbind", test_results)
de_results <- de_results[order(abs(de_results$log2FoldChange), decreasing = TRUE), ]
de_results <- de_results %>% filter(!(is.infinite(log2FoldChange)))
de_results <- de_results %>% mutate(target = target)
de_results <- de_results %>% dplyr::select(gene_id, target, everything())
write_tsv(de_results, sprintf("%s_CROPseq_DE.tsv", target))
