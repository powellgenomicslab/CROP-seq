library(tidyr)
library(dplyr)
library(readr)
library(tibble)

grna_count_df <- read_csv("Pooled1_gRNA_UMI_CountTable.csv")
summarized_reads <- grna_count_df %>% group_by(cell_barcode, gRNA) %>% summarize(total_reads = sum(reads))
summarized_reads <- ungroup(summarized_reads)

selectGuide <- function(x, y){
  if (nrow(x) > 1){
    # 1. If all of them are 1, we will not assign it something
    # 2. Otherwise, if all guides are equal, we will not assign it something
    # 3. If one is larger than the rest, and it is more than (sum(others) * 3) we will keep it
    if (all(x["total_reads"] == 1)){
      return(NULL)
    } else if(length(unique(x["total_reads"])) == 1){
      return(NULL)
    } else{
      max_val <- max(x["total_reads"])
      max_guide <- x["gRNA"][which(x["total_reads"] == max_val)]
      other_val <- sum(x["total_reads"][which(x["gRNA"] != max_guide)])
      
      if(max_val > (other_val*3)){
        output <- data.frame(gRNA = max_guide, cell_barcode = y)
        return(output) 
      } else{
        return(NULL)
      }
    }
  } else{
    output <- data.frame(gRNA = x["gRNA"], cell_barcode = y)
    return(output)
  }
}

# Remove null values


