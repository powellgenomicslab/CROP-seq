#!/usr/bin/env Rscript
# Scripts for calculating combined likelihood ratio

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
  likB <- length(x = x2) *
    log(x = xal)
  
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
