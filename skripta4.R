library(data.table)
library(ggplot2)


GetPopulation <- function(type, probability, size) {
  # Computes the sample covariance between two vectors.
  #
  # Args:
  #   type: type of independent variable
  #   n: number of simulations
  #   prob: popultion probability of each independent variable
  #   verbose: If TRUE, prints sample covariance; if not, not. Default is TRUE.
  #
  # Returns:
  #   The sample covariance between x and y.
  
  if (type == 'stage') {
    values <- 1:length(probability)
    stage <- sample(values, size = size, prob = probability, replace = TRUE)
    population <- model.matrix(~as.factor(stage)-1)
    colnames(population) <- paste0('S', values)
  }
  
  if (type == 'complication') {
    values <- 1:length(probability)
    population <- t(replicate(size, rbinom(length(probability), 1, probability)))
    colnames(population) <- paste0('Z', values)
  }
  population
}

GetHypothesisData <- function(population, Ha = NULL) {
  
  #H0
  if (is.null(Ha)) {
    # z = 0 + c(0, 0, 0, 0) %*% t(population)
    # linear combination with bias
    z = 0 + rep(0, dim(population)[2]) %*% t(population)
  }
  # probability for each sample based on linearn combination
  Pr = 1/(1+exp(-z))
  
  y = rbinom(dim(population)[1], 1, Pr)
  
  data <- list(y = y, X = population)
}
 
size <- 10000

# STADIJ
probability  <- c(0.60, 0.25, 0.10, 0.05)


# OPERACIJA
probability <- c(0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.3, 0.4, 0.5, 0.6)

GetSimulation <- function(pop.matrix, n) {
  # Computes the sample covariance between two vectors.
  #
  # Args:
  #   k: population size
  #   n: number of simulations
  #   prob: popultion probability of each independent variable
  #   verbose: If TRUE, prints sample covariance; if not, not. Default is TRUE.
  #
  # Returns:
  #   The sample covariance between x and y.
  set.seed(24)
  
  
  
}

stadij.verjetnosti  <- c(0.60, 0.25, 0.10, 0.05)
