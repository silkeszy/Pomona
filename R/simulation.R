
#' Simulate gene expression data.
#'
#' Generates simulated gene expression data in a regression
#' setting.
#'
#' The underlying simulation model is described in detail in the paper XXX. In
#' brief, a nonlinear regression model based on three uniformly distributed
#' variables is used. Predictor variables are simulated to be correlated with one
#' of those functional variables. In addition, independent, uniformly distributed
#' predictor variables are simulated.
#'
#' @param no.samples number of samples.
#' @param group.size number of variables in each of the six groups of correlated 
#' variables.
#' @param no.var.total total number of variables. 
#' @param null simulate null model (using independent functional variables).
#' 
#' @return A data.frame with samples in rows and variables in columns (Note: first
#'   column contains simulated phenotype). Variables are named as y (= phenotype),
#'   g.i.j (= variable j in group i) and ind.k (= k-th independent variable).
#'
#' @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'                      
#' @export                      

simulation.data.cor <- function(no.samples, group.size, no.var.total, null = FALSE) {
  
  if (length(group.size) != 6) {
    stop("group.size needs to be a vector of length 6!")
  }
  ## simulate functional variables
  x = sapply(1:6, function(i) {
    #    rnorm(n = no.samples, mean = 0, sd = 1)})
    runif(n = no.samples, min = 0, max = 1)})
  colnames(x) = paste("x", 1:6, sep = ".")
  
  ## simulate groups of correlated variables based on functional variables
  g = lapply(1:3, function(i) {
    v = sapply(1:group.size[i], function(j) {
      x[, i] + 0.01 + 0.5 * (j - 1) / (group.size[i] - 1)  *
        rnorm(n = no.samples, mean = 0, sd = 0.3)
    })
    colnames(v) = paste("g", i, 1:group.size[i], sep = ".")
    return(v)
  })
  g = do.call(cbind, g)
  
  ## simulate groups of correlated variables based on null variables
  g.0 = lapply(1:3, function(i) {
    #     v = sapply(1:group.size[i], function(j) {
    #       x[, i+3] + 0.2 * rnorm(n = no.samples, mean = 0, sd = 0.2)
    #     })
    v = sapply(1:group.size[i], function(j) {
      x[, i+3] + 0.01 + 0.5 * (j - 1) / (group.size[i] - 1)  *
        rnorm(n = no.samples, mean = 0, sd = 0.3)
      
    })
    colnames(v) = paste("g.0", i, 1:group.size[i], sep = ".")
    return(v)
  })
  g.0 = do.call(cbind, g.0)
  
  if (ncol(g) + ncol(g.0) >= no.var.total) {
    warning("No additional independent variables are simulated!")
    ind = NULL
  } else {
    ind = matrix(runif(no.samples * (no.var.total - (ncol(g) + ncol(g.0))),
                       min = 0, max = 1),
                 nrow = no.samples)
    colnames(ind) = paste("ind", 1:ncol(ind), sep = ".")
  }
  
  ## phenotype
  if (null) {
    x.use = sapply(1:6, function(i) {
      #    rnorm(n = no.samples, mean = 0, sd = 1)})
      runif(n = no.samples, min = 0, max = 1)})
    
  } else {
    x.use = x[, 1:3]  
  }
  
  #  y = 0.25 * exp(4 * x[, 1]) + 4 / (1 + exp(-20 * (x[, 2] - 0.5))) + 3 * x[, 3] +
  #    rnorm(n = no.samples, mean = 0, sd = 0.2)
  y = 0.25 * exp(4 * x.use[, 1]) + 4 / (1 + exp(-20 * (x.use[, 2] - 0.5))) + 3 * x.use[, 3] +
    rnorm(n = no.samples, mean = 0, sd = 0.2)
  
  ## data
  return(data.frame(y = y, g, g.0, ind))
  
}
