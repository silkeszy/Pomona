
#' Variable selection using Vita approach.
#'
#' This function calculates p-values based on the empirical null distribution from non-positive VIMs as
#' described in Janitza et al. (2015). Note that this function uses the \code{importance_pvalues} function in the R package
#' \code{\link[ranger]{ranger}}.
#'
#' @inheritParams wrapper.rf
#' @param p.t threshold for p-values (all variables with a p-value = 0  or < p.t will be selected)
#'
#' @return List with the following components:
#'   \itemize{
#'   \item \code{info} data.frame
#'   with information for each variable
#'   \itemize{
#'   \item vim = variable importance (VIM)
#'   \item CI_lower = lower confidence interval boundary
#'   \item CI_upper = upper confidence interval boundary
#'   \item pvalue = empirical p-value
#'   \item selected = variable has been selected
#'   }
#'   \item \code{var} vector of selected variables
#'   }
#'  @references
#'  Janitza, S., Celik, E. & Boulesteix, A.-L., (2015). A computationally fast variable importance test for random forest for high dimensional data, Technical Report 185, University of Munich, https://epub.ub.uni-muenchen.de/25587.
#'
#'  @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 500)
#'
#' # select variables
#' res = var.sel.vita(x = data[, -1], y = data[, 1], p.t = 0.05)
#' res$var
#'
#' @export

var.sel.vita <- function(x, y, p.t = 0.05,
                            ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1,
                            no.threads = 1, method = "ranger", type = "regression", importance = "impurity_corrected") {

  ## train holdout RFs
  res.holdout = holdout.rf(x = x, y = y,
                           ntree = ntree, mtry.prop = mtry.prop, nodesize.prop = nodesize.prop,
                           no.threads = no.threads, type = type,  importance = importance)

  ## variable selection using importance_pvalues function
  res.janitza = ranger::importance_pvalues(x = res.holdout,
                                           method = "janitza",
                                           conf.level = 0.95)
  res.janitza = as.data.frame(res.janitza)
  colnames(res.janitza)[1] = "vim"

  ## select variables
  ind.sel = as.numeric(res.janitza$pvalue == 0 | res.janitza$pvalue < p.t)

  ## info about variables
  info = data.frame(res.janitza, selected = ind.sel)
  return(list(info = info, var = sort(rownames(info)[info$selected == 1])))
}

#' Helper function for variable selection using Vita approach.
#'
#' This function calculates a modified version of the permutation importance using two cross-validation folds (holdout folds)
#' as described in Janitza et al. (2015). Note that this function is a reimplementation of the \code{holdoutRF} function in the
#' R package \code{\link[ranger]{ranger}}.
#'
#' @param x matrix or data.frame of predictor variables with variables in
#'   columns and samples in rows (Note: missing values are not allowed).
#' @param y vector with values of phenotype variable (Note: will be converted to factor if
#'   classification mode is used).
#' @param ntree number of trees.
#' @param mtry.prop proportion of variables that should be used at each split.
#' @param nodesize.prop proportion of minimal number of samples in terminal
#'   nodes.
#' @param no.threads number of threads used for parallel execution.
#' @param type mode of prediction ("regression", "classification" or "probability").
#' @param importance See \code{\link[ranger]{ranger}} for details.
#'
#' @return Hold-out random forests with variable importance
#'
#' @references
#' Janitza, S., Celik, E. & Boulesteix, A.-L., (2015). A computationally fast variable importance test for random forest for high dimensional data, Technical Report 185, University of Munich, https://epub.ub.uni-muenchen.de/25587.

holdout.rf <- function(x, y, ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1, no.threads = 1,
                       type = "regression", importance = importance) {

  ## define two cross-validation folds
  n = nrow(x)
  weights = rbinom(n, 1, 0.5)

  ## train two RFs
  res = list(rf1 = wrapper.rf(x = x, y = y,
                              ntree = ntree, mtry.prop = mtry.prop,
                              nodesize.prop = nodesize.prop, no.threads = no.threads,
                              method = "ranger", type = type,
                              case.weights = weights, replace = FALSE,
                              holdout = TRUE, importance = importance),
             rf2 = wrapper.rf(x = x, y = y,
                              ntree = ntree, mtry.prop = mtry.prop,
                              nodesize.prop = nodesize.prop, no.threads = no.threads,
                              method = "ranger", type = type,
                              case.weights = 1 - weights, replace = FALSE,
                              holdout = TRUE,  importance = importance))

  ## calculate mean VIM
  res$variable.importance = (res$rf1$variable.importance +
                               res$rf2$variable.importance)/2
  res$treetype = res$rf1$treetype
  res$importance.mode = res$rf1$importance.mode
  class(res) = "holdoutRF"
  return(res)
}
