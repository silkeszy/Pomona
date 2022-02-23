
#' Variable selection using hybrid function.
#'
#' Variable selection using the hybrid function and combining Vita and hybrid.
#'
#' This function selects only variables that are confirmed based on hybrid implementation.
#' Note that this function uses the ranger implementation for variable selection.
#'

#' @inheritParams wrapper.rf
#' @param pValue confidence level (default: 0.01 based on \code{\link[Boruta]{Boruta}} package)
#' @param maxRuns maximal number of importance source runs (default: 100 based on \code{\link[Boruta]{Boruta}} package)
#' @param alpha significance threshold used by Vita
#' @param seed Seed
#'
#' @return List with the following components:
#'   \itemize{
#'   \item \code{info} data.frame
#'   with information of each variable
#'   \itemize{
#'   \item run.x = original variable importance (VIM) in run x
#'   (includes min, mean and max of VIM of shadow variables)
#'   \item decision = Boruta decision (Confirmed, Rejected or Tentative)
#'   \item selected = variable has been selected
#'   }
#'   \item \code{var} vector of selected variables
#'   \item \code{info.shadow.var} data.frame with information about
#'   minimal, mean and maximal shadow variables of each run
#'   }
#'
#'  @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # select variables
#' res = var.sel.hybrid(x = data[, -1], y = data[, 1])
#' res$var
#'
#' @export


var.sel.hybrid <- function(x, y,
                           pValue = 0.01,
                           maxRuns = 100,
                           alpha = 0.05,
                           ntree = 500,
                           mtry.prop = 0.2,
                           nodesize.prop = 0.1,
                           no.threads = 1,
                           method = "ranger",
                           type = "regression",
                           importance = "impurity_corrected",
                           seed = 123,
                           case.weights = NULL) {

  ## variable selection using Boruta function
  ## ----------------------------------------
  # modified version of getImpRfRaw function to enable user defined mtry
  # values
  get_imp_ranger <- function(x, y, ...){
    x <- data.frame(x)
    imp <- wrapper.rf(x = x,
                      y = y,
                      ...)$variable.importance
    return(imp)
  }
  res.boruta = hybrid(x = x, y = y,
                      pValue = pValue,
                      maxRuns = maxRuns,
                      alpha = alpha,
                      ntree = ntree,
                      nodesize.prop = nodesize.prop,
                      no.threads = no.threads,
                      mtry.prop = mtry.prop,
                      getImp = get_imp_ranger,
                      importance = importance,
                      seed = seed,
                      type = type,
                      case.weights = case.weights)

  ## select variables
  dec = res.boruta$finalDecision
  ind.sel = rep(0, ncol(x))
  ind.sel[dec == "Confirmed"] = 1
  info.sel = data.frame(decision = dec, selected = ind.sel)

  #   ## info about variables
  #   info.var = t(res.boruta$ImpHistory)
  #   colnames(info.var) = paste("run", 1:ncol(info.var), sep = ".")
  #   info.var = merge(info.var, info.sel, all.x = TRUE, by.x = "row.names", by.y = "row.names")
  #   rownames(info.var) = info.var[, 1]
  #   info.var = info.var[, -1]
  #   info.var = info.var[c(colnames(x), "shadowMax", "shadowMean", "shadowMin"), ]
  #
  #   ind.shadow = grep("shadow", rownames(info.var))
  #   return(list(info = info.var[-ind.shadow, ],
  #               var = sort(rownames(info.var)[info.var$selected == 1]),
  #               info.shadow.var = info.var[ind.shadow, -which(colnames(info.var) %in% c("decision", "selected"))]))

  ## info about variables
  info.var = t(res.boruta$ImpHistory)
  colnames(info.var) = paste("run", 1:ncol(info.var), sep = ".")
  info.shadow.var = info.var[grep("shadow", rownames(info.var)),]
  info.var = info.var[-grep("shadow", rownames(info.var)),]
  if (all.equal(rownames(info.var), rownames(info.sel))) {
    info.var = cbind(info.var, info.sel)
  } else {
    info.var = merge(info.var, info.sel, by.x = "row.names", by.y = "row.names")
  }

  return(list(info = info.var,
              var = sort(rownames(info.var)[info.var$selected == 1]),
              info.shadow.var = info.shadow.var))

}
