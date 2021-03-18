
#' Variable selection using recursive feature elimination.
#'
#' Compares random forests based on nested subsets of the variables and selects
#' those variables leading to the forest with the smallest prediction error within a tolerance.
#'
#' Note: This function differs from the approach implemented in the R package
#' \code{\link[varSelRF]{varSelRF}} because it recalculates importance scores in each step. The tolerance step is based on the
#' \code{pickSizeTolerance} function in the R package \code{caret}.
#'
#' @inheritParams wrapper.rf
#' @param prop.rm proportion of variables removed at each step (default value of \code{\link[varSelRF]{varSelRF}})
#' @param tol acceptable difference in optimal performance (finds the smallest subset size that has a percent loss less than tol)
#' @param recalculate logical stating if importance should be recalculated at each iteration (default: TRUE)
#'
#' @return List with the following components:
#'   \itemize{
#'   \item \code{info} data.frame
#'   with information for each variable
#'   \itemize{
#'   \item included.until.subset = number of smallest subset which contains variable
#'   \item selected = variable has been selected
#'   }
#'   \item \code{var} vector of selected variables
#'   \item \code{info.runs} data.frame with information for each run
#'   \itemize{
#'   \item n = number of variables
#'   \item mse = mean squared error
#'   \item rsq = R^2
#'   }}
#'
#'  @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # select variables
#' res = var.sel.rfe(x = data[, -1], y = data[, 1], prop.rm = 0.2, recalculate = TRUE)
#' res$var
#'
#' @export

var.sel.rfe <- function(x, y, prop.rm = 0.2, recalculate = TRUE, tol = 10, ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1,
                        no.threads = 1, method = "ranger", type = "regression", importance = "impurity_corrected", case.weights = NULL) {

  ## importance for nested subsets
  infos = NULL
  info.var = matrix(ncol = 0, nrow = ncol(x),
                    dimnames = list(colnames(x), NULL))

  var = colnames(x)
  imp = NULL
  while (length(var) >= 2) {
    print(paste(length(var), "variables"))

    ## train RF
    if (length(var) == 1) {
      x.sub = matrix(x[, var], ncol = 1)
    } else {
      x.sub = x[, var]
    }
    rf = wrapper.rf(x = x.sub, y = y,
                    ntree = ntree, mtry.prop = mtry.prop, nodesize.prop = nodesize.prop,
                    no.threads = no.threads,
                    method = method, type = type, importance = importance, case.weights = case.weights)

    if (is.null(imp) | recalculate) {
      #            print("recalculating importance ...")
      imp = get.vim(rf)
    }
    imp = sort(imp[var])
    #        print(paste("length imp:", length(imp)))

    ## save information
    error = calculate.error(rf = rf, true = y)
    infos = rbind(infos, c(n = length(var), error))
    temp = rep(0, nrow(info.var))
    temp[rownames(info.var) %in% var] = 1
    info.var = cbind(info.var, temp)

    ## remove variables
#    no.rm = ceiling(prop.rm * length(var))
    no.rm = round(prop.rm * length(var))
    var = names(imp)[-(1:no.rm)]
  }

  ## get number of variables with minimal error within tolerance
  ## (select smallest set if several RFs have smallest error)
  infos = data.frame(infos)
  if (type == "regression") {
    error.name = "rmse"
  } else {
    error.name = "err"
  }
  best = min(infos[, error.name])
  if (best == 0) {
    ind.min = max(which(infos[, error.name] == best))
  } else {
    error.prop = (infos[, error.name] - best)/best * 100
    ind.min = max(which(error.prop <= tol))
  }

  ## calculate smallest subset for each variable
  sum = apply(info.var, 1, sum)

  ## select variables
  ind.sel = as.numeric(info.var[, ind.min] == 1)

  ## info about variables
  info = data.frame(sum, ind.sel)
  colnames(info) = c("included.until.subset", "selected")
  return(list(info = info, var = sort(rownames(info)[info$selected == 1]), info.runs = infos))

}
