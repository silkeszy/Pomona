
#' Variable selection using recurrent relative variable importance (r2VIM).
#'
#' Generates several random forests using all variables and different random
#' number seeds. For each run, the importance score is divided by the (absolute)
#' minimal importance score (relative importance scores). Variables are selected
#' if the minimal relative importance score is >= factor.
#'
#' Note: This function is a reimplementation of the R package \code{RFVarSelGWAS}.
#'
#' @inheritParams wrapper.rf
#' @param no.runs number of random forests to be generated
#' @param factor minimal relative importance score for a variable to be selected
#'
#' @return List with the following components:
#'   \itemize{
#'   \item \code{info} data.frame
#'   with information for each variable
#'   \itemize{
#'   \item vim.run.x = original variable importance (VIM) in run x
#'   \item rel.vim.run.x = relative VIM in run x
#'   \item rel.vim.min = minimal relative VIM over all runs
#'   \item rel.vim.med = median relative VIM over all runs
#'   \item selected = variable has been selected
#'   }
#'   \item \code{var} vector of selected variables
#'   }
#'
#'  @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # select variables
#' res = var.sel.r2vim(x = data[, -1], y = data[, 1], no.runs = 5, factor = 1)
#' res$var
#'
#' @export

var.sel.r2vim <- function(x, y, no.runs = 10, factor = 1, ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1,
                          no.threads = 1, method = "ranger", type = "regression", importance = "impurity_corrected",
                          case.weights = NULL) {

  ## importance for each run
  imp.all = NULL
  for (r in 1:no.runs) {
    print(paste("run", r))
    rf = wrapper.rf(x = x, y = y,
                    ntree = ntree, mtry.prop = mtry.prop, nodesize.prop = nodesize.prop, no.threads = no.threads,
                    method = method, type = type, importance = importance, case.weights = case.weights)
    imp.all = cbind(imp.all, get.vim(rf))
  }

  ## factors
  min.global = min(imp.all)
  if (min.global >= 0) {
    stop("Global minimal importance score is not negative!")
  }
  no.neg.min = 0
  fac = matrix(nrow = nrow(imp.all), ncol = ncol(imp.all),
               dimnames = dimnames(imp.all))
  for (i in 1:ncol(imp.all)) {
    x = imp.all[,i]
    min = min(x)
    if (min >= 0) {
      no.neg.min = no.neg.min + 1
      fac[,i] = x / abs(min.global)
    } else {
      fac[, i] = x / abs(min)
    }
  }
  if (no.neg.min > 0) {
    print(paste(no.neg.min, "runs with no negative importance score!"))
  }
  fac.min = apply(fac, 1, min)
  fac.med = apply(fac, 1, median)

  ## select variables
  ind.sel = as.numeric(fac.min >= factor)

  ## info about variables
  info = data.frame(imp.all, fac, fac.min, fac.med, ind.sel)
  colnames(info) = c(paste("vim.run.", 1:no.runs, sep = ""),
                     paste("rel.vim.run.", 1:no.runs, sep = ""),
                     "rel.vim.min", "rel.vim.median", "selected")
  return(list(info = info, var = sort(rownames(info)[info$selected == 1])))
}
