
#' Variable selection using a permutation approach.
#'
#' Selects variables which importance scores are larger than scores calculated
#' after permuting the phenotype. Output is a p-value calculated as the
#' proportion of permutations with an equal or larger importance score.
#'
#' Note:
#' This function is a reimplementation of the approach in the R package \code{rfPermute} and the
#' parametric permutation approach by Altmann et al. (2010).
#'
#' @inheritParams wrapper.rf
#' @param no.perm number of permutations
#' @param p.t threshold for p-values (all variables with a p-value = 0  or < p.t will be selected)
#' @param parametric logical stating if parametric permutation approach of Altmann et al. 2010 (based on normal
#' distribution) should be used (default: FALSE)
#'
#' @return List with the following components:
#'   \itemize{
#'   \item \code{info} data.frame
#'   with information for each variable
#'   \itemize{
#'   \item vim.original = original variable importance (VIM)
#'   \item vim.perm.x = VIM in permutation x
#'   \item pvalue = proportion of permutations with a larger VIM than original VIM (nonparametric) or probability of observing the original or a larger VIM,
#'   given the fitted null importance distribution based on normal distributions (parametric)
#'   \item selected = variable has been selected
#'   }
#'   \item \code{var} vector of selected variables
#'   }
#'
#' @references
#'   Altmann, A., Tolosi, L., Sander, O. & Lengauer, T. (2010). Permutation importance: a corrected feature importance measure, Bioinformatics 26:1340-1347.
#'
#'  @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # select variables based on nonparametric permutation approach
#' res = var.sel.perm(x = data[, -1], y = data[, 1], no.perm = 10, p.t = 0)
#' res$var
#'
#' # select variables based on parametric permutation approach
#' res.par = var.sel.perm(x = data[, -1], y = data[, 1], no.perm = 10, p.t = 0.05, parametric = TRUE)
#' res.par$var
#'
#' @export

var.sel.perm <- function(x, y, no.perm = 100, p.t = 0, ntree = 500, mtry.prop = 0.2, nodesize.prop = 0.1,
                         no.threads = 1, method = "ranger", type = "regression",
                         parametric = FALSE, importance = "impurity_corrected", case.weights = NULL) {

  ## importance for original data
  rf = wrapper.rf(x = x, y = y,
                  ntree = ntree, mtry.prop = mtry.prop, nodesize.prop = nodesize.prop, no.threads = no.threads,
                  method = method, type = type, importance = importance, case.weights = case.weights)
  imp = get.vim(rf)

  ## importance for permutated data
  imp.perm = sapply(1:no.perm, function(i) {
    if (i %% 50 == 0) print(paste("permutation", i))

    ## permute phenotype
    y.perm = sample(y)

    ## RF for permuted data
    rf.perm = wrapper.rf(x = x, y = y.perm,
                         ntree = ntree, mtry.prop = mtry.prop, nodesize.prop = nodesize.prop,
                         no.threads = no.threads,
                         method = method, type = type, importance = importance)
    get.vim(rf.perm)
  })

  ## calculate p-value
  if (!parametric) {
    pval = sapply(1:nrow(imp.perm), function(i) {
      sum(imp.perm[i,] >= imp[i]) / no.perm})
  } else {
    # estimate mean and variance
    mean = apply(imp.perm, 1, mean)
    sd = apply(imp.perm, 1, sd)

    # p-value based on normal distribution
    pval = sapply(1:nrow(imp.perm), function(i) {
      pnorm(imp[i], mean[i], sd[i], lower.tail = FALSE)})
  }
  names(pval) = rownames(imp)

  ## select variables
  ind.sel = as.numeric(pval == 0 | pval < p.t)

  ## info about variables
  info = data.frame(imp, imp.perm, pval, ind.sel)
  colnames(info) = c("vim.original",
                     paste("vim.perm.", 1:no.perm, sep = ""),
                     "pvalue", "selected")
  return(list(info = info, var = sort(rownames(info)[info$selected == 1])))
}


