
#' Wrapper function to call random forests function.
#'
#' Provides an interface to different parallel implementations of the random
#' forest algorithm. Currently, only the \code{ranger} package is
#' supported.
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
#' @param method implementation to be used ("ranger").
#' @param type mode of prediction ("regression", "classification" or "probability").
#' @param importance Variable importance mode ('none', 'impurity',
#' 'impurity_corrected' or 'permutation'). Default is 'impurity_corrected'.
#' @param case.weights Weights for sampling of training observations. Observations with larger weights will be selected with higher probability in the bootstrap (or subsampled) samples for the trees.
#' @param ... further arguments needed for \code{\link[Pomona]{holdout.rf}} function only.
#'
#' @return An object of class \code{\link[ranger]{ranger}}.
#'
#' @import methods stats
#'
#' @export
#'
#' @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # regression
#' wrapper.rf(x = data[, -1], y = data[, 1],
#'            type = "regression", method = "ranger")

wrapper.rf <- function(x, y, ntree = 500,
                       mtry.prop = 0.2,
                       nodesize.prop = 0.1,
                       no.threads = 1,
                       method = "ranger",
                       type = "regression",
                       importance = "impurity_corrected",
                       case.weights = NULL, ...) {

  ## check data
  if (length(y) != nrow(x)) {
    stop("length of y and number of rows in x are different")
  }

  if (any(is.na(x))) {
    stop("missing values are not allowed")
  }

  if (type %in% c("probability", "regression") & (is.character(y) | is.factor(y))) {
    stop("only numeric y allowed for probability or regression mode")
  }

  ## set global parameters
  nodesize = floor(nodesize.prop * nrow(x))
  mtry = floor(mtry.prop * ncol(x))
  if (mtry == 0) mtry = 1

  if (type == "classification") {
    #    print("in classification")
    y = as.factor(y)
  }

  ## run RF
  if (method == "ranger") {
    if (type == "probability") {
      y = as.factor(y)
      prob = TRUE
    } else {
      prob = FALSE
    }

    rf = ranger::ranger(data = data.frame(y, x),
                        dependent.variable.name = "y",
                        probability = prob,
                        importance = importance,
                        scale.permutation.importance = FALSE,
                        num.trees = ntree,
                        mtry = mtry,
                        case.weights = case.weights,
                        min.node.size = nodesize,
                        num.threads = no.threads,
                        write.forest = TRUE,
                        ...)
  } else {
    stop(paste("method", method, "undefined. Use 'ranger'."))
  }

  return(rf)
}


#' Error calculation.
#'
#' Calculates errors by comparing predictions with the true values. For
#' regression and probability mode, it will give root mean squared error (rmse) and
#' pseudo R-squared (rsq). For classification mode, overall accuracy (acc), overall
#' error (err), Matthews correlation coefficient (mcc), sensitivity (sens) and
#' specificity (spec) are returned.
#'
#' @param rf Object of class \code{\link[ranger]{ranger}}
#' @param true vector with true value for each sample
#' @param test.set matrix or data.frame of predictor variables for test set with variables in
#'   columns and samples in rows (Note: missing values are not allowed)
#'
#' @return numeric vector with two elements for regression and probability estimation (rmse, rsq) and
#' five elements for classification (acc, err, mcc, sens, spec)
#'
#' @export
#'
#' @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # random forest
#' rf = wrapper.rf(x = data[, -1], y = data[, 1],
#'                 type = "regression")
#'
#' # error
#' calculate.error(rf = rf, true = data[, 1])

calculate.error <- function(rf, true, test.set = NULL) {

  if (is(rf, "ranger")) {
    if (!is.null(test.set)) {
      pred = predict(rf, data = test.set)$predictions
    } else {
      pred = rf$predictions
    }
    if (rf$treetype == "Probability estimation") {
      pred = pred[, 2]
    }
  } else {
    stop(paste("rf needs to be of class ranger"))
  }

    if ((is(rf, "randomForest") && rf$type == "classification") |
        (is(rf, "ranger") && rf$treetype == "Classification")) {
        conf.matrix = table(pred = pred, true = true)
        tp = conf.matrix[2, 2]
        tn = conf.matrix[1, 1]
        fn = conf.matrix[2, 1]
        fp = conf.matrix[1, 2]

        ## accuracy
        acc = (tp + tn) / sum(conf.matrix)

        ## Matthews correlation coefficient
        mcc = (tp * tn - fp * fn) /
          sqrt( (tp + fn) * (tn + fp) * (tp + fp) * (tn + fn))

        ## sensitivity
        sens = tp / (tp + fn)

        ## specificity
        spec = tn / (fp + tn)

        error = c(err = 1 - acc, acc = acc, mcc = mcc, sens = sens, spec = spec)
    } else {
      mse = sum((pred - true)^2, na.rm = TRUE) / sum(!is.na(pred))

      ## pseudo R-squared uses sum of squared differences divided by n instead of variance!
      v = sum((true - mean(true))^2) / length(true)
      rsq = 1 - mse/v
      error = c(rmse = sqrt(mse), rsq = rsq)
    }

    return(error)
}


#' Get variable importance.
#'
#' Extracts variable importance depending on class of random forest object.
#'
#' @param rf Object of class \code{\link[ranger]{ranger}}
#'
#' @return numeric vector with importance value for each variable (in original order)
#'
#' @export

get.vim <- function(rf) {
  if (is(rf, "ranger")) {
    vim = ranger::importance(rf)
  } else {
    stop(paste("rf needs to be of class ranger"))
  }
  return(vim)
}
