#' Weighted Quadratic Discriminant Analysis
#'
#' @author Solveig Pospiech, package 'MASS'
# #' Raimon Tolosana-Delgado, K. Gerald v.d. Boogaart
#'
#' @description Extension of the qda() of package 'MASS' to calculate a QDA incorporating individual, cell-wise uncertainties, e.g. variances per measuring point.
#'
#' @details Uncertainties can be considered in a statistical analysis either by each measured variable, by each observation or by using the individual, cell-wise uncertainties.
#' There are several methods for incorporating variable-wise or observation-wise uncertainties into a QDA, most of them using the uncertainties as weights for the variables or observations of the data set.
#' The term 'cell-wise uncertainties' describe a data set of $d$ analysed variables where each observation has an individual uncertainty for each of the $d$ variables conforming it.
#' Hence, a data set of $n \\times d$ data values has associated a data set of $n \\times d$ individual uncertainties.
#' Instead of weighting the columns or rows of the data set, the vqda() function uses uncertainties to recalculate better estimates of the group variances and group means.
#' If the presence of uncertainties is not accounted for, the decision rules  are based on the group variances calculated by the given data set.
#' But this observed group variance might deviate notably from the group variance, which can be estimated including the uncertainties.
#' This methodological framework does not only allow to incorporate cell-wise uncertainties, but also would largely be valid if the information about the co-dependency between uncertainties within each observation would be reported.
#'
#' @param x data frame or matrix
#' @param uncertainties data frame or matrix, values for uncertainties per cell. Uncertainties should be relative errors, e.g. relative standard deviation
#' @param grouping a factor specifying the group for each observation.
#' @param prior the prior probabilities of class membership. If unspecified, the class proportions for the training set are used. If present, the probabilities should be specified in the order of the factor levels.
#'
#' @examples
#' # for non-compositional data:
#' data("dataobs")
#' data("uncertainties")
#' myqda = vqda(x = dataobs[, 1:2], uncertainties = uncertainties[, 1:2], grouping = dataobs$Group)
#' mypred = predict(myqda, newdata = dataobs[, 1:2], newerror = uncertainties[, 1:2])
#' forplot = cbind(dataobs, LG1 = mypred$posterior[,1])
#' if (require("ggplot2")) {
#'   scatter_plot = ggplot(data = forplot, aes(x = Var1, y = Var2)) +
#'     geom_point(aes(shape = Group, color = LG1))
#'   if (require("ggthemes")) {
#'     scatter_plot = scatter_plot +
#'         scale_color_gradientn(colours = colorblind_pal()(5))
#'   }
#'   scatter_plot
#' }
#'
#' # for compositional data
#' data("dataobs_coda")
#' data("uncertainties_coda")
#' require(compositions)
#' # generate ilr-transformation
#' data_ilr = ilr(dataobs_coda[, 1:3])
#' uncert_ilr = t(simplify2array(apply(uncertainties_coda[, 1:3],1,
#'                        function(Delta) clrvar2ilr(diag(Delta)))))
#' attr(x = uncert_ilr, which = "class") <- "rmult"  # correct the class
#' myqda_coda = vqda(x = data_ilr, uncertainties = uncert_ilr, grouping = dataobs_coda$Group)
#' mypred_coda = predict(myqda_coda, newdata = data_ilr, newerror = uncert_ilr)
#' forplot_coda = cbind(dataobs_coda, LG1 = mypred_coda$posterior[,1])
#' # if 'ggtern' is installed, you can plot via ggtern:
#' # if (require("ggtern")) {
#' #   ternary_plot = ggtern(data = forplot_coda, aes(x = Var1, y = Var2, z = Var3)) +
#' #     geom_point(aes(shape = Group, color = LG1))
#' #   if (require("ggthemes")) {
#' #     ternary_plot = ternary_plot +
#' #         scale_color_gradientn(colours = colorblind_pal()(5))
#' #   }
#' #   ternary_plot
#' # }
#'
#' @return object of class 'vqda' containing the following components:
#' \code{prior} the prior probabilities used.
#' \code{counts} counts per group.
#' \code{means} the group means.
#' \code{generalizedMeans} the group means calculated by the function \code{\link{generalized_mean}}
#' \code{groupVarCorrected} the group variances calculated by the function \code{\link{calc_estimate_true_var}}
#' \code{lev} the levels of the grouping factor.
#' \code{grouping} the factor specifying the class for each observation.
#'
#' @export
vqda <- function(x,
                 uncertainties,
                 grouping,
                 prior = proportions) {

  # missing: formula, ....
  # check if uncertainties and x have same class
  if (class(x) != class(uncertainties)) {
    if ("rmult" %in% class(x) & !"rmult" %in% class(uncertainties)) warning("x has class'rmult' but uncertainties has not. Are you sure this is correct?")
    if (!"rmult" %in% class(x) & "rmult" %in% class(uncertainties)) warning("uncertainties has class'rmult' but x has not. Are you sure this is correct?")
    if (!"rmult" %in% class(x) & !"rmult" %in% class(uncertainties)) warning("x and uncertainties have different classes. This might cause problems at a later stage.")
  }

  # prepare data ------------------------------------------------------------
  # copied from MASS::qda.default
  if (is.null(dim(x)))
    stop("'x' is not a matrix")
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite, NA or NaN values in 'x'")
  n <- nrow(x)
  p <- ncol(x)
  # --
  if (class(uncertainties) == 'rmult') {
    p = p^2 # because then the uncertainties are stored as arrays in the rows.
    if (ncol(uncertainties) != p) stop("transformed uncertainties are expected to be the outcome of 't(apply( <original_uncertainties> , 1, function(Delta) clrvar2ilr(diag(Delta))))'. Please make sure the input of the transformed uncertainties is correct.")
  }
  # --
  if (n != length(grouping))
    stop("nrow(x) and length(grouping) are different")
  g <- as.factor(grouping)
  lev <- levels(g)
  counts <- as.vector(table(g))
  names(counts) <- lev
  if (any(counts < p + 1))
    stop("some group is too small for 'qda'")
  proportions <- counts/length(g)
  ng <- length(proportions)
  if (any(prior < 0) || round(sum(prior), 5) != 1)
    stop("invalid 'prior'")
  if (length(prior) != ng)
    stop("'prior' is of incorrect length")
  names(prior) <- lev
  group.means <- tapply(unclass(x), list(rep(g, ncol(x)), col(x)), mean)
  # scaling <- array(dim = c(p, p, ng))
  # ldet <- numeric(ng)

  # copy end

  # own code --------

  Zg = split(x, grouping) # split works perfectly fine for rmult-class to get in the lists matrix.
  # if class is not rmult, split needs a data.frame to put it into matrix-like entries in the list-entries
  # to avoid the costs of data.frame but also to avoid that rmult is mandatory for running this function, here comes the reforming of the list- entries:
  if (!"rmult" %in% class(x)) Zg = lapply(Zg, function(y) matrix(y, ncol = ncol(x)))
  # split errors into the groups, now they are vectors
  sigmaIg <- split(compositions::rmult(uncertainties), grouping)
  if (!"rmult" %in% class(uncertainties)) sigmaIg = lapply(sigmaIg, function(y) matrix(y, ncol = ncol(uncertainties))) # for this line see comment for Zg
  # generate a new sigma: the sigmasums have to be normalized to nrow(per group) and subtracted:
  sigmacorrected <- mapply(calc_estimate_true_var, Zg, sigmaIg, SIMPLIFY = FALSE)
  meancorrected <- mapply(generalized_mean, Zg, sigmacorrected, sigmaIg, SIMPLIFY = FALSE)
  structure(list(
    prior = prior,
    counts = counts,
    means = group.means,
    generalizedMeans = meancorrected,
    groupVarCorrected = sigmacorrected,
    lev = lev,
    grouping = grouping
  )
  ,class = "vqda")
}



#' Weighted Linear Discriminant Analysis
#'
#' @author Solveig Pospiech, package 'MASS'
# #' Raimon Tolosana-Delgado, K. Gerald v.d. Boogaart
#'
#' @description Extension of the qda() of package 'MASS' (not the lda() function) to calculate a LDA incorporating individual, cell-wise uncertainties, e.g. variances per measuring point.
#'
#' @details Uncertainties can be considered in a statistical analysis either by each measured variable, by each observation or by using the individual, cell-wise uncertainties.
#' There are several methods for incorporating variable-wise or observation-wise uncertainties into a QDA, most of them using the uncertainties as weights for the variables or observations of the data set.
#' The term 'cell-wise uncertainties' describe a data set of $d$ analysed variables where each observation has an individual uncertainty for each of the $d$ variables conforming it.
#' Hence, a data set of $n \\times d$ data values has associated a data set of $n \\times d$ individual uncertainties.
#' Instead of weighting the columns or rows of the data set, the vlda() function uses uncertainties to recalculate better estimates of the group variances and group means.
#' It is internally very similar to the vqda() function, but with an averaged group variance for all groups.
#' If the presence of uncertainties is not accounted for, the decision rules  are based on the group variances calculated by the given data set.
#' But this observed group variance might deviate notably from the group variance, which can be estimated including the uncertainties.
#' This methodological framework does not only allow to incorporate cell-wise uncertainties, but also would largely be valid if the information about the co-dependency between uncertainties within each observation would be reported.
#'
#' @param x data frame or matrix
#' @param uncertainties data frame or matrix, values for uncertainties per cell. Uncertainties should be relative errors, e.g. relative standard deviation#'
#' @param grouping a factor specifying the group for each observation.
#' @param prior the prior probabilities of class membership. If unspecified, the class proportions for the training set are used. If present, the probabilities should be specified in the order of the factor levels.
#'
#' @examples
#' # for non-compositional data:
#' data("dataobs")
#' data("uncertainties")
#' mylda = vlda(x = dataobs[, 1:2], uncertainties = uncertainties[, 1:2], grouping = dataobs$Group)
#' mypred = predict(mylda, newdata = dataobs[, 1:2], newerror = uncertainties[, 1:2])
#' forplot = cbind(dataobs, LG1 = mypred$posterior[,1])
#' if (require("ggplot2")) {
#'   scatter_plot = ggplot(data = forplot, aes(x = Var1, y = Var2)) +
#'     geom_point(aes(shape = Group, color = LG1))
#'   if (require("ggthemes")) {
#'     scatter_plot = scatter_plot +
#'         scale_color_gradientn(colours = colorblind_pal()(5))
#'   }
#'   scatter_plot
#' }
#'
#' # for compositional data
#' data("dataobs_coda")
#' data("uncertainties_coda")
#' require(compositions)
#' # generate ilr-transformation
#' data_ilr = ilr(dataobs_coda[, 1:3])
#' uncert_ilr = t(simplify2array(apply(uncertainties_coda[, 1:3],1,
#'                        function(Delta) clrvar2ilr(diag(Delta)))))
#' attr(x = uncert_ilr, which = "class") <- "rmult"  # correct the class
#' mylda_coda = vlda(x = data_ilr, uncertainties = uncert_ilr, grouping = dataobs_coda$Group)
#' mypred_coda = predict(mylda_coda, newdata = data_ilr, newerror = uncert_ilr)
#' forplot_coda = cbind(dataobs_coda, LG1 = mypred_coda$posterior[,1])
#' # if 'ggtern' is installed, you can plot via ggtern:
#' # if (require("ggtern")) {
#' #   ternary_plot = ggtern(data = forplot_coda, aes(x = Var1, y = Var2, z = Var3)) +
#' #     geom_point(aes(shape = Group, color = LG1))
#' #   if (require("ggthemes")) {
#' #     ternary_plot = ternary_plot +
#' #         scale_color_gradientn(colours = colorblind_pal()(5))
#' #   }
#' #   ternary_plot
#' # }
#'
#' @return object of class 'vlda' containing the following components:
#' \code{prior} the prior probabilities used.
#' \code{counts} counts per group.
#' \code{means} the group means.
#' \code{generalizedMeans} the group means calculated by the function \code{\link{generalized_mean}}
#' \code{groupVarCorrected} the group variances calculated by the function \code{\link{calc_estimate_true_var}}
#' \code{lev} the levels of the grouping factor.
#' \code{grouping} the factor specifying the class for each observation.
#'
#' @export
vlda <- function(x,
                 uncertainties,
                 grouping,
                 prior = proportions) {

  # missing: formula, ....
  # check if uncertainties and x have same class
  if (class(x) != class(uncertainties)) {
    if ("rmult" %in% class(x) & !"rmult" %in% class(uncertainties)) warning("x has class'rmult' but uncertainties has not. Are you sure this is correct?")
    if (!"rmult" %in% class(x) & "rmult" %in% class(uncertainties)) warning("uncertainties has class'rmult' but x has not. Are you sure this is correct?")
    if (!"rmult" %in% class(x) & !"rmult" %in% class(uncertainties)) warning("x and uncertainties have different classes. This might cause problems at a later stage.")
  }

  # prepare data ------------------------------------------------------------
  # copied from MASS::qda.default
  if (is.null(dim(x)))
    stop("'x' is not a matrix")
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite, NA or NaN values in 'x'")
  n <- nrow(x)
  p <- ncol(x)
  # --
  if (class(uncertainties) == 'rmult') {
    p = p^2 # because then the uncertainties are stored as arrays in the rows.
    if (ncol(uncertainties) != p) stop("transformed uncertainties are expected to be the outcome of 't(apply( <original_uncertainties> , 1, function(Delta) clrvar2ilr(diag(Delta))))'. Please make sure the input of the transformed uncertainties is correct.")
  }
  # --
  if (n != length(grouping))
    stop("nrow(x) and length(grouping) are different")
  g <- as.factor(grouping)
  lev <- levels(g)
  counts <- as.vector(table(g))
  names(counts) <- lev
  if (any(counts < p + 1))
    stop("some group is too small for 'lda'")
  proportions <- counts/length(g)
  ng <- length(proportions)
  if (any(prior < 0) || round(sum(prior), 5) != 1)
    stop("invalid 'prior'")
  if (length(prior) != ng)
    stop("'prior' is of incorrect length")
  names(prior) <- lev
  group.means <- tapply(unclass(x), list(rep(g, ncol(x)), col(x)), mean)
  # scaling <- array(dim = c(p, p, ng))
  # ldet <- numeric(ng)

  # copy end

  # own code --------

  Zg = split(x, grouping) # split works perfectly fine for rmult-class to get in the lists matrix.
  # if class is not rmult, split needs a data.frame to put it into matrix-like entries in the list-entries
  # to avoid the costs of data.frame but also to avoid that rmult is mandatory for running this function, here comes the reforming of the list- entries:
  if (!"rmult" %in% class(x)) Zg = lapply(Zg, function(y) matrix(y, ncol = ncol(x)))
  # split errors into the groups, now they are vectors
  sigmaIg <- split(compositions::rmult(uncertainties), grouping)
  if (!"rmult" %in% class(uncertainties)) sigmaIg = lapply(sigmaIg, function(y) matrix(y, ncol = ncol(uncertainties))) # for this line see comment for Zg
  # sum up all group variances and normalize by DF (maybe not the cleanest code, because the var is still there)
  averaged_variance = Reduce("+", lapply(Zg, function(y) compositions::var(y)*(nrow(y) - 1)))/(n - ng)
  # sum all uncertainties by group
  averaged_uncertainties = Reduce("+", lapply(sigmaIg, colSums))/n
  message("Checking positive definiteness of corrected variance for all groups...")
  if ("rmult" %in% class(uncertainties)) {
    sigmacorrected_t <- force_posdef(averaged_variance - matrix(averaged_uncertainties, ncol = ncol(averaged_variance)))
  } else {
    sigmacorrected_t <- force_posdef(averaged_variance - diag(averaged_uncertainties))
  }
  # generate a list with the sigmacorrected for each group, to mimick the group variances of QDA
  sigmacorrected = rep(list(sigmacorrected_t), ng)
  meancorrected <- mapply(generalized_mean, Zg, sigmacorrected, sigmaIg, SIMPLIFY = FALSE)
  structure(list(
    prior = prior,
    counts = counts,
    means = group.means,
    generalizedMeans = meancorrected,
    groupVarCorrected = sigmacorrected,
    lev = lev,
    grouping = grouping
  )
  ,class = "vlda")
}


#' predict.vqda
#'
#' @author Solveig Pospiech, package 'MASS'
#'
#' @description Classify multivariate observations in conjunction with qda() or lda() of class 'vqda' or 'vlda'.
#'
#' @param object object of class 'vqda' or 'vlda'.
#' @param ... additional arguments affecting the predictions produced.
#'
# predict <- function(object, ...) {
#   UseMethod("predict", object)
# }
#' @describeIn predict predict() for class 'vqda'
#' @param newdata data frame or matrix of cases to be classified or, if object has a formula, a data frame with columns of the same names as the variables used.
#' A vector will be interpreted as a row vector. If newdata is missing, an attempt will be made to retrieve the data used to fit the qda object.
#' @param newerror data frame or matrix of uncertainties corresponding to the cases in 'newdata'.
#' @param prior the prior probabilities of group membership. If unspecified, the prior of the object are used.
#' @param ... ...
#'
#' @return list containing the following components:
#' \code{class} factor containing the predicted group
#' \code{likelihood} matrix of dimension 'number of samples' x 'number of groups', containing the likelihood for each sample to belong to one of the groups
#' \code{grouping} original grouping of the samples, copied from the input object
#'
#' @export
predict.vqda <- function(object,
                                 newdata,
                                 newerror,
                                 prior = object$prior, ...) {

  # copied and slightly adjusted from predict.qda
  if (!inherits(object, "vqda"))
    stop("object not of class \"vqda\"")

  ngroup <- length(object$prior)
  if (!missing(prior)) {
    if (any(prior < 0) || round(sum(prior), 5) != 1)
      stop("invalid 'prior'")
    if (length(prior) != ngroup)
      stop("'prior' is of incorrect length")
  }

  if (!is.null(Terms <- object$terms)) {
    if (missing(newdata))
      newdata <- model.frame(object)
    else {
      newdata <- model.frame(as.formula(delete.response(Terms)),
                             newdata, na.action = function(x) x, xlev = object$xlevels)
    }
    x <- model.matrix(delete.response(Terms), newdata, contrasts = object$contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0L)
    if (xint > 0)
      x <- x[, -xint, drop = FALSE]
    # if (method == "looCV")
    #   g <- model.response(newdata)
  }
  else {
    if (missing(newdata)) {
      if (!is.null(sub <- object$call$subset)) {
        newdata <- eval.parent(parse(text = paste(deparse(object$call$x,
                                                          backtick = TRUE), "[", deparse(sub, backtick = TRUE),
                                                  ",]")))
        g <- eval.parent(parse(text = paste(deparse(object$call[[3L]],
                                                    backtick = TRUE), "[", deparse(sub, backtick = TRUE),
                                            "]")))
      }
      else {
        newdata <- eval.parent(object$call$x)
        g <- eval.parent(object$call[[3L]])
      }
      if (!is.null(nas <- object$call$na.action)) {
        df <- data.frame(g = g, X = newdata)
        df <- eval(call(nas, df))
        g <- df$g
        newdata <- df$X
      }
      g <- as.factor(g)
    }
    if (is.null(dim(newdata)))
      dim(newdata) <- c(1, length(newdata))
    x <- as.matrix(newdata)
  }

  # own code ----

  # check if uncertainties and x have same class
  if (class(newdata) != class(newerror)) {
    if ("rmult" %in% class(newdata) & !"rmult" %in% class(newerror)) warning("newdata has class'rmult' but newerror has not. Are you sure this is correct?")
    if (!"rmult" %in% class(newdata) & "rmult" %in% class(newerror)) warning("newerror has class'rmult' but newdata has not. Are you sure this is correct?")
    if (!"rmult" %in% class(newdata) & !"rmult" %in% class(newerror)) warning("newdata and newerror have different classes. This might cause problems at a later stage.")
  }

  newSigmaIs <- if ("rmult" %in% class(newerror)) lapply(1:nrow(newerror), function(i) matrix(newerror[i,], ncol = ncol(newdata)))
  else lapply(1:nrow(newerror), function(i) diag(newerror[i,]))

  newZs <- lapply(1:nrow(newdata), function(i) newdata[i,]) # packt alle Beobachtungen in Listen

  barmug = lapply(object$generalizedMeans, c)
  groupVarCorrected = object$groupVarCorrected

  LgsFunc <- function(Z0, Sigma0) {
    AuxFunc <- function(mug, barSigma) {
      SigmaSimp = barSigma + Sigma0
      -0.5*determinant(SigmaSimp, logarithm = T)$modulus - 0.5*c(as.numeric(Z0 - mug) %*% solve(SigmaSimp) %*% as.numeric(Z0 - mug))
    }
    mapply(AuxFunc, barmug, groupVarCorrected)
  }

  L <- mapply(LgsFunc, newZs, newSigmaIs) # Ergebnis von dieser Funktion sollte matrix mit [n x Anzahl der Gruppen]

  p = t(exp(L))
  posterior = compositions::clo(p)
  cl <- factor(max.col(posterior), levels = seq_along(object$lev),
               labels = object$lev)
  dimnames(posterior) <- list(rownames(x), object$lev)
  #structure( # not used at the moment
  return(list(class = cl,
       posterior = posterior,
       likelihood = L,
       grouping = object$grouping)
  #,class = "predicted_weighted")
  )
}

#' @describeIn predict predict() for class 'vlda'
#' @param newdata data frame or matrix of cases to be classified or, if object has a formula, a data frame with columns of the same names as the variables used.
#' A vector will be interpreted as a row vector. If newdata is missing, an attempt will be made to retrieve the data used to fit the qda object.
#' @param newerror data frame or matrix of uncertainties corresponding to the cases in 'newdata'.
#' @param prior the prior probabilities of group membership. If unspecified, the prior of the object are used.
#' @param ... ...
#'
#' @return list containing the following components:
#' \code{class} factor containing the predicted group
#' \code{likelihood} matrix of dimension 'number of samples' x 'number of groups', containing the likelihood for each sample to belong to one of the groups
#' \code{grouping} original grouping of the samples, copied from the input object
#'
#' @export
predict.vlda <- function(object,
                                 newdata,
                                 newerror,
                                 prior = object$prior, ...) {

  # copied and slightly adjusted from predict.qda
  if (!inherits(object, "vlda"))
    stop("object not of class \"vlda\"")

  ngroup <- length(object$prior)
  if (!missing(prior)) {
    if (any(prior < 0) || round(sum(prior), 5) != 1)
      stop("invalid 'prior'")
    if (length(prior) != ngroup)
      stop("'prior' is of incorrect length")
  }

  if (!is.null(Terms <- object$terms)) {
    if (missing(newdata))
      newdata <- model.frame(object)
    else {
      newdata <- model.frame(as.formula(delete.response(Terms)),
                             newdata, na.action = function(x) x, xlev = object$xlevels)
    }
    x <- model.matrix(delete.response(Terms), newdata, contrasts = object$contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0L)
    if (xint > 0)
      x <- x[, -xint, drop = FALSE]
    # if (method == "looCV")
    #   g <- model.response(newdata)
  }
  else {
    if (missing(newdata)) {
      if (!is.null(sub <- object$call$subset)) {
        newdata <- eval.parent(parse(text = paste(deparse(object$call$x,
                                                          backtick = TRUE), "[", deparse(sub, backtick = TRUE),
                                                  ",]")))
        g <- eval.parent(parse(text = paste(deparse(object$call[[3L]],
                                                    backtick = TRUE), "[", deparse(sub, backtick = TRUE),
                                            "]")))
      }
      else {
        newdata <- eval.parent(object$call$x)
        g <- eval.parent(object$call[[3L]])
      }
      if (!is.null(nas <- object$call$na.action)) {
        df <- data.frame(g = g, X = newdata)
        df <- eval(call(nas, df))
        g <- df$g
        newdata <- df$X
      }
      g <- as.factor(g)
    }
    if (is.null(dim(newdata)))
      dim(newdata) <- c(1, length(newdata))
    x <- as.matrix(newdata)
  }

  # own code ----

  # check if uncertainties and x have same class
  if (class(newdata) != class(newerror)) {
    if ("rmult" %in% class(newdata) & !"rmult" %in% class(newerror)) warning("newdata has class'rmult' but newerror has not. Are you sure this is correct?")
    if (!"rmult" %in% class(newdata) & "rmult" %in% class(newerror)) warning("newerror has class'rmult' but newdata has not. Are you sure this is correct?")
    if (!"rmult" %in% class(newdata) & !"rmult" %in% class(newerror)) warning("newdata and newerror have different classes. This might cause problems at a later stage.")
  }

  newSigmaIs <- if ("rmult" %in% class(newerror)) lapply(1:nrow(newerror), function(i) matrix(newerror[i,], ncol = ncol(newdata)))
  else lapply(1:nrow(newerror), function(i) diag(newerror[i,]))

  newZs <- lapply(1:nrow(newdata), function(i) newdata[i,]) # packt alle Beobachtungen in Listen

  barmug = lapply(object$generalizedMeans, c)
  groupVarCorrected = object$groupVarCorrected

  LgsFunc <- function(Z0, Sigma0) {
    AuxFunc <- function(mug, barSigma) {
      SigmaSimp = barSigma + Sigma0
      as.numeric(Z0) %*% solve(SigmaSimp) %*% mug - 0.5*c(mug %*% solve(SigmaSimp) %*% mug)
    }
    mapply(AuxFunc, barmug, groupVarCorrected)
  }

  L <- mapply(LgsFunc, newZs, newSigmaIs) # Ergebnis von dieser Funktion sollte matrix mit [n x Anzahl der Gruppen]

  p = t(exp(L)) * object$prior
  posterior = compositions::clo(p)
  cl <- factor(max.col(posterior), levels = seq_along(object$lev),
               labels = object$lev)
  dimnames(posterior) <- list(rownames(x), object$lev)
  #structure( # not used at the moment
  return(list(class = cl,
       posterior = posterior,
       likelihood = L,
       grouping = object$grouping)
  #,class = "predicted_weighted")
  )
}
