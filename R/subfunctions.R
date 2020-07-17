#' Generalized mean
#'
#' @author Solveig Pospiech, K. Gerald v.d. Boogaart
#'
#' @description Calculates the generalized mean of a data set by using a given group variance and individual, observation-wise variances for each observation of the data set
#'
#' @param x a matrix containing the data for which the mean should be calculated
#' @param ... not implemented
#'
#' @return vector of lenght of ncol(x) of generalized means
#'
#' @export
generalized_mean <- function(x, ...) {
  UseMethod("generalized_mean", x)
}

#' @describeIn generalized_mean for class matrix or data.frame
#' @param var a matrix containing the corrected (estimated true) variance of the data set
#' @param individual_var a matrix containing individual variances. Default is a 0 - matrix with the dimensions of x, can be used for implementing the individual uncertainties of each observation
#' @export
generalized_mean.default <- function(x, var, individual_var = matrix(0, nrow = nrow(x), ncol = ncol(x)), ...) {
  # these two functions are equivalent in all generalized_mean functionos and to gsi.Diag resp. gsi.Inv of the pkg 'compositions'.
  # they are here newly defined because until now gsi-functions from 'compositions' cannot be used outside the package
  gsiDiag <- function(d) {
    if (length(d) > 1)
      return(diag(d))
    return( structure(d,dim = c(1,1)))
  }
  gsiInv <- function(A,tol=1E-15) {
    with(svd(A),v %*% gsiDiag(ifelse(abs(d)/max(d) > tol,1/d,0)) %*% t(u))
  }

  if (is.null(dim(x)))
    stop("'x' is not a matrix")
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite, NA or NaN values in 'x'")
  if (is.null(dim(individual_var)))
    stop("'individual_var' is not a matrix")
  individual_var <- as.matrix(individual_var)
  if (any(!is.finite(individual_var)))
    stop("infinite, NA or NaN values in 'individual_var'")
  if (dim(individual_var)[2] == 1) { # if only one variable exists:
    sigmasum <- apply(individual_var, 1, function(s) var + s)
    sigmaInv = 1/sigmasum
    normalization = sum(sigmaInv)
    aux = sapply(1:nrow(x), function(i) sigmaInv[i] * x[i,] )
    aux_sum = sum(aux)
  } else {# for more than one variable
    sigmasum <- apply(individual_var, 1, function(s) force_posdef(var + diag(s))) # diag because uncertainties are expected to have only entries in the diagonal
    sigmaInv <- lapply(1:ncol(sigmasum), function(i) gsiInv(matrix(sigmasum[,i], ncol = ncol(x))))
    zw <- unlist(sigmaInv)
    zw[is.na(zw)] <- 0
    normalization = matrix(rowSums(matrix(zw, ncol = nrow(x), byrow = F)), ncol = ncol(x))
    aux = sapply(1:nrow(x), function(i) sigmaInv[[i]] %*% x[i,] )
    aux_sum = rowSums(aux)
  }
  # erg<-svdInv(normalisation)  %*% tensorA::margin(mapply(invMul,Sigmas,xi),1) # wird anders programmiert
  erg <- solve(normalization) %*% aux_sum
  attr(erg,"Sigma") <- normalization
  return(erg)
}


#' @describeIn generalized_mean for class rmult of package 'compositions'
#' @param var a matrix containing the corrected (estimated true) group variances
#' @param individual_var a matrix containing individual variances. Default is a 0 - matrix with the dimensions of x, can be used for implementing the individual uncertainties
#' @export
generalized_mean.rmult <- function(x, var, individual_var = matrix(0, nrow = nrow(x), ncol = ncol(x)^2), ...) {
  # these two functions are equivalent in all generalized_mean functionos and to gsi.Diag resp. gsi.Inv of the pkg 'compositions'.
  # they are here newly defined because until now gsi-functions from 'compositions' cannot be used outside the package
  gsiDiag <- function(d) {
    if (length(d) > 1)
      return(diag(d))
    return( structure(d,dim = c(1,1)))
  }
  gsiInv <- function(A,tol=1E-15) {
    with(svd(A),v %*% gsiDiag(ifelse(abs(d)/max(d) > tol,1/d,0)) %*% t(u))
  }

  if (is.null(dim(x)))
    stop("'x' is not a matrix")
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite, NA or NaN values in 'x'")
  if (is.null(dim(individual_var)))
    stop("'individual_var' is not a matrix")
  individual_var <- as.matrix(individual_var)
  if (ncol(individual_var) != ncol(x)^2) stop("The individual variances seem to have not the right format.
                                                Please make sure that your uncertainties are also converted into the respective log-ratio space and that all entries of the resulting variances-covariance matrix are in one row.")
  if (any(!is.finite(individual_var)))
    stop("infinite, NA or NaN values in 'individual_var'")
  if (dim(individual_var)[2] == 1) { # if only one variable exists:
    sigmasum <- apply(individual_var, 1, function(s) var + s)
    sigmaInv = 1/sigmasum
    normalization = sum(sigmaInv)
    aux = sapply(1:nrow(x), function(i) sigmaInv[i] * as.vector(x[i,] )) # has to be "as.vector" because x is still rmult class
    aux_sum = sum(aux)
  } else {# for more than one variable
    sigmasum <- apply(individual_var, 1, function(s) force_posdef(var + matrix(s, ncol = ncol(x)))) # matrix because the uncertainties are expected to also have covariance entries
    sigmaInv <- lapply(1:ncol(sigmasum), function(i) gsiInv(matrix(sigmasum[,i], ncol = ncol(x))))
    zw <- unlist(sigmaInv)
    zw[is.na(zw)] <- 0
    normalization = matrix(rowSums(matrix(zw, ncol = nrow(x), byrow = F)), ncol = ncol(x))
    aux = sapply(1:nrow(x), function(i) sigmaInv[[i]] %*% as.vector(x[i,] )) # has to be "as.vector" because x is still rmult class
    aux_sum = rowSums(aux)
  }
  erg <- gsiInv(normalization) %*% aux_sum
  attr(erg,"Sigma") <- normalization
  return(erg)
}

#' Estimate true group variance
#'
#' @author Solveig Pospiech, K. Gerald v.d. Boogaart
#'
#' @description Estimation of true group variance incorporating observation wise variances.
#' The function uses the data from x and the individual variances for each observation, for example derived from uncertainties, to calculate a 'true' group variance.
#' The variance of the matrix is corrected for the sum of the individual variances of the data set, which is normalized to the number of rows of the matrix.
#'
#' @param x a matrix of data
#' @param ... ...
#'
#' @return matrix of corrected group variance
#'
#' @export
calc_estimate_true_var <- function(x, ...) {
  UseMethod("calc_estimate_true_var", x)
}

#' @describeIn calc_estimate_true_var for class matrix or data.frame
#' @param individual_var a matrix of cell-wise uncertainties, corresponding to the entries of 'x'
#' @param force_pos_def force positive definiteness of the new group variances, default TRUE
#' @export
calc_estimate_true_var.default <- function(x,
                                          individual_var,
                                          force_pos_def = T, ...) {
  if (is.null(dim(x)))
    stop("'x' is not a matrix")
  if (is.null(dim(individual_var)))
    stop("'individual_var' is not a matrix")
  # average the individual individual_var:
  averaged_sigmas = colSums(individual_var)/nrow(x)
  newgroupsigma = compositions::var(x) - diag(averaged_sigmas) # diag because uncertainties are expected to have only entries in the diagonal
  if (force_pos_def) {
    message("Checking positive definiteness of corrected group variances...")
    newgroupsigma = force_posdef(newgroupsigma)
    message("Checking done")
  }
  return(newgroupsigma)
}


#' @describeIn calc_estimate_true_var for class rmult
#' @param individual_var a matrix of cell-wise uncertainties, corresponding to the entries of 'x'
#' @param force_pos_def force positive definiteness of the new group variances, default TRUE
#' @export
calc_estimate_true_var.rmult <- function(x, individual_var, force_pos_def = T, ...) {
  if (is.null(dim(x)))
    stop("'x' is not a matrix")
  if (is.null(dim(individual_var)))
    stop("'individual_var' is not a matrix")
  # average the individual individual_var:
  averaged_sigmas = colSums(individual_var)/nrow(x)
  newgroupsigma = compositions::var(x) - matrix(averaged_sigmas, ncol = ncol(x)) # matrix because the uncertainties are expected to also have covariance entries
  if (force_pos_def) {
    message("Checking positive definiteness of corrected group variances...")
    newgroupsigma = force_posdef(newgroupsigma)
    message("Checking done")
  }
  return(newgroupsigma)
}

#' Force positive definiteness
#'
#' @description Function to force positive definiteness on a matrix.
#'
#' @author Solveig Pospiech
#'
#' @param x matrix
#' @param verbose logical, default TRUE. Should the function print the corrected eigenvalues?
#'
#' @return positive definite matrix
#'
# #' @export
force_posdef <- function(x, verbose = T) {
  # test for positive definiteness
  myeigen = eigen(x)
  if (sum(zw <- myeigen$values < 0) > 0 ) {
    message("Matrix is not positive definite. Eigenvalues are forced to non-negativeness.")
    if (verbose) print.noquote(paste("The eigenvalues:", myeigen$values))
    myeigen$values[zw] <- 0
    x = myeigen$vectors %*% diag(myeigen$values) %*% t(myeigen$vectors)
  }
  return(x)
}

