#' Update Phi, a parameter in the hyperprior of Lambda
#'
#' @description
#' `updatePhi` updates Phi in the Dirichlet-Laplace prior of Lambda.
#'
#' @param a The parameter in the Dirichlet-Laplace prior for Lambda,
#' typically set at 0.5
#' @param lambda1 A matrix of p by k, the factor loading matrix
#'
#' @returns updated Phi parameter

#' @importFrom GIGrvg rgig
updatePhi <- function(a, lambda1) {
  temp <- abs(lambda1);
  ts <- sapply(c(temp), function(x) GIGrvg::rgig(n=1,  lambda = a - 1, chi = 2 * x, psi = 1))
 # ts <- sapply(c(temp), function(x) rmutil::rginvgauss(n = 1, m = sqrt(2 * x), s = 0.5 / x, f = a - 1))
#  ts <- sapply(c(temp), function(x) .Call("rgig",n=1,  a - 1, 2 * x, 1, PACKAGE = 'GIGrvg'))
  ts <- matrix(ts, nrow = dim(lambda1)[1],byrow = FALSE)
  res <- sweep(ts, 1, Matrix::rowSums(ts), '/')
  return(res)
}

#' Update Kap, a parameter in the hyperprior of Lambda
#'
#' @description
#' `updateKap` updates the kappa parameter in the Dirichlet-Laplace prior of Lambda.
#'
#'
#' @param a The parameter in the Dirichlet-Laplace prior for Lambda.
#' @param phi1 The parameter phi in the Dirichlet-Laplace prior.
#' @param lambda1 A matrix of p by k, the factor loading matrix.
#'
#' @returns Updated Kap parameter
#' @importFrom GIGrvg rgig
updateKap <- function(a, phi1, lambda1) {
  k <- dim(phi1)[2]; p <- dim(phi1)[1]
  para3 <- rowSums(abs(lambda1) / phi1)
 # res <- sapply(1:p, function(j) rmutil::rginvgauss(n = 1, m = sqrt(2 * para3[j]),
         #                                           s = 0.5 / para3[j],f = k * a - k))
 res <- sapply(1:p, function(j) GIGrvg::rgig(n = 1, lambda = k * a - k,  chi = 2 * para3[j],psi = 1))
#  res <- sapply(1:p, function(j) .Call("rgig",n = 1, k * a - k,  2 * para3[j],1, PACKAGE = 'GIGrvg'))

  return(res)
}

#' Update Psi, a parameter in the hyperprior of Lambda
#'
#' @description
#' `updatePsi` updates the psi used in the Dirichlet-Laplace prior of Lambda.
#'
#' @param kappa1 The kappa parameter in the Dirichlet-Laplace prior of Lambda.
#' @param phi1 The phi parameter in the Dirichlet-Laplace prior of Lambda.
#' @param lambda1 A matrix of p by k, the factor loading matrix.
#'
#' @returns Updated Psi parameter
#' @importFrom GIGrvg rgig
updatePsi <- function(kappa1, phi1, lambda1) {
 para <-  sweep(phi1, 1, kappa1, FUN = '*')
 para <- para / abs(lambda1)
  k <- ncol(para); p <- nrow(para)

 # res <- apply(para, 2, function(x) statmod::rinvgauss(p, mean = x, shape = 1))
#  res <- apply(para, 2, function(x) rmutil::rinvgauss(n = p, m = x, s = 1))
  res <- apply(para, 2, function(x) GIGrvg::rgig(n = p, lambda = -0.5, chi = 1, psi = 1 / (x^2)))
  res <- 1 / res
  return(res)
}

