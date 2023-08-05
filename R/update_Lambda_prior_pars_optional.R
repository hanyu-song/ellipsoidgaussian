#' Update Phi, a parameter in the hyperprior of Lambda
#'
#' @importFrom GIGrvg rgig
updatePhi <- function(a, lambda1) {
  temp <- abs(lambda1);
  ts <- sapply(c(temp), function(x) GIGrvg::rgig(n=1,  a - 1, 2 * x, 1))
#  ts <- sapply(c(temp), function(x) .Call("rgig",n=1,  a - 1, 2 * x, 1, PACKAGE = 'GIGrvg'))
  ts <- matrix(ts, nrow = dim(lambda1)[1],byrow = F)
  res <- sweep(ts, 1, Matrix::rowSums(ts), '/')
  return(res)
}
#' Update Kap, a parameter in the hyperprior of Lambda
#' @importFrom GIGrvg rgig
updateKap <- function(a, phi1, lambda1) {
  k <- dim(phi1)[2]; p <- dim(phi1)[1]
  para3 <- rowSums(abs(lambda1) / phi1)
 res <- sapply(1:p, function(j) GIGrvg::rgig(n = 1, k * a - k,  2 * para3[j],1))
#  res <- sapply(1:p, function(j) .Call("rgig",n = 1, k * a - k,  2 * para3[j],1, PACKAGE = 'GIGrvg'))
# res <- sapply(1:p, function(j) GIGrvg::rgig(n = 1, 1 - k, 1, 2 * para3[j]))
  return(res)
}
#' Update Psi, a parameter in the hyperprior of Lambda
#' @importFrom statmod rinvgauss
updatePsi <- function(kappa1, phi1, lambda1) {
 para <-  sweep(phi1, 1, kappa1, FUN = '*')
 para <- para / abs(lambda1)
  k <- ncol(para); p <- nrow(para)
 #res <- apply(para, 2, function(x) statmod::rinvgauss(p, 1, x))
  res <- apply(para, 2, function(x) statmod::rinvgauss(p, mean = x, shape = 1))
  res <- 1 / res
  return(res)
}

