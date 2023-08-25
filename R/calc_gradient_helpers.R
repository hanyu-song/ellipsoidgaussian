### log likelihood derivative w.r.t mu, center and Lambda - for speed up purpose
#' Gradient w.r.t the first parameter of a pseudo-normalising constant
#'
#' @description
#' `calclogPseudoconst_VecParGrad_fast` calculates the gradient of the
#' log-likelihood w.r.t. the vector parameter of the pseudo-normalising constant
#'
#' @param vec_par The vector parameter value, a vector of length k.
#' @param eigVals Eigenvalues of the matrix parameter, a vector of length k.
#' @param eigVecs Eigenvectors of the matrix parameter, a matrix.
#' @param ord_eigVal XXX
#' @param order_  The order of saddlepoint approximation, 1: first-order, 2:
#' second-order first version, 3: second-order second version. Default is 3 for
#' better accuracy. See \insertCite{KumeWood05}{ellipsoidgaussian}
#' P975 for more details.
#'
#' @importFrom numDeriv grad
#' @importFrom matrixcalc vech
calclogPseudoconst_VecParGrad_fast <- function(vec_par,
                                               eigVals,
                                               eigVecs,
                                               ord_eigVal,
                                               order_ = 3) {
  # vec_par: length k
############## These three lines are replaced by the two lines below them ######
  k <- length(vec_par)
  para1 <- matrix(vec_par, nrow = k, ncol = 1)
  Vec <- t(eigVecs) %*% para1
#################
##############  These two lines are replaced by the 3 lines above them ######
 # Vec2 <- sweep(eigVecs, MARGIN = 1, STATS = vec_par, FUN = '*')
 # Vec <- colSums_zheyuan(Vec2, 2)
######################3
  Vec <- Vec[ord_eigVal]
 # res <- findFBconst(Vec, eigVals,order_, ordered = T)
  res <- findFBconst_cpp(Vec, eigVals,order_, ordered = T)
 # res <- res[order_]
  return(res)
}




#' Derivative w.r.t the transformed tau parameter
#'
#' @description
#' `calc_transTau_grad` calculates the derivative w.r.t. the transformed tau
#' parameter, namely the parameter that is being directly updated in the sampling
#' procedure.
#'
#' @param tau The value of tau.
#' @param vec_grad_sum XXXX
#' @param mu A vector of length k.
#' @param tauBounds A list of tau bounds, with names \code{lowerB} and
#' \code{upperB}.
#'
#' @returns The log-likelihood derivative w.r.t tau
calc_transTau_grad <- function(tau, vec_grad_sum, mu, tauBounds) {
  grad_tau <- sum(vec_grad_sum * mu);
  # change more
  if (is.null(tauBounds$upperB)) {
    return(logTransGradAdjust(tau, grad_tau, lowerB = tauBounds$lowerB))
  }
  res <- grad_tau * dtau_dlogit(tau, tauBounds$lowerB, tauBounds$upperB)
  res <- res / tauBounds$fac
  return(res)
}

#' Gradient w.r.t the first parameter in the pseudo-normalising constant
#'
#' @description
#' `calcVecGrad_all` calculates the gradient w.r.t. the first parameter in the
#' pseudo-normalising constant for all the data in the minibatch.
#'
#' @param para1 The first parameter in the pseudo-normalising constant,
#' a matrix of k by minibatch size.
#' @param ord_eigVals XXXX
#' @param eigVecs The eigenvectors of the second (or matrix) parameter in the
#' pseudo-normalising constant
#' @param ord The order of saddlepoint approximation, 1: first-order, 2:
#' second-order first version, 3: second-order second version. Default is 3 for
#' better accuracy. See \insertCite{KumeWood05}{ellipsoidgaussian}
#' @param keyword 'Richardson' or 'simple', the method used for gradient
#' approximation.
#'
#' @returns The gradient w.r.t. the first parameter in the pseudo-normalising
#' constant.
calcVecGrad_all <- function(para1, ord_eigVals, eigVecs, ord, keyword) {
  # k by minibatch size
vec_grad <- apply(para1, 2,
                  function(y) numDeriv::grad(function(par)
                    {calclogPseudoconst_VecParGrad_fast(par,ord_eigVals, eigVecs,ord)},
                                             y, method = keyword))
return(vec_grad);

}




#' Gradient change due to log transformation
#'
#' @description
#' `logTransGradAdjust` computes the derivative w.r.t log(\code{var}) given
#'  the derivative w.r.t \code{var}.
#'
#' @param var Variable before taking the logarithm.
#' @param orig_deri The derivative w.r.t the variable.
#' @param lowerB The lower bound of the variable \code{var}.
#'
#' @returns Derivative w.r.t log(\code{var}).
logTransGradAdjust <- function(var, orig_deri, lowerB = 0) {
  # var: var before logged (has to be positive)
  res <- (var - lowerB) * orig_deri
  return(res)
}

#' Calculate derivative for tau
#'
#' @description
#' `dtau_dlogit` computes the derivative of tau w.r.t
#' logit((tau - lowerB) / (upperB - lowerB)).
#'
#' @param ltau The parameter that is being updated in gSGNHT, namely
#' logit((tau - lowerB) / (upperB - lowerB)).
#' @param lowerB Lower bound of tau
#' @param upperB Upper bound of tau
#'
dtau_dlogit <- function(ltau, lowerB, upperB) {
  # p = (tau - lowerB) / (UpperB - lowerB)
  # logit(p) = target var
  res <- exp(ltau) * (upperB - lowerB) / (1 + exp(ltau))^2
  #  res <- (tau - lowerB) * (upperB - tau) / (upperB - lowerB)
  return(res)
}


#' Calculate derivative
#'
#' @description
#' `dlogDtauDlogit_dlogit` calculates log(DtauDlogit) w.r.t the
#' logit function.
#'
#' @param logit_val The transformed tau gSGNHT is updating.
#'
dlogDtauDlogit_dlogit <- function(logit_val) {
  #logit_val: the transformed tau we are updating
  #  logit_val - 2*log(1 + exp(logit_val))
  exp_ <- exp(logit_val)
  res <- (1 - exp_) /(1 + exp_)
  return(res)
}


#' Truncated normal gradient
#'
#' @description
#' `calcHalfNormalPriorGrad` calculates the gradient w.r.t. \eqn{\log(1 / \sigma^2)},
#' with \eqn{\sigma^2} following a prior distribution N(\code{mean}, \code{asig}\eqn{^2}).
#'
#' @param linvSigma A vector of length p, \eqn{\log(1 / \sigma_j^2), j = 1, \ldots,p}.
#' @param logInv Logical, if logInv, compute the gradient w.r.t. \eqn{\log(1 / \text{var})};
#' else, compute the gradient w.r.t. \eqn{\log(\text{var})}, where \eqn{\text{var}}
#' follows a truncated normal prior in both cases.
#' @param asig The standard deviation parameter in the normal distribution.
#' @param mean The mean parameter in the normal distribution.
#'
#' @returns A vector of gradient
calcHalfNormalPriorGrad <- function(linvSigma, logInv, asig, mean = 0) {
  #calcHalfNormalPriorGrad is gradien w.r.t linvSigma, where sigma^2 \sim halfnormal(mean, asig^2)
  # warning: have not verified
  if (logInv) {
    sign <- -1
  }
  else {
    sign <- 1
  }
  if (mean == 0) {
    res <- exp(sign * 2 * linvSigma) / asig^2 - 1
  }
  else {
    res <- exp(sign * linvSigma) * (exp(sign * linvSigma) - mean) / asig^2 - 1
  }
  res <- -res * sign

  return(res)
}


#' Derivative w.r.t. the function of interest of tau
#'
#' @description
#' `calcgradTau_truncatedNormal_prior` calculates derivative w.r.t the function of tau,
#' that is actually being updated in the sampler, with tau following a truncated
#' normal prior \eqn{\text{N}(\text{mean}, \text{sd}^2)\mathbm{1}(l,u).}
#'
#' @details
#' The parameter that is being updated in the sampler is
#' \eqn{\text{logit}\{(\tau - l)/(u - l)\}.}
#'
#' @param ltau The parameter that is being updated in the sampling procedure.
#' @param prior_par A list of prior parameters, ordered by names \code{lowerB},
#' \code{upperB}, \code{cvmf} (mean), \code{dvmf} (sd).
#'
calcgradTau_truncatedNormal_prior <- function(ltau, prior_par) {
  # dvmf: standard deviation
  # cvmf: mean of the normal
  stopifnot(prior_par$fac == 1)
  tau <- transformedTau2Tau(ltau, prior_par$lowerB, prior_par$upperB, fac = 1)
  deri_normal_prior <- (prior_par$cvmf - tau) / prior_par$dvmf^2
  ans <- deri_normal_prior * dtau_dlogit(ltau, prior_par$lowerB, prior_par$upperB) +
    dlogDtauDlogit_dlogit(ltau)
  return(ans)
}

#' Derivative of Laplace prior
#'
#' @description
#' `calcloglaplaceDeri` calculates the derivative w.r.t. \code{var}.
#'
#' @param var The variable.
#' @param scale Scale of the distribution.
#' @param mean Mean of the distribution.
#'
calcloglaplaceDeri <- function(var, scale, mean = 0) {
  res <- 1 / scale
  if (var > mean) {
    res <- -res
  }
  return(res)
}



#' Approximate derivative w.r.t.
#' \eqn{\text{logit}\{(\tau - \text{lowerB}) / (\text{upperB} - \text{lowerB})\}.}
#'
#' @description
#' `approxDeriTransformedTau` approximates the derivative of [calclogC_up2const()]
#' w.r.t.  \eqn{\text{logit}\{(\tau - \text{lowerB}) / (\text{upperB} - \text{lowerB})\}.}
#'
#' @param transformedPara \eqn{\text{logit}\{(\tau - \text{lowerB}) / (\text{upperB} - \text{lowerB})\}}.
#' @param lowerB The lower bound of tau.
#' @param upperB The upper bound of tau.
#' @param minibatchSize The minibatch size.
#' @param k The latent dimension.
#' @param fac Default is 1.
#'
approxDeriTransformedTau <- function(transformedPara, lowerB, upperB, minibatchSize, k, fac) {
  # tau has truncated gamma distribution
  # p = (tau - lowerB) / (upperB - lowerB);
  # logit(p) = p / (1 - p) = temp
  # target var = temp * fac
  res <- numDeriv::grad(function(par) {
    x <- transformedTau2Tau(par, lowerB, upperB, fac)
    calclogC_up2const(x, k, logged = F)},
    transformedPara, method = 'Richardson')
  res <- res * minibatchSize
  return(res)

}
#' Approximate derivative w.r.t. \eqn{\log(\tau)}
#'
#' @description
#' `approxDeriltau` calculates derivative w.r.t. \eqn{\log(\tau)}.
#'
#' @param ltau \eqn{\log(\tau)}.
#' @param minibatchSize Minibatch size.
#' @param k The latent dimension.
#'
approxDeriltau <- function(ltau, minibatchSize,k) {
  # the current version does not allow for using a lowerB
  # check old code for the older version
  res <- numDeriv::grad(function(par) {calclogC_up2const(par, k, logged = T)},
                        ltau, method = 'Richardson')
  #
  res <- res * minibatchSize
  return(res)
}

#' Calculate derivative w.r.t \eqn{\log(1 / \sigma^2)}
#'
#' @description
#' `calcDerilinvSig` calculates the derivative w.r.t. \eqn{\log(1 / \sigma^2)},
#' where \eqn{\sigma^2} follows a normal prior.
#'
#' @param invSig1 The precision.
#' @param centered_sdat Centered minibatch, a numeric matrix.
#'
#' @returns A vector, gradient w.r.t \eqn{\log(1 / \sigma^2)}.
calcDerilinvSig <- function(invSig1, centered_sdat) { # correct!
  minibatchSize <- nrow(centered_sdat)
  term1 <- 0.5 / invSig1 * minibatchSize
  #centered_sdat <- sweep(sdat, 2, center1, FUN = '-')
  term2 <- -colSums(centered_sdat^2) / 2
  res <- term1 + term2
  res <- logTransGradAdjust(invSig1, res)
  # cat('part 1', res, '\n')
  return(res)
}

#' Calculate derivative for parameters with a normal prior
#'
#' @description
#' `calcDericenter` calculates derivative w.r.t parameters (e.g. center and )
#'
#' @param invSig1 The inverse of noise variances, precision.
#' @param centered_sdat Centered minibatch, a numeric matrix.
calcDericenter <- function(invSig1, centered_sdat) { # correct!
  # (x - c)^T diag(invSig)
  temp <- sweep(centered_sdat, 2, invSig1, FUN = '*') # n by p
  res <- colSums(temp)
  return(res)
}


