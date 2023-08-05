### log likelihood derivative w.r.t mu, center and Lambda - for speed up purpose

#' @importFrom numDeriv grad
#' @importFrom matrixcalc vech
calclogPseudoconst_VecParGrad_fast <- function(vec_par, eigVals, eigVecs, ord_eigVal, order_ = 3) {
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





calcGrad_llh_fast <- function(par_names,invSigma, lambda, tau, mu, centered_sdat,
                              tauBounds,accurate_grad) {
  # para_interest: parameter of interests, a vector of names; names should be consistent with
  # stochatsic gradient MCMC one at a time (currently not set; assuming all parameters are of interest)
  p <- ncol(centered_sdat); k <- length(mu);   nsdat <- nrow(centered_sdat)
# temp <- t(lambda) %*% diag(invSigma)
  lambinvSig = sweep(t(lambda),MARGIN = 2, STATS = invSigma, FUN = '*') # equivalent to  t(lambda) %*% diag(invSigma)
 # para12 <- apply(centered_sdat, 1, function(x) t(lambda) %*% diag(invSigma) %*% x)
  para12 <- lambinvSig %*% t(centered_sdat)
  para1 <- sweep(para12, 1, tau * mu, FUN = '+')
  #para2 <- t(lambda) %*% diag(invSigma) %*% lambda / 2
  para2 <- lambinvSig %*% lambda / 2
  lowerTri <- t(chol(para2)); lowerTriVech <- lowerTri[lower.tri(lowerTri, diag = T)];
  para2_vech <- c(matrixcalc::vech(para2))
  eig_ <- eigen(para2, symmetric = T)
  eigVals <- eig_$values; eigVecs <- eig_$vectors;
  ord <- order(eigVals); ord_eigVals <- eigVals[ord]
  if(accurate_grad) {
    keyword = 'Richardson'
  }
  else {
    keyword = 'simple'
  }
 # mat_grad <- numDeriv::grad(function(par) {calclogPseudoconst_MatParGrad2(par, para1)},
      #                   lowerTriVech, method = 'Richardson') # length k(k+1)/2 vec
 # mat_grad2 <- numDeriv::jacobian(function(par) {A <- rockchalk::vech2mat(par, lowerOnly = F);
 # c(vech(t(chol(A))))},
 # c(vech(para2)), method = 'Richardson')
 # psdMat_grad  <- numDeriv::grad(function(par) {calclogPseudoconst_MatParGrad4(par, para1)}, para2_vech, method = keyword)
  psdMat_grad  <- numDeriv::grad(function(par) {calclogPseudoconst_MatParGrad4_cpp(par, para1)}, para2_vech, method = keyword)
#  vec_grad <- apply(para1, 2, function(y) numDeriv::grad(function(par) {calclogPseudoconst_VecParGrad(par, para2)},
                        #      y, method = keyword)) # k by n
  vec_grad <- calcVecGrad_all(para1, ord_eigVals, eigVecs, ord, keyword)
  #vec_grad <- apply(para1, 2,
                  #  function(y) numDeriv::grad(function(par) {calclogPseudoconst_VecParGrad_fast(par,ord_eigVals, eigVecs,ord)},
                                                    #     y, method = keyword)) # k by n

 # psdMat_grad <- mat_grad %*% mat_grad2 # 1 by (k)(k+1)/2
  vec_grad_sum <- rowSums(vec_grad)
  res <- vector("list", length(par_names)); names(res) <- par_names
  for(i in seq_along(res)) {
    if(par_names[i] %in% 'mu') {
      res[[i]] <- vec_grad_sum * tau
    }
    else if(par_names[i] %in% 'ltau') {
      # CHANGE !!!
      res[[i]] <- calc_transTau_grad(tau, vec_grad_sum, mu, tauBounds)
    }
    else if(par_names[i] %in% 'center') {
      res[[i]] <-  -c(vec_grad_sum %*% lambinvSig)
    }
    else if(par_names[i] %in% 'linvSig') {
      res[[i]] <- calc_linvSig_grad_cpp(invSigma, psdMat_grad, vec_grad, lambda, centered_sdat)
    }
    else if(par_names[i] %in% 'lambda') {
      temp <- calc_lambda_grad_cpp(lambda,psdMat_grad, vec_grad,invSigma, centered_sdat)
      res[[i]] <- matrix(temp, nrow = p, ncol = k)
    }
    else {
      stop('The parameter names is not expected!')
    }
  }
 # if ('mu' %in% para_interest) {
 # }
  return(res)
}
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
calcVecGrad_all <- function(para1, ord_eigVals, eigVecs, ord, keyword) {
  # k by minibatch size
vec_grad <- apply(para1, 2,
                  function(y) numDeriv::grad(function(par) {calclogPseudoconst_VecParGrad_fast(par,ord_eigVals, eigVecs,ord)},
                                             y, method = keyword))
return(vec_grad);

}

