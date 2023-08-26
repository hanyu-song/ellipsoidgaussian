#' Computes stochastic gradient
#'
#' @description
#' `calcSGgrad` computes the stochastic gradient for \code{dat} at parameter
#'  values \code{theta1}.
#'
#' @param theta1 A list of parameter values at which gradient is evaluated.
#' @param dat The minibatch, a numeric matrix of n by p
#' @param k The latent dimension.
#' @param prior_para A list of prior parameters.
#' @param accurate_grad Logical, whether to evaluate the accurate gradient.
#' @param minibatchSize Minibatch size.
#' @param initialization_period Logical, whether it is the initialization period.
#' @param prior_exclude Logical, whether to ignore the priors in evaluating
#' the gradient.
#'
#' @importFrom numDeriv grad
calcSGgrad <- function(theta1,
                       dat,
                       k,
                       prior_para,
                       accurate_grad,
                       minibatchSize,
                       initialization_period,
                       prior_exclude) {
# add two input par - initialization_period - percentage
  # initialization_period: index of points to be replicated
  n <- nrow(dat);p <- ncol(dat)
 # stopifnot(minibatchSize <= n)
  if(minibatchSize < n) {
  ##  intersection <- intersect(initialization_period, idx)
  #   if (print_status) {
  #   df <- data.frame(dat); colnames(df) <- paste0('X',1:ncol(df)); df$ind <- 1; df$ind[intersection] <- 0
  #   df$ind <- factor(df$ind)
  #   q <- ggplot(df, aes(x = X1, y = X2, col = ind, group = ind)) + geom_point()
  #   print(q)
  # #  q <- ggplot(df, aes(x = X1, y = X3, col = ind, group = ind)) + geom_point()
  # #  print(q)
  # }
    if(is.null(initialization_period) ) {
      idx <- sample(1:n, minibatchSize)
    }
    else {
      # oversample the small density data points
      prob <- base::rep(1, n)
      prob[initialization_period] <- 10
      idx <- sample(1:n, minibatchSize, prob = prob)
    #  print_status <- F
     #    if (print_status) {
      #   df <- data.frame(dat); colnames(df) <- paste0('X',1:ncol(df)); df$ind <- 1; df$ind[idx] <- 0
      #   df$ind <- factor(df$ind)
      #   q <- ggplot(df, aes(x = X1, y = X2, col = ind, group = ind)) + geom_point()
      #   print(q)
      # #  q <- ggplot(df, aes(x = X1, y = X3, col = ind, group = ind)) + geom_point()
      # #  print(q)
    #   }
    }
    dat <- dat[idx,, drop = FALSE]
  }
 # if (initialization_period) {
 #   # evaluate likelihood
 #   den <- calcMargDatallh_onebyone_cpp(exp(-theta1$linvSig), dat, theta1$lambda,
 #                                transformedTau2Tau(theta1$ltau,
 #                                                   prior_para$ltau$lowerB, prior_para$ltau$upperB,
 #                                                   prior_para$ltau$fac),
 #                                theta1$mu, theta1$center)
 #   idx <- order(den, decreasing = F)[1:ceiling(0.1 * minibatchSize)]
 #   ######### print - visualization of the subset
 # # if (print_status) {
 # # df <- data.frame(dat); colnames(df) <- paste0('X',1:ncol(df)); df$ind <- 1; df$ind[idx] <- 0
 # # df$ind <- factor(df$ind)
 # # q <- ggplot(df, aes(x = X1, y = X2, col = ind, group = ind)) + geom_point()
 # # print(q)
 # # q <- ggplot(df, aes(x = X1, y = X3, col = ind, group = ind)) + geom_point()
 # # print(q)
 # # }
 # #   ###########
 #   # 10% upweighted by 4 + 1 times
 #   dat <- rbind(dat, apply(dat[idx,,drop = F], 2, rep, 3))
 # }
 # else {
 # sdat <- dat
 # }
  ltau_ <- ifelse(prior_para$ltau$update,'non','ltau')
  par_names <- setdiff(names(prior_para), ltau_)
#  par_names <- names(prior_para)
  fac <- -n / minibatchSize
  grad <-  calcllhgrad(par_names,theta1, dat, k,
                       prior_para$ltau, # needs modification
                       accurate_grad)
  #dUtilde <- -n / minibatchSize * calcllhgrad(par_names,theta1, dat, k,
      #                                        prior_para$ltau[c('lowerB','upperB')],
                                      #        accurate_grad)
  dUtilde <- lapply(grad, '*', fac)
  if (!prior_exclude) {
 # dUtilde <- dUtilde - calcpriorgrad(theta1, p, k,prior_para)
    prior_grad <- calcpriorgrad(theta1, p, k,prior_para[par_names]) ## (possibly change) this line
    temp <- intersect(names(dUtilde),names(prior_grad))
    dUtilde <- Map('-',dUtilde[temp], prior_grad[temp] )
   # dUtilde <- dUtilde - calcpriorgrad(theta1, p, k,prior_para)
  }

  return(dUtilde)
}



#' Prior gradient
#'
#' @description
#' `calcpriorgrad` computes the prior gradient w.r.t the parameters in \code{theta1}.
#'
#' @param theta1 A list of the parameter values, ordered by their names
#' @param p The ambient dimension
#' @param k The latent dimension
#' @param prior_para A list of the prior parameters
#'
calcpriorgrad <- function(theta1, p, k, prior_para) {
  # para_idx: a list, each element has the name of the parameter,
  # and the content is the para index
  res <- vector('list', length = length(prior_para))
  names(res) <- names(prior_para)
  for (i in seq_along(res)) {
    par_name <- names(res)[i]
    if (par_name %in% 'ltau') {
      # too clucky, need to change!
      #  CHANGE!!!
      res[[i]] <- calcgradTau_truncatedNormal_prior(theta1$ltau, prior_para$ltau)

    #  res[[i]] <- calcGrad_transformedTau_prior(theta1$ltau, prior_para$ltau)
      }
    else if (par_name %in% 'lambda') {
      if (length(prior_para$lambda) == 0) {
        res[[i]] <- matrix(0, p, k)
      }
      else {
   #   lambda <- matrix(theta1[para_idx[[par_name]]], nrow = p, ncol = k, byrow = F)
    #  res[[i]] <- calcgradLambdaprior(Lambda, prior_para$lambda$phi, prior_para$lambda$kappa)
       # print('here')
      res[[i]] <- calcgradLambdaprior_cpp(theta1$lambda, prior_para$lambda$phi, prior_para$lambda$kappa)
      res[[i]] <- matrix(res[[i]], nrow = p, ncol = k)
     # print('prior')

      }
      }
    else if (par_name %in% 'linvSig') {
     # linvSigma <- theta1[para_idx[[par_name]]]
      res[[i]] <- calcHalfNormalPriorGrad(theta1$linvSig,
                                          logInv = TRUE,
                                          prior_para$linvSig$asig)  # not verified yet
      }
    else if (par_name %in% 'center') {
    #  center1 <- theta1[para_idx[[par_name]]]
      res[[i]] <- -(theta1$center - prior_para$center$cenmu) / prior_para$center$censig^2
      }
    else if (par_name %in% 'mu') {
        # this is in fact reporting the grad. of f(\Lambda\mu)
        # w.r.t \Lambda\mu, meaning \Lambda\mu is updated in this algorihtm
        # lambdamu: p by 1
        if (!is.null(prior_para$mu$type) && (prior_para$mu$type %in% 'unconstrained')) {
          # mu = y / ||y||, derivative here is w.r.t y, y ~ N(0,1)
          #res[[i]] <- -theta1[para_idx[[par_name]]]
          res[[i]] <- -theta1$mu
        }
        else {
          # mu ~ vmf(tau = 0), derivative here is w.r.t. mu
          # res[[i]] <- rep(0, k)
          # mu ~ vmf(tau = tau0, mu = mu0)
          res[[i]] <- prior_para$mu$tau0 * prior_para$mu$mu0
        }
        }
    else if (par_name %in% 'axesLen') {
      # not checked
      #temp1 <- calcHalfNormalPriorGrad(-theta1$axesLen[1], logInv = F,
                      #        asig = prior_para$axesLen$sig1,
                       #       mean = prior_para$axesLen$mu1)
      #temp2 <-calcHalfNormalPriorGrad(-theta1$axesLen[-1], logInv = F,
                        #      asig = prior_para$axesLen$sig_incre,
                         #     mean = prior_para$axesLen$mu_incre)
      #res[[i]] <- c(temp1, temp2)
      #### new part - reparametrize (N(0, \sigma^2))
     #res[[i]] <- -(theta1$axesLen - prior_para$axesLen$mu1) / prior_para$axesLen$sig1^2
      ### new part 2 - reparmetrize (N(0, \sigma^2)1(0, \infty))
      res[[i]] <- calcHalfNormalPriorGrad(theta1$axesLen,
                                          logInv = FALSE,
                                          prior_para$axesLen$sig1,
                                          mean = prior_para$axesLen$mu1)  # not verified yet

    }
    else if (par_name %in% 'axesDir') {
      ########### original
      res[[i]] <- matrix(0, p, k)
      #############
      ########### experiment - Bayesian constraint relaxation
      #res[[i]] <- calcgradLambdaprior_cpp(theta1$axesDir, prior_para$axesDir$phi, prior_para$axesDir$kappa)
      #temp <- constraint_fac_grad(theta1$axesDir, p, k, 2)
      #res[[i]] <- matrix(res[[i]] + temp, nrow = p, ncol = k)
      #############################
    }
    else {
        stop("The parameter name is not expected.")
      }
    }
 # res <- unlist(res)
  return(res)
}



#' Gradient of the log-likelihood
#'
#' @description
#' `calcllhgrad` calculates gradient of log-likelihood w.r.t. the parameters in
#' \code{par_names}.
#'
#' @param par_names A list of parameters to be updated
#' @param theta1 A list of parameter values
#' @param sdat The minibatch data.
#' @param k The latent dimension.
#' @param tauBounds A list of tau information, lowerB, upperB and fac.
#' @param accurate_grad Logical, whether accurate or simple gradient is desired.
#'
calcllhgrad <- function(par_names, theta1, sdat, k, tauBounds, accurate_grad) {
  # par_names: a list of parameters to be updated
  # center_prior_var: vector of variances
  n <- nrow(sdat)
  p <- ncol(sdat)
  stopifnot(!is.null(p))
 # par_names <- unique(c(par_names,'lambda'))
  gn2 <- vector('list', length = length(par_names))
  names(gn2) <- par_names
  for (i in seq_along(gn2)) {
    par_name <- names(gn2)[i]
    if (par_name %in% 'ltau') {
    #  ltau <- theta1[para_idx[[par_name]]]
     # print(transformedTau2Tau(ltau, tauBounds$lowerB, tauBounds$upperB))
      #  (CHANGE!!!)
      if (is.null(tauBounds$upperB)) {
        gn2[[i]] <-  approxDeriltau(theta1$ltau, n, k)
      }
      else {
        gn2[[i]] <- approxDeriTransformedTau(theta1$ltau, tauBounds$lowerB,
                                             tauBounds$upperB,n,k,
                                           tauBounds$fac)
      }
    }
    else if (par_name %in% c('lambda')) {
    #  gn2[[i]] <- rep(0, length(para_idx[[par_name]]))
      gn2[[i]] <- array(0, dim = dim(theta1$lambda))
    }
    else if (par_name %in% c('mu')) {
      gn2[[i]] <- base::rep(0, length = length(theta1$mu))
    }
    else if (par_name %in% 'linvSig') {
      #center1<- theta1[para_idx[['center']]] # CHANGED, did not check
      centered_sdat <- sweep(sdat, 2, theta1$center, FUN = '-')
     # linvSigma <- theta1[para_idx[[par_name]]]
      gn2[[i]] <- calcDerilinvSig(exp(theta1$linvSig), centered_sdat)
    }
    else if (par_name %in% 'center') {
      centered_sdat <- sweep(sdat, 2, theta1$center, FUN = '-')
      gn2[[i]] <- calcDericenter(exp(theta1$linvSig) , centered_sdat)  # theta1[1:p], parameter for linvSigma
    }
    else if (par_name %in% 'axesLen') {
      gn2[[i]] <- base::rep(0, length = length(theta1$axesLen))
    }
    else if (par_name %in% c('axesDir')) {
      #  gn2[[i]] <- rep(0, length(para_idx[[par_name]]))
      gn2[[i]] <- array(0, dim = dim(theta1$axesDir))
    }
    else {
      stop("The parameter name is not expected.")
    }
  }

 # para_idx_vec <- unlist(para_idx)
  #if (accurate_grad) {
  # gn <- numDeriv::grad(function(par) {
   #   theta1[para_idx_vec] <- par
  #    calclogPseudoconst(theta1, sdat, k, sampler_type)},
  #                     theta1[para_idx_vec], method = 'Richardson')
 #   }
 # else {
   # gn <- numDeriv::grad(function(par) {
 #  theta1[para_idx_vec] <- par
  #  calclogPseudoconst(theta1, sdat, k, sampler_type)},
  #  theta1[para_idx_vec], method = 'simple')
  #}
  ### RCPP version
  gn <- calcGrad_llh_fast_wrapper(par_names, theta1, sdat, tauBounds,accurate_grad)
  # needs modification
  temp <- intersect(names(gn),names(gn2))
  ###############3
  res <- Map("+", gn[temp], gn2[temp])
  return(res)
}


# calcGradlogPseudoconst <- function(theta1, sdat, k) {
#   n <- nrow(sdat)
#   p <- ncol(sdat)
#   len <- length(theta1)
#   bd_idx <- p + p * k
#   linvSigma <- theta1[1:p]
#   idx <- (p + 1) : (bd_idx)
#   lambda <- matrix(theta1[idx], nrow = p, ncol = k, byrow = F)
#   #  ltau <- theta1[bd_idx + 1]
#   ltau <- theta1[bd_idx + 1]
#   idx2 <- (bd_idx + 2):(len - p)
#   mu <- theta1[idx2]
#   center1 <- theta1[(len - p + 1):len]
#   centered_dat <- sweep(sdat, 2, center1, FUN = '-')
#   # res <- dEllipsoidGauss_psedo_const(invSigma, lambda, exp(ltau), mu, centered_dat)
#   res <- dEllipsoidGauss_psedo_const(exp(linvSigma), lambda, exp(ltau), mu, centered_dat)
#   return(res)
# }



#lambda <- matrix(theta1[par_idx$lambda], nrow = p, ncol = k, byrow = F)
#tau <- transformedTau2Tau(theta1[par_idx$ltau], tauBounds$lowerB, tauBounds$upperB)
#mu <- theta1[par_idx$mu] # unconstrained mu
#invSigma <- exp(theta1[par_idx$linvSig])



#' Wrapper of [calcGrad_llh_fast()]
#'
#' @description
#' `calcGrad_llh_fast_wrapper` is wrapper to [calcGrad_llh_fast()].
#'
#' @param par_names Parameter names.
#' @param theta1 The parameter values.
#' @param sdat The minibatch.
#' @param tauBounds A list of bounds for tau, lowerB and upperB.
#' @param accurate_grad Logical, whether accurate gradient is desired.
#'
calcGrad_llh_fast_wrapper <- function(par_names,theta1, sdat, tauBounds, accurate_grad) {
  #paras <- separateParas(para_idx, theta1, sdat, tauBounds)
#  res <- calcGrad_llh_fast(par_names, exp(paras$linvSig, paras$lambda, paras$tau, paras$mu, paras$centered_dat, tauBounds, accurate_grad)
#  if(!(is.null(tauBounds$upperB))) {
    tau <- transformedTau2Tau(theta1$ltau, tauBounds$lowerB, tauBounds$upperB, tauBounds$fac)
 # }
 # else {
   # tau <- exp(theta1$ltau) + tauBounds$lowerB
  #}
  centered_dat <- sweep(sdat, 2, theta1$center, FUN = '-')
  if (any(c('axesLen', 'axesDir') %in% par_names)) {
    par_names <- unique(c(par_names, 'lambda'))
  }
  par_names2 <- par_names[!par_names %in% c('axesLen','axesDir')]

  res <- calcGrad_llh_fast(par_names2, exp(theta1$linvSig), theta1$lambda,
                           tau, c(theta1$mu), centered_dat,
                           tauBounds, accurate_grad)
  if (any(c('axesDir', 'axesLen') %in% par_names)) {
  increment <- exp(theta1$axesLen)
 # len <- cumsum(increment)

  if ('axesDir'%in% par_names) {
    # axesLen on the log. scale
    res$axesDir <- sweep(res$lambda, MARGIN = 2, FUN = '*',STATS = increment)
    # axeLen oriignal scale
    #res$axesDir <- sweep(res$lambda, MARGIN = 2, FUN = '*',STATS = theta1$axesLen)
  }
  if ('axesLen' %in% par_names) {
    ##### old parametrization
   # temp <- colSums(res$lambda * theta1$axesDir)
   # revcumSums <- spatstat.utils::revcumsum(temp)
   # res$axesLen <- revcumSums * increment
    #### new parametrization - no order
    res$axesLen <- colSums(res$lambda * theta1$axesDir)
    res$axesLen <- res$axesLen * increment
    #* theta1$axesLen
  }
  }


  return(res)
}


#' Calculates the gradient of the log-likelihood w.r.t each parameter
#'
#' @description
#' `calcGrad_llh_fast`is a high-level function that calculates the gradient of the log-
#' likelihood w.r.t parameters that are in \code{par_names}.
#'
#' @param par_names A list of parameter names.
#' @param invSigma A vector of the noise precision.
#' @param lambda The factor loading matrix, p by k.
#' @param tau The concentration parameter
#' @param mu The direction parameter, a vector of length k.
#' @param centered_sdat A matrix of centered data that are in the minibatch.
#' @param tauBounds The lower and upper bounds of tau
#' @param accurate_grad Logical, whether accurate gradient is desired or not.
#'
#' @returns A list of gradients for the parameters in \code{par_names}.
calcGrad_llh_fast <- function(par_names,
                              invSigma,
                              lambda,
                              tau,
                              mu,
                              centered_sdat,
                              tauBounds,
                              accurate_grad) {
  # para_interest: parameter of interests, a vector of names; names should be consistent with
  # stochatsic gradient MCMC one at a time (currently not set; assuming all parameters are of interest)
  p <- ncol(centered_sdat); k <- length(mu);   nsdat <- nrow(centered_sdat)
  # temp <- t(lambda) %*% diag(invSigma)
  lambinvSig <- sweep(t(lambda),MARGIN = 2, STATS = invSigma, FUN = '*') # equivalent to  t(lambda) %*% diag(invSigma)
  # para12 <- apply(centered_sdat, 1, function(x) t(lambda) %*% diag(invSigma) %*% x)
  para12 <- lambinvSig %*% t(centered_sdat)
  para1 <- sweep(para12, 1, tau * mu, FUN = '+')
  #para2 <- t(lambda) %*% diag(invSigma) %*% lambda / 2
  para2 <- lambinvSig %*% lambda / 2
  lowerTri <- t(chol(para2))
  lowerTriVech <- lowerTri[lower.tri(lowerTri, diag = TRUE)]
  para2_vech <- c(matrixcalc::vech(para2))
  eig_ <- eigen(para2, symmetric = TRUE)
  eigVals <- eig_$values
  eigVecs <- eig_$vectors
  ord <- order(eigVals)
  ord_eigVals <- eigVals[ord]
  if(accurate_grad) {
    keyword <- 'Richardson'
  }
  else {
    keyword <- 'simple'
  }
  # mat_grad <- numDeriv::grad(function(par) {calclogPseudoconst_MatParGrad2(par, para1)},
  #                   lowerTriVech, method = 'Richardson') # length k(k+1)/2 vec
  # mat_grad2 <- numDeriv::jacobian(function(par) {A <- rockchalk::vech2mat(par, lowerOnly = F);
  # c(vech(t(chol(A))))},
  # c(vech(para2)), method = 'Richardson')
  # psdMat_grad  <- numDeriv::grad(function(par) {calclogPseudoconst_MatParGrad4(par, para1)}, para2_vech, method = keyword)
  psdMat_grad  <- numDeriv::grad(function(par)
    {calclogPseudoconst_MatParGrad4_cpp(par, para1)},
    para2_vech, method = keyword)
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
      res[[i]] <- calc_linvSig_grad_cpp(invSigma, psdMat_grad,
                                        vec_grad, lambda, centered_sdat)
    }
    else if(par_names[i] %in% 'lambda') {
      temp <- calc_lambda_grad_cpp(lambda,psdMat_grad, vec_grad,invSigma,
                                   centered_sdat)
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
