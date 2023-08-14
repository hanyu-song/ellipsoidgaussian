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
    dat <- dat[idx,, drop = F]
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




#calcGradCov2 <- function(grad_2mom, avgGrad, nentry) {
  # grad_2mom: a list
  # avgGrad: a list of the same shape as grad_2mom
  # nentry: minibatchsize
 # mapply(function(X, Y) {X / (nentry - 1) - nentry / (nentry - 1) * Y^2},
   #      X = grad_2mom, Y = avgGrad)
#}



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

#' Prior gradient
#'
#' @description
#' `calcpriorgrad` computes the prior gradient w.r.t the parameters in \code{theta1}.
#'
#' @param theta1 A list of the parameter values, ordere by their names
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
                                          logInv = T,
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
                                          logInv = F,
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
        gn2[[i]] <- approxDeriTransformedTau(theta1$ltau, tauBounds$lowerB, tauBounds$upperB,n,k,
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

