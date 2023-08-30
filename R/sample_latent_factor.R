#' Fisher-Bingham parameters
#'
#' @description
#' `calcFB_pars` calculates the parameters for the Fisher-Bingham distribution
#' that the latent factors follow a posteriori.
#'
#' @param invSigma The precision.
#' @param lambda A matrix, p by k, factor loadings.
#' @param tau The concentration parameter.
#' @param mu A vector of length k.
#' @param center A vector of length p.
#' @param dat A numeric matrix.
#' @param norm_ Logical, whether normalized \code{para1} or not.
#'
calcFB_pars <- function(invSigma, lambda, tau, mu, center,dat, norm_) {
  # -log(normalising)+ mu_kap_prod * x - x^T A x
  centered_dat <- sweep(dat, MARGIN = 2, STATS = center, FUN = '-')
  lambinvSig <- sweep(t(lambda),MARGIN = 2, STATS = invSigma, FUN = '*') # equivalent to  t(lambda) %*% diag(invSigma)
  para12 <- lambinvSig %*% t(centered_dat)
  para1 <- sweep(para12, 1, tau * mu, FUN = '+')
  para2 <- lambinvSig %*% lambda / 2
  if(!norm_) {
    return(list(para1 = para1, para2 = para2))
  }

  #  kap <- apply(para1,2,function(y) pracma::Norm(y,2))
  kap <- apply(para1,2,function(y) sqrt(sum(abs(y^2))))
  para1 <- sweep(para1, 2, kap, FUN = '/')
  res <- list(para1 = para1, para2 = para2, kap = kap)
  return(res)
}




#' Simulation from Fisher-Bingham distribution
#'
#' @description
#' `rFB_MH` simulates from a Fisher-Bingham distribution using
#' Hoff (2009).
#'
#' @param eta1 A vector of length k
#' @param para1 A vector of length k
#' @param A A square matrix, k by k.
#'
#' @importFrom Rfast rvmf
#' @importFrom stats runif
#'
#' @references
#' Peter D. Hoff (2009) Simulation of the Matrix Bingham–von Mises–Fisher
#' Distribution, With Applications to Multivariate and Relational Data,
#' Journal of Computational and Graphical Statistics, 18:2, 438-456,
#' DOI: 10.1198/jcgs.2009.07177
rFB_MH <- function(eta1, para1, A) {
  # tau1: L2 norm of para1; if null, have to compute it
  # para1: a vector
  # para2: k by k
  if (is.null(eta1)) {
    eta1 <- para1
  }
  eta1 <- c(eta1)
  # k <- length(eta1)
  cur <- eta1 / sum(eta1^2)
  # proposal:
  #  prop <- one_rvmf(eta1,40)  # can change to 40
  prop <- Rfast::rvmf(1, eta1, 40)
  #  Omega <- diag(k) + 2 * A + tau1 * (diag(k) - tcrossprod(para1))
  # Sigma <- chol2inv(chol(Omega))
  #chol_Sigma <- chol(Sigma)
  # prop <- my_racg(1, chol_Sigma) # draw n of them at the same time
  logAcceptRatio <- dFB_cpp_log(c(prop), para1, A, 0) -
    dFB_cpp_log(cur, para1, A, 0)
  # dacg_up2const(eta1, Omega, log_ = T) - dacg_up2const(prop, Omega, log_ = T)
  u <- log(stats::runif(1))
  acpt_status <- (u<= logAcceptRatio)
  if (acpt_status) {
    cur <- prop
  }

  return(list(acpt_status = acpt_status, res = cur)) # n by k
}



#' Impute the missing values
#'
#' @description
#' `predictY` imputes the missing entries in \code{dat} by assuming that
#' the data comes from an ellipsoid-Gaussian distribution.
#'
#' @param temp A list of the distribution parameters.
#' @param tau_par A list of parameters regarding tau.
#' @param dat A numeric matrix, data.
#' @param missing_idx Indices of the missing entries.
#' @param mis_pos Positions of the missing entries.
#'
#' @importFrom condMVNorm condMVN
#' @importFrom Rfast rmvnorm
predictY <- function(temp, tau_par, dat,missing_idx,mis_pos) {
  # has issue with tau
  # dat: does not contain missingness
  # missing_idx: missing idx of the original missing data subset
  invSig <- exp(temp$linvSig)
  tau <- transformedTau2Tau(temp$ltau,
                            tau_par$lowerB, tau_par$upperB,
                            tau_par$fac)
  # lat_fac <- drawLatFac(invSig, temp$lambda, tau, temp$mu, temp$center, dat)
  lat_fac <- drawLatFac_chain(invSig, temp$lambda, tau, temp$mu, temp$center, dat,30)
  # print(lat_fac$acpt_rate) #all 1, may not be good
  # if there is arbitrary missing pattern, this step can be easily adapted for that
  # assume that the last column is the missing y

  #mean <- sweep(temp$lambda[ncol(dat),,drop = F] %*% t(lat_fac$latFacs), MARGIN = 1,
  #  STATS = temp$center, FUN = '+') # p by n
  # problem? should be conditonal distribution rather than marginal?
  ## may be wrong
  ## mean <- apply(lat_fac$latFacs, 1, function(x) sum(x * temp$lambda[ncol(dat),])) +
  ##  temp$center[ncol(dat)]
  ## sd <- 1 / sqrt(invSig[ncol(dat)])
  ##
  sigma <- diag(1 / sqrt(invSig))
  for (i in seq_along(missing_idx)) {
    #    row <- as.numeric(names(missing_idx)[i])
    dep <- missing_idx[[i]]
    given <- setdiff(seq_len(ncol(dat)), dep)
    # should speed this function up
    par <- condMVNorm::condMVN(temp$lambda %*% t(lat_fac$latFacs[i,,drop = FALSE]) +
                                 temp$center, sigma = sigma, dep = dep,
                               given =  given, X.given = dat[i, given])
    # should speed this function up
    dat[i, dep] <- Rfast::rmvnorm(1, par$condMean, sigma = par$condVar)
  }
  # browser(expr = {any(dat[mis_pos] > 9)})

  # return(list(dat = dat, acpt_rate = lat_fac$acpt_rate)) # may change this line to remove acpt_rate
  return(dat[mis_pos])
}

#' Simulate the latent factors.
#'
#' @description
#' `drawLatFac_chain` simulates the latent factors in the von-Mises Fisher linear
#' factor model. They follow a Fisher-Bingham distribution a posteriori.
#'
#' @param invsig The precision.
#' @param lambda A matrix, the factor loading matrix.
#' @param tau The concentration parameter.
#' @param mu A vector.
#' @param center A vector.
#' @param dat A numeric matrix.
#' @param nsample Number of samples to draw from the Metropolis-Hasting procedure.
#'
drawLatFac_chain <- function(invsig, lambda, tau, mu, center, dat, nsample = 20) {
  latFacs <- matrix(NA, nrow = nrow(dat), ncol = length(mu))
  acpt <- 0
  pars <- calcFB_pars(invsig,
                      lambda,
                      tau,mu, center,dat,norm_ = FALSE)
  for (l in seq_len(nrow(dat))) {
    # print('new')
    # print(pars$para1[,l])
    temp <- rFB_MH(NULL, pars$para1[,l], pars$para2)
    for (ss in seq_len(nsample)) {
      #   print(temp$res)
      temp <- rFB_MH(temp$res, pars$para1[,l], pars$para2)
      #print(apply(pars$para1,2,function(y) pracma::Norm(y,2)))
      acpt <- acpt + temp$acpt_status
      # temp <- rLatFac_MH(center, invsig, lambda, pars$para1[,l] / pars$kap[l],
      #  mu, tau, dat[l,,drop = F])
      #   print(temp$res)
    }
    latFacs[l,] <- temp$res
    #  print(temp$acpt_status)
    #  print(anyNA(temp$res))

    # ct <- ct + 1
    #  cat('i = ',i, 'j =',j, 'l =',l,'\n')
  }
  acpt_rate <- acpt / nrow(dat) / nsample
  # print(acpt_rate)
  #   latFacs[i,,j] <- updateLatFac_MH(latFacs[i,,j - 1], res_gibbs$lambda[ii,,], res_gibbs$invsig[ii,],
  #   res_gibbs$tau[ii,1], res_gibbs$mu[ii,], res_gibbs$center[ii,], dat)
  return(list(latFacs = latFacs, acpt_rate = acpt_rate))
}


#' Log density of the Fisher-Bingham distribution
#'
#' @description
#' `dFB_cpp_log` calculates log density of the Fisher-Bingham distribution.
#'
#' @param X A matrix, 1 by n.
#' @param para1 A matrix, 1 by n.
#' @param para2 The matrix parameter, n by n.
#' @param fb_const_part The Fisher-Bingham constant. If NULL, the
#' function will compute it.
dFB_cpp_log <- function(X,para1, para2, fb_const_part) {
  # X: 1 by n
  # para1: 1 by n
  if (is.null(fb_const_part)) {
    fb_const_part <- approxFBconst_cpp(para1, para2, idx_ = 3) # n by 1
  }
  # term2 <- apply(X, 2, function(y) t(y) %*% para2 %*% y)
  temp1 <-  para2 %*% X # k by n

  if (class(X) %in% 'matrix') {
    term2 <- colSums(X * temp1)
    res <- -fb_const_part + colSums(X * para1) - term2
  }
  else {
    term2 <- sum(X * temp1)
    res <- -fb_const_part + sum(X * para1) - term2
  }
  return(res)
}

#' Sample latent factors
#'
#' @description
#' `draw_latent_factors` simulates \eqn{\min(500,
#' \text{number of post burn-in samples})} latent factors for each observation from
#' their posterior distribution, a Fisher-Bingham distribution. Simulating from
#' a Fisher-Bingham is challenging and we use the method by Hoff (2009).
#' @param dat A matrix of data, n by p.
#' @param samples A list of posterior samples, output of [ellipsoid_gaussian()].
#' @param burnin Number of burn-in samples.
#' @param num_samp Number of samples for each latent factor.
#' @export
#' @references
#' Peter D. Hoff (2009) Simulation of the Matrix Bingham–von Mises–Fisher
#' Distribution, With Applications to Multivariate and Relational Data,
#' Journal of Computational and Graphical Statistics, 18:2, 438-456,
#' DOI: 10.1198/jcgs.2009.07177
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' lat_facs <- draw_latent_factors(shell, samples, burnin = 2500, num_samp = 5)
#'
draw_latent_factors <- function(dat, samples, burnin = NULL, num_samp = 200) {
  #  res <- matrix(NA, 1000, nrow(dat))

  if (is.null(burnin)) {
    burnin <- nrow(samples$tau) / 2
  }
  niter <- nrow(samples$tau) - burnin + 1
  num_samp <- min(num_samp, niter)
  idx <- sample((burnin + 1):nrow(samples$tau), num_samp)
  iterations <- numeric(num_samp)
  lat_facs <- array(NA, c(nrow(dat), dim(samples$lambda)[3], num_samp))
  progress_bar <- utils::txtProgressBar(min=0, max=num_samp, style = 3, char="=")

  for (i in seq_along(idx)) {
    iter <- idx[i]
    iterations[i] <- iter
    temp2<- drawLatFac_chain(samples$invsig[iter,], samples$lambda[iter,,],
                             samples$tau[iter,1], samples$mu[iter,],
                             samples$center[iter,], dat,30)
    lat_facs[,,i] <- temp2$latFacs
    #   temp <- predictY2(samples$invsig[iter,], samples$lambda[iter,,],
    #                samples$center[iter,], lat_facs[,,i],dat)
    #print(temp$acpt_rate)
    #  res[i,] <- temp$y
    #  cat('iteration',i,'completed.\n')
    utils::setTxtProgressBar(progress_bar, value = i)
  }
  close(progress_bar)

  return(list(lat_facs = lat_facs, iterations = iterations))
}

#' Resolve rotational ambiguity in samples of factor loadings and factors jointly
#'
#' @description
#' `joint_rot_samples` performs the varimax rotation on the factor loadings samples
#' and column-based matching to resolve sign and label switching; see the help page of
#' [infinitefactor::jointRot()] for more details.
#'
#' @param lambda_samps The lambda samples in the list of samples as output of
#' [ellipsoid_gaussian()].
#' @param lat_facs The latent factors, an array of samples, number of observations
#' x latent dimensions x number of samples
#' @param iterations Iteration indices of the samples that are used for drawing
#' the latent factors, a vector.
#' @export
#'
joint_rot_samples <- function(lambda_samps, lat_facs, iterations) {
  if (!requireNamespace("infinitefactor", quietly = TRUE)) {
    stop(
      "Package \"infinitefactor\" must be installed to use this function.",
      call. = FALSE
    )
  }
  num_samp <- dim(lat_facs)[3]
  etaSamples <- vector('list', length = num_samp)
  lambSamples <- vector('list', length = num_samp)
  for (i in seq_len(num_samp)) {
    iter <- iterations[i]
    etaSamples[[i]] <- lat_facs[,,i]
    lambSamples[[i]] <- lambda_samps[iter,,]
  }
  aligned <- infinitefactor::jointRot(lambSamples, etaSamples)# a list of niter
  return(aligned)
}
#' Post-process the factor loadings
#'
#' @description
#' `postprocess` post process the samples of the factor loadings to resolve
#' rotational ambiguity. It performs the varimax rotation on the factor loadings samples
#' and column-based matching to resolve sign and label switching; see the help page of
#' [infinitefactor::jointRot()] for more details.
#'
#' @param dat A matrix of data, n by p.
#' @param samples A list of samples from the posterior distribution, output of
#' [ellipsoid_gaussian()].
#' @param burnin Number of burn-in samples.
#' @param num_samps Number of samples for each latent factor
#' @export
#' @examples
#' aligned <- postprocess(shell, samples, burnin = 2500, num_samps = 5)
#'
postprocess <- function(dat, samples, burnin = NULL, num_samps = 200) {
  lat_facs <- draw_latent_factors(dat, samples, burnin, num_samp = num_samps)
  res <- joint_rot_samples(samples$lambda, lat_facs$lat_facs, lat_facs$iterations)
  return(res)
}
#' Visualize the factor loadings
#'
#' @description
#' `plot_factor_loadings` visualizes the posterior mean of the factor loadings matrix. We
#' recommend post processing the samples first using [draw_latent_factors()] and
#' [postprocess()].
#'
#' @param lambda_samples A list of the samples, length = number of samples,
#' each element is p by k.
#' @param row_idx A vector row indices of the factor loadings to be plotted. Default
#' is all the rows are included.
#' @export
plot_factor_loadings <- function(lambda_samples, row_idx = NULL) {
  if (!requireNamespace("infinitefactor", quietly = TRUE)) {
    stop(
      "Package \"infinitefactor\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (is.null(row_idx)) {
    row_idx <- seq_len(nrow(lambda_samples[[1]]))
  }
  mat <- infinitefactor::lmean(lambda_samples)[row_idx,]
  p1 <- infinitefactor::plotmat(mat)
  return(p1)
}

