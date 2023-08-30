#' Sample from the posterior distribution of an ellipsoid-Gaussian likelihood
#' @description
#' `gSGNHT_EG` draws from the posterior distribution associated with an ellipsoid-Gaussian
#' likelihood using geodesic Stochastic Gradient Nose-Hoover Thermostats (gSGNHT)
#' (Liu et al. 2016)
#'
#' @param dat A numeric matrix of data (n by p).
#' @param k The latent dimension.
#' @param init A list of initial values to be used for the sampler.
#' @param prior_par A list of prior parameters
#' @param alg_par A list of algorithm-related parameters.
#' @param niter The total number of iterations.
#' @param updateCenter Logical, whether to update the center or not during sampling.
#' @returns A list of samples, including samples for \code{axesDir} (axes directions),
#' \code{mu} (\eqn{\mu}),\code{center} (\eqn{\mathbf{c}}),
#'     \code{axesLen} (axes lengths), \code{Lambda} (\eqn{\Lambda}), \code{tau} (\eqn{\tau}) and
#'      \code{invsig} (\eqn{1 / \sigma_j^2, j = 1, \ldots, p}). For each parameter,
#'      the samples are saved as either as an array (e.g. \code{Lambda}) or as
#'      a matrix (e.g. \code{axesLen} and \code{tau}), with the first dimension being
#'      the total number of iterations.
#'
#' @export
#' @importFrom extraDistr rdirichlet
#' @importFrom extraDistr rlaplace
#' @importFrom expm expm
#' @importFrom R.utils extract
#' @importFrom narray rep
#' @importFrom reshape2 dcast melt
#' @importFrom mice mice complete
#' @importFrom ramcmc adapt_S
#' @importFrom stats rgamma
#' @examplesIf reticulate::py_module_available('ctef')
#'  parList <- gen_input_eg(shell,k = 3, updateCenter = TRUE)
#'  res <- gSGNHT_EG(shell, parList$k, parList$init, parList$prior_par,
#'                  parList$alg_par, 100, parList$updateCenter)
#' @references
#' Liu, C., Zhu, J., and Song, Y. (2016). Stochastic gradient geodesic MCMC
#' methods. In Advances in Neural Information Processing Systems 29, pages
#' 3009–3017.
gSGNHT_EG <- function(dat,k, init, prior_par, alg_par, niter, updateCenter) {
  #### initialize ####
  stopifnot(length(init$mu) == k)
  stopifnot(ncol(init$axesDir) == k)
  n <- nrow(dat)
  p <- ncol(dat)
  if(!is.null(prior_par$lambda_prior)) {
    phis <- matrix(NA, p, k)
    lambda <- matrix(NA, nrow = p, ncol = k)
    ts <- stats::rgamma(p,k* prior_par$lambda_prior, 0.5)
    para_dim <- p + p * k + 1 + k + p
    for (j in 1:p) {
    #  phis[j, ] <- MCMCpack::rdirichlet(1, base::rep(prior_par$lambda_prior, k));
      phis[j, ] <- extraDistr::rdirichlet(1, base::rep(prior_par$lambda_prior, k))
    #  lambda[j, ] <- rmutil::rlaplace(k, m=0, s= phis[j,] * ts[j])
      for (ll in 1:k) {
        lambda[j,ll] <- extraDistr::rlaplace(n = 1, mu = 0, sigma = phis[j,ll] * ts[j])
      }

    }
    kappa <- updateKap(prior_par$lambda_prior, phis, lambda) # positive - length = p
  }


  #### prior ####
  prior_para <- list(
    ltau = list(cvmf = init$tau,
                dvmf = prior_par$tau_sd,
                lowerB = prior_par$tau_lowerB,
                upperB = prior_par$tau_upperB,
                fac = 1,
                update = 1), # cvmf:mu; dvmf:standard dev
    mu = list(tau0 = 7, mu0 = init$mu),
    linvSig = list(asig = prior_par$noise_sd), # implicit: linvSig - lowerB = 0, upperB = 1
    #        axesDir = list(phi = phis, kappa = kappa, a = prior_par$lambda_prior),
    axesDir = list( a = prior_par$lambda_prior),
    axesLen = list(mu1 = 0, sig1 = 1)
    #    lambda = list(phi = phis, kappa = kappa, a = prior_par$lambda_prior)
    #center = list(cenmu = init$center, censig = prior_par$cen_sd) # experimental: Omar, removed this line, fixed center
  )
  #   lkappa = list(a = a))
  if (updateCenter) {
    prior_para$center <-  list(cenmu = init$center, censig = prior_par$cen_sd)
  }
  par_type_notEuc<- list(mu = 'sphere',axesDir = 'stiefel')
  par_type_notEuc <- par_type_notEuc[names(par_type_notEuc) %in% names(prior_para)]

  initial_val <- list(linvSig = log(init$noise_prec),
                      # axesDir = diag(1,p,k),
                      #  axesDir  = pca_$rotation[,1:k],
                      axesDir = init$axesDir,
                      ltau = tau2transformedTau(init$tau, prior_para$ltau),
                      #  ltau =  car::logit((tau_true - prior_para$ltau$lowerB) / (prior_para$ltau$upperB - prior_para$ltau$lowerB)),
                      mu = init$mu,
                      center = c(init$center),
                      # sph_r
                      axesLen =  log(init$sph_r)) # model the log of the increment

  if ('axesDir' %in% names(initial_val)) {
    initial_val$lambda <- initial_val$axesDir %*% diag(exp(initial_val$axesLen))
    # initial_val$lambda = initial_val$axesDir %*% diag(initial_val$axesLen)
  }
  if ('lambda' %in% names(prior_para)) {
    initial_val$axesDir <- NULL
    initial_val$axesLen <- NULL
  }
  #prior_para$ltau <- NULL
  # axesLen = axesLen_true)
  # axesDir = axesDir_tru,
  #   mu = c(mu_true2),
  # axesLen = c(log(censig * 2), rep(0, k - 1))
  # prior_para <- matchtheList(prior_para, par_idx)
  # initial_val <- matchtheList(initial_val, par_idx)



  #### sampling ####
  # stiefel_update: X and V associated with U, where Lambda = US.
  # , stiefel_update = NULL
  res <- gSGNHT_EG_sampling(dat, k,initial_val, prior_para,
                            par_type_notEuc,alg_par$epsilon,
                            alg_par$scalar_para, niter,
                            #NULL,
                            alg_par$minibatchSize)

  return(res)
}


#' Posterior sampler
#'
#' @description
#' `gSGNHT_EG_sampling` draws posterior samples using the geodesic stochastic
#' gradient Nose-Hoover Thermostat algorithm.
#' @param dat The data set, a numeric matrix of N by p.
#' @param latent_dim The latent dimension.
#' @param variable_interest The variables of interest.
#' @param prior_para The list of prior parameters.
#' @param par_type The parameter types.
#' @param epsilon1 The step size used in the gSGNHT sampler.
#' @param Apara The scalar parameter used in the gSGNHT sampler.
#' @param niter The total number of iterations.
#' @param minibatchSize The minibatch size.
#'
#' @returns A list of the samples, by names of the parameters.
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
gSGNHT_EG_sampling <- function(dat,
                               latent_dim,
                               variable_interest,
                               prior_para,
                               par_type,
                               epsilon1,
                               Apara,
                               niter,
                               #mode_grad,
                               minibatchSize = 50) {

  # this function uses output from gSGNHT to
  # start a sampling process with control variates for gradwaient estimation
  # dat:
  # variable_interest: initial values
  # prior_para: a list of prior parameters
  # initialized parameter (in the same order as zMat[1,])
  # tuning parameter: epsilon1, Apara, minibatchSize = 50

  #### setting up for sampling ####
  n <- nrow(dat)
  p <- ncol(dat)
  ltau_ <- ifelse(prior_para$ltau$update, 'non','ltau')
  update_par_name <- setdiff(names(prior_para), ltau_)
  para_dims <- findParDims(update_par_name, variable_interest)
  para_dim_updated <- findParaDims(para_dims)
  if(any(c('axesLen','axesDir') %in% update_par_name)) {
    par_name_lambda <- unique(c(update_par_name, 'lambda'))
  } else {
    par_name_lambda <- update_par_name
  }
  if(all(c('ltau','linvSig') %in% update_par_name)) {
    SS <- diag(0.1,  p + 1)
    #  SS <- Matrix::Diagonal(x  = rep(0.1,p + 1))
  }
  else if (c('linvSig') %in% update_par_name) {
    SS <- diag(0.1, p)
  }
  else {
    SS <- NULL
  }
  aux_par <- storeUpdates(para_dims, 1, 1)
  accurate_grad <- FALSE
  xi <- Apara
  aux_mat <- storeUpdates(para_dims, niter + 1, NA)
  grad <- storeUpdates(para_dims, niter + 1, NA)
  zMat <- storeUpdates(para_dims, niter + 1, variable_interest) # these are the ones being updated
  notUpdate <- notUpdatedPar(zMat, variable_interest)
  if ('lambda' %in% names(prior_para) ) {
    # have to CHANGE
    phis <- prior_para$lambda$phi
  }
  ########### experiments: bayesian constrain relaxation
  #if ('axesDir' %in% names(prior_para) ) {
  # have to CHANGE
  # phis <- prior_para$axesDir$phi
  # }
  #################
  covPara <- NULL
  samples <- extractSamples(zMat, 1)
  samples <- c(samples, notUpdate)
  temp <- list(par = samples, aux = aux_par, xi = xi, SS = SS, covPara = covPara)
  # SS_array <- array(NA, c(dim(SS)[1], dim(SS)[2], niter)) may need later
  acpt_ct <- 0
  mis_res <- misIdx(dat)
  dat <-  simpleImpute(dat)
  fixcenter <- FALSE
  dist <- NULL
  cur_elbow <- 0
  acpt_rate <- numeric(niter)
  if (!is.null(mis_res)) {
    predictedY <- matrix(NA, nrow = floor(niter / 2), ncol = nrow(mis_res$misPos))
  impute_idx <- 1}
  else {predictedY <- NULL}
  upweight_idx <- NULL

  #### sampling ####
  progress_bar <- utils::txtProgressBar(min=0, max=niter, style = 3, char="=")
  for (it in 2:(niter + 1)) {
    initialization_period <- FALSE

    ########## the code below is correct, but we don't need it now #######
    # if (it %% 100 == 0) {
    # print(temp$par$center)
    #  cat('Iteration', it, 'completed!\n')
    # evaluate the likelihood & create a replicate index
    # if (wt_adjust & it < 0.4 * niter) {
    #  print('adjustment!')
    #  den <- calcMargDatallh_onebyone_cpp(exp(-temp$par$linvSig), dat, temp$par$lambda,
    #  transformedTau2Tau(temp$par$ltau,
    #    prior_para$ltau$lowerB, prior_para$ltau$upperB,
    #     prior_para$ltau$fac),
    # temp$par$mu, temp$par$center)
    ## upweight_idx <- order(den, decreasing = F)[1:ceiling(0.2 * n)] # indices of to be upweighed points
    # }
    # else {upweight_idx <- NULL}
    # }
    ########################
    # if ((it < 0.3 * niter)) {
    #  initialization_period  <- T
    #}
    temp <- updateSample_gSGNHT(temp$par,par_type, para_dims, para_dim_updated,
                                temp$aux, temp$xi, dat, latent_dim, epsilon1, Apara,
                                prior_para, accurate_grad,
                                #mode_grad,
                                minibatchSize, upweight_idx, temp$SS, it - 1,
                                #temp$covPara,
                                Lsteps = 10) # (CHANGE) returning scaled para

    # grad[it,] <- temp$grad
    #  zMat[it, ] <- unlist(temp$par)
    #samples <- temp$par; aux_par <- temp$aux; xi <- temp$xi; SS <- temp$SS
    #if (controlVar) {covPara <- temp$covPara;
    #}
    #  if (initialization_period) {temp$par$linvSig <- variable_interest$linvSig; temp$par$ltau <- variable_interest$ltau}


    if (!is.null(mis_res) && it > 100) {
      dat[mis_res$misPos] <- predictY(temp$par,
                                      prior_para$ltau,
                                      dat[mis_res$na_yes,,drop = FALSE],
                                      mis_res$misIdx,
                                      mis_res$misPos_relative2Mis)
      if (it > (ceiling(niter / 2) + 1)) {
        predictedY[impute_idx,] <- dat[mis_res$misPos]
        impute_idx <- impute_idx + 1
      }
    }
    #   predicted <- predictY(temp$par, prior_para$ltau, dat[na_yes,,drop = F])  # acpt_rate
    #  predictedY[it - 1,] <- dat[na_yes, 'y'] <- predicted$y; acpt_rate[it-1] <- predicted$acpt_rate
    zMat <- assignSamples(zMat, it, temp$par, par_name_lambda)
    grad <- assignSamples(grad, it, temp$grad, update_par_name)
    # aux_mat <- assignSamples(aux_mat, it, temp$aux, update_par_name)
    ### SS_array[,,it - 1] <- temp$SS may need later
    acpt_ct <- acpt_ct + temp$accept_status
    if('lambda' %in% names(prior_para)) {
      #lambda <- matrix(zMat[it, par_idx$lambda], p, latent_dim)
      lambda <- temp$par$lambda
      psis <- updatePsi(prior_para$lambda$kappa, phis, lambda) # not neccessary, omit?
      # kappa <- updateKap(prior_para$lkappa$a, phis, lambda) # positive - length = p
      # phis <- updatePhi(prior_para$lkappa$a, lambda) # sum to 1 - keep it this way
      kappa <- updateKap(prior_para$lambda$a, phis, lambda) # positive - length = p
      phis <- updatePhi(prior_para$lambda$a, lambda) # sum to 1 - keep it this way
      prior_para$lambda$phi <- phis
      prior_para$lambda$kappa <- kappa
    }
    ########## experimental - bayesian-constrain-relaxation
    # if('axesDir' %in% names(prior_para)) {
    #   lambda <- temp$par$axesDir
    #   psis <- updatePsi(prior_para$axesDir$kappa, phis, lambda) # not neccessary, omit?
    #   kappa <- updateKap(prior_para$axesDir$a, phis, lambda) # positive - length = p
    #   phis <- updatePhi(prior_para$axesDir$a, lambda) # sum to 1 - keep it this way
    #   prior_para$axesDir$phi <- phis;
    #   prior_para$axesDir$kappa <- kappa
    # }
    ######
    utils::setTxtProgressBar(progress_bar, value = it - 1)
  }
  close(progress_bar)


  #### saving the samples ####
  res <- collectSamples(zMat, notUpdate,niter + 1, prior_para$ltau)
  # res$grad <- grad
  res$imputed <- predictedY
  #   res$acpt_rate <- acpt_ct / niter

  return(res)
}


#' Step B in the gSGNHT sampler
#'
#' @description
#' `updateV_stepB` updates the step B in the sampling procedure of gSGNHT.
#' See Algorithm 2 in the Appendix of Liu et al. (2016)
#'
#' @param epsilon The step size.
#' @param xi XXX.
#' @param v A list of the auxiliary parameters.
#'
#' @references
#' Liu, C., Zhu, J., and Song, Y. (2016). Stochastic gradient geodesic MCMC
#' methods. In Advances in Neural Information Processing Systems 29, pages
#' 3009–3017.
#' @returns A list of the updated auxiliary parameters.
updateV_stepB <- function(epsilon, xi, v) {
  # v: a list of aux par
  # res: a list of aux par
  # aux par goes to infinity
  temp <- exp(-xi * epsilon)
  res <- lapply(v, '*',temp)

  return(res)
}
#' step B in the gSGNHT sampler
#'
#' @description
#' `updateV_stepB_ctrlCov` runs step B in the sampling procedure with
#' controlled variance option. The controlled variance option is currently
#' disabled in the sampler since it does not help much with the performance.
#'
#' @param epsilon Step size.
#' @param xi Thermostat.
#' @param v The auxiliary parameter.
#' @param covContrlPar If not controlling covariance, then set to NULL; else,
#' covariance control parameter.
#'
#' @returns A list of updated V
#'
updateV_stepB_ctrlCov <- function(epsilon, xi, v, covContrlPar = NULL)  {
  # d by 1 vec
  if (!is.null(covContrlPar)) {
    sqrtN <- covContrlPar$N^2
    fac <- epsilon / 2 * sqrtN / covContrlPar$minibatchsize
    xi2 <- lapply(covContrlPar$covPara, function(x) xi + fac * x)
    #  xi <- xi + epsilon / 2 * sqrtN / covContrlPar$minibatchsize * covContrlPar$covPara # not sure whether divide by 4 or 2
    temp <- lapply(xi2, function(x) exp(-x*epsilon))
    name_ <- intersect(names(v),names(temp))
    res <- Map('*',v[name_], temp[name_])
  }
  else {
    res <- updateV_stepB(epsilon, xi, v)
  }

  return(res)
}
#' Sampling procedure for one iteration of update
#'
#' @description
#' `updateSample_gSGNHT` is the sampling procedure for one iteration of update
#' in gSGNHT.
#' See Algorithm 2 in the Appendix of Liu et al. (2016)
#' for more details.
#'
#' @param all_theta A list of the parameter values
#' @param par_type Parameter types
#' @param para_dims Parameter dimensions
#' @param para_dim_updated The dimensions of the parameters being updated.
#' @param aux_par auxiliary parameters
#' @param xi The thermostat parameter
#' @param dat A numerical matrix of p by k, the data
#' @param k The latent dimension.
#' @param epsilon step size
#' @param Apara the scalar parameter
#' @param prior_para A list of prior parameters
#' @param accurate_grad Logical, whether to use accurate or simple gradient.
#' @param minibatchSize Minibatch size.
#' @param initialization_period If NULL, then not initialization period. Upweighing
#' certain observations is disabled.
#' @param SS The covariance matrix parameter used in [proposeSigTau_ram2()].
#' @param iteration Total number of iterations.
#' @param Lsteps Number of the repetition of the ABOBA steps. Default is 10.
#'
#' @references
#' Liu, C., Zhu, J., and Song, Y. (2016). Stochastic gradient geodesic MCMC
#' methods. In Advances in Neural Information Processing Systems 29, pages
#' 3009–3017.
updateSample_gSGNHT <- function(all_theta, par_type, para_dims, para_dim_updated,
                                aux_par, xi, dat, k,
                                epsilon, Apara, prior_para, accurate_grad,
                                #mode_grad,
                                minibatchSize, initialization_period, SS = NULL,
                                iteration = NULL,
                               # covPara = NULL,
                                Lsteps = 10) {
  N <- nrow(dat)
 # if (!is.null(covPara)) {
  #  covCtrl <- T; covCtrlPar <- list(N = N, minibatchsize = minibatchSize, covPara = covPara)
#  }
#  else {
    covCtrl <- FALSE
    covCtrlPar <- NULL
 # }
  ltau_ <- ifelse(prior_para$ltau$update, 'non','ltau')
  par_names <- setdiff(names(prior_para), ltau_)

  ## experimenting updating ltau and linvSig in a AM only every 50 iterations or so##
  if (all(c('ltau','linvSig') %in% par_names)) {
    # update both tau and precision when they are present
    update_ltau_linvSig <- 2
    ##### removing ltau and linvSig from graident update
    #  par_names <- setdiff(par_names, c('ltau','linvSig'))
    # para_dim_updated <- para_dim_updated - 1 - ncol(dat)
    ##################
  }
  else if ('linvSig' %in% par_names){
    # update precision when only precision is present
    update_ltau_linvSig <- 1
  }
  else {
    # when neither are present, do not do ramcmc
    update_ltau_linvSig <- 0
  }
  par_names2 <- unique(c(par_names, 'lambda'))

  for (l in 1:Lsteps) {
    res <- updateSample_stepA(par_names, para_dim_updated, par_type,
                              all_theta, aux_par, xi, epsilon / 2)
    # v <- unlist(res$aux)
    all_theta[par_names2] <- res$par ## (possible change) double check this line
    # B.
    # cautious! testing run!
    #v <-  updateV_stepB_ctrlCov(epsilon / 2, res$xi, res$aux, NULL)
    v <-  updateV_stepB_ctrlCov(epsilon / 2, res$xi, res$aux, covCtrlPar)

    # O.

    # v <- updateV_stepO_fast(dat, k,par_idx, par_type, all_theta, v, epsilon, Apara, prior_para, accurate_grad, minibatchSize)
    v_list <- updateV_stepO_fast(dat, para_dims, k, par_type, all_theta,
                                 v, epsilon, Apara, prior_para,
                                 accurate_grad,
                                 #mode_grad,
                                 minibatchSize, initialization_period)


    v <- v_list$new_v

    # covPara
   # if (covCtrl) {
    #  covCtrlPar$covPara <- updateCovPara(covCtrlPar$covPara, iteration, v_list$gradCov)
   # }
    # B
    # cautious! testing run!
    # v <-  updateV_stepB_ctrlCov(epsilon / 2, res$xi, v, NULL)
    v <-  updateV_stepB_ctrlCov(epsilon / 2, res$xi, v, covCtrlPar)

    # A
    #aux_par <- convtVec2List(par_idx, v)
    #aux_par <- convtVec2List(par_idx[par_names], v)
    res <- updateSample_stepA(par_names, para_dim_updated,
                              par_type, all_theta, v, res$xi, epsilon / 2)
    aux_par <- res$aux
    all_theta[par_names2] <- res$par
    xi <- res$xi # (possibly change) this line

  }
  res$par <- all_theta
  res$grad <- v_list$grad # (possibly change) this line
  #cat('SG-tau', transformedTau2Tau(res$par$ltau,prior_para$ltau$lowerB, prior_para$ltau$upperB,
  #    prior_para$ltau$fac),' ')
  # add a transition kernel of ltau

  # if (all(c('ltau','linvSig','lambda') %in% par_names)) {
  #   temp <- proposeSigTauLambda_ram2(all_theta, prior_para, iteration, SS, dat)
  #   res$par$linvSig <- temp$theta$linvSig;
  #   res$par$ltau <- tau2transformedTau(temp$theta$ltau, prior_para$ltau);
  #   res$par$lambda <- temp$theta$lambda;
  #   res$SS <- temp$SS;res$accept_status <- temp$accept_status
  # }
  #  else if (all(c('ltau','linvSig','axesLen') %in% par_names)) {
  #  temp <- proposeSigTauAxes_ram(all_theta, prior_para, iteration, SS, dat)
  # res$par$ltau <- tau2transformedTau(temp$theta$ltau, prior_para$ltau);
  #res$par$ltau <- temp$theta$ltau;
  # res$par$linvSig <- temp$theta$linvSig;
  # res$par$axesLen <- temp$theta$axesLen;
  # res$SS <- temp$SS;res$accept_status <- temp$accept_status
  # }

  #    if (update_ltau_linvSig == 2) { # both tau and precision are present
  # if (iteration %% 2 == 0) {### experimenting only updating every 50 iterations $$$
  #if (all(c('ltau','linvSig') %in% par_names)) { # (possibly change) this line
  #temp <- proposeSigTauCenter(all_theta, prior_para, iteration, SS, dat)
  temp <- proposeSigTau_ram2(all_theta, prior_para, iteration, SS, dat, case = update_ltau_linvSig)
  if (!is.null(temp)) {
    if (update_ltau_linvSig == 2) {

      res$par$ltau <- tau2transformedTau(temp$theta$ltau, prior_para$ltau)
    }

    res$par$linvSig <- temp$theta$linvSig
    res$SS <- temp$SS
    res$accept_status <- temp$accept_status
  }
  # if ('mu' %in% par_names) {
  # temp <- proposeMu(res$par, prior_para$mu, prior_para$ltau,dat)
  # res$accept_status <- c(res$accept_status, temp$accept_status)
  #res$par$mu <- temp$mu
  # }
  #  }
  # else {
  #   res$SS <- SS
  # }


#  if (covCtrl) {
  #  res$covPara <- covCtrlPar$covPara
 # }
  return(res)
}

#' Update V in the gSGNHT sampling step O
#'
#' @description
#' `updateV_step0_fast` updates the auxiliary parameter V in step O of the
#' sampling procedure. See Algorithm 2 in the Appendix of Liu et al. (2016)
#' for more details.
#'
#' @param dat The minibatch, a numerical matrix of n by p.
#' @param para_dim Parameter dimensions.
#' @param k The latent dimension.
#' @param par_notin_euclidean Parameters that are not in the Euclidean space,
#' e.g. mu.
#' @param all_theta A list of the parameter values.
#' @param v The auxiliary parameter.
#' @param epsilon The step size.
#' @param Apara The scalar parameter.
#' @param prior_para A list of the prior parameters.
#' @param accurate_grad Logical, whether to use a simple or accurate gradient.
#' @param minibatchSize minibatch size.
#' @param initialization_period If NULL, then not initialization period.
#'
#' @returns A list of updated V and gradient
#' @references
#' Liu, C., Zhu, J., and Song, Y. (2016). Stochastic gradient geodesic MCMC
#' methods. In Advances in Neural Information Processing Systems 29, pages
#' 3009–3017.
updateV_stepO_fast <- function(dat,
                               para_dim,
                               k,
                               par_notin_euclidean,
                               all_theta,
                               v,
                               epsilon,
                               Apara,
                               prior_para,
                               accurate_grad,
                               #mode_grad,
                               minibatchSize,
                               initialization_period) {
  # here, we are assuming there is only one parameter not in euclidean

  #####
  # (CHANGE) just scale this one - var(theta) \geq 1 / I(\theta) - var(theta)^{-1} \leq I(\theta)
 # if (controlVar) {
  #  dUtilde <- calcSGgrad_adaptiveLangevin(all_theta, dat, k, prior_para, accurate_grad = accurate_grad,
       #                                    minibatchSize, prior_exclude = F)
#  }
 # else {
    # warning!
    dUtilde <- calcSGgrad(all_theta, dat, k, prior_para, accurate_grad = accurate_grad,
                          minibatchSize, initialization_period, prior_exclude = FALSE)
 # }

 # term1 <- if(controlVar) lapply(dUtilde$grad, '*', -epsilon) else
  term1 <- lapply(dUtilde, '*', -epsilon)

  term2 <- randomDrift_fast(Apara, epsilon, para_dim, Bmat = 0)
  names_ <- intersect(names(term1),names(term2))
  new_v <- Map('+', term1[names_], term2[names_])
  # new_v <- term1 + term2
  # new_v <- convtVec2List(par_idx[names(prior_para)], new_v)
  if (length(par_notin_euclidean) != 0) {
    if ('mu' %in% names(par_notin_euclidean)) {
      new_v$mu <- projectionTanget(all_theta$mu,new_v$mu,'sphere')
    }
    if ('axesDir' %in% names(par_notin_euclidean)) {
      #  temp <- matrix(all_theta[par_idx$lambda], nrow = ncol(dat), ncol = k)
      temp2 <- matrix(new_v$axesDir, nrow = ncol(dat), ncol = k )
      new_v$axesDir <- projectionTanget(all_theta$axesDir, temp2, 'stiefel')
    }
  }
  names_ <- intersect(names(new_v),names(v))
  new_v <- Map("+", new_v[names_], v[names_]) # have not being tested

 # if (controlVar) {
 #   return(list(new_v = new_v, grad = dUtilde$grad, gradCov = dUtilde$gradCov))
 # }
  return(list(new_v = new_v, grad = dUtilde))
}

#' Update the sample in step A
#'
#' @description
#' `updateSample_stepA` updates the samples in step A of the gSGNHT sampling
#' procedure in Algorithm 2  in the Appendix of Liu et al. (2016).
#'
#' @param par_name parameter names
#' @param para_dim parameter dimensions
#' @param par_type parameter types (in Euclidean space or not)
#' @param all_theta all the parameter values, saved in a list
#' @param aux_par auxiliary parameter V
#' @param xi The thermostats parameter
#' @param epsilon step size
#'
#' @references
#' Liu, C., Zhu, J., and Song, Y. (2016). Stochastic gradient geodesic MCMC
#' methods. In Advances in Neural Information Processing Systems 29, pages
#' 3009–3017.
updateSample_stepA <- function(par_name,  para_dim,
                               par_type, all_theta,
                               aux_par, xi, epsilon) {
  # add an additional par that contains theta&aux_par-matrix form
  # para_dim <- sum(sapply(par_idx, length))
  #  aux_par_vec <- unlist(aux_par[names(par_idx)])
  new_aux <- list()
  new_theta <- list()
  for (i in seq_along(par_name)) {
    var_ <- par_name[i]
    # (CHANGE) add if to check whether it is U and S - if so, convert the dimesnioan and call geodesic!
    if (!is.null(par_type[[var_]]) ) {
      gdf <- geodesicFlow(epsilon,all_theta[[var_]], aux_par[[var_]], type = par_type[[var_]])
    }
    else {
      gdf <- geodesicFlow(epsilon,all_theta[[var_]], aux_par[[var_]], type = 'euclidean')
    }
    new_aux[[var_]] <- gdf$v
    new_theta[[var_]] <- gdf$theta
  }

  inner_prod <- sum(sapply(aux_par, function(x) sum(x^2)))

  xi_new <- xi + (inner_prod / para_dim - 1) * epsilon


  if (any(c('axesLen','axesDir') %in% par_name)) {
    new_theta$lambda <- new_theta$axesDir %*% diag(exp(new_theta$axesLen))
    # new_theta$lambda <- new_theta$axesDir %*% diag(new_theta$axesLen)
  }
  else if (!'lambda' %in% par_name){
    new_theta$lambda <- all_theta$lambda
  }
  # add if the new par is not NULL
  res <- list(par = new_theta, aux = new_aux, xi = xi_new)
  return(res)
}

#' Transform tau
#'
#' @description
#' `transformedTau2Tau` transforms logit((tau - l) / (u - l)), the parameter
#' that is being updated to tau, where u and l are the upper and lower bounds of
#' tau respectively.
#'
#' @param par logit((tau - l) / (u - l)).
#' @param lowerB The lower bound of tau.
#' @param upperB The upper bound of tau.
#' @param fac Default is 1.
#'
transformedTau2Tau <- function(par, lowerB, upperB,fac) {
  # this function belongs to truncated tau prior case!
  # temp = logit((tau - l) / (u - l)); ltau = temp * fac
  if (is.null(upperB)) {
    return(exp(par) + lowerB)
  }
  par <- par / fac
  p <- 1 / (1 + exp(-par))
  x <- p * (upperB - lowerB) + lowerB
  return(x)
}

#' Compute stiefel flow
#'
#' @description
#' `compute_stiefel_flow` computes the stiefel flow
#'
#' @param X The orthogonal matrix, a p by k matrix \eqn{p \geq k}.
#' @param V A p by k matrix \eqn{p \geq k}.
#' @param h step size
#'
compute_stiefel_flow <- function(X, V, h) {
  # h: step size
  # V is someimtes NaN
  # print(V)
  p <- ncol(X)
  Id <- diag(p)
  Z <- matrix(0, nrow = p, ncol = p)
  A <- t(X) %*% V

  S0 <- t(V) %*% V
  temp <- h * rbind(cbind(A, -S0), cbind(Id, A))

  # stopifnot(!is.nan(A))
  temp[temp < 1e-100] <- 0
  X_V <- cbind(X, V) %*%
    expm::expm(temp) %*%
    rbind(cbind(expm::expm(-h*A), Z), cbind(Z, expm::expm(-h*A)))

  X <- X_V[,1:p, drop = FALSE]
  # browser(expr={any(is.nan(X))})
  svd_ <- svd(X)
  X <- X %*% svd_$v %*% diag(1 / (svd_$d)) %*% t(svd_$v)
  V <- X_V[,(p+1):(2*p), drop = FALSE]
  #browser(expr={all(V == 0)})
  # if p is 1 we need to convert back to matrix
  list(theta =X, v =V)
}


#' Find geodesic flow
#'
#' @description
#' `geodesicFlow` calculates the geodesic flow
#'
#' @param epsilon XXX
#' @param theta XXX
#' @param v  XXX
#' @param type The type of the surface, either 'euclidean', 'sphere' or 'stiefel'.
#'
geodesicFlow <- function(epsilon,theta, v, type) {
  if (is.null(type) || type %in% 'euclidean') {
    new_theta <- theta + v * epsilon
    res <- list(theta = new_theta, v = v)
  }
  else if (type %in% 'sphere') {
    alpha <- sqrt(sum(v^2))
    new_theta <- theta * cos(alpha * epsilon) + (v / alpha) * sin(alpha * epsilon)
    new_v <- -alpha*theta * sin(alpha * epsilon) + v * cos(alpha * epsilon)
    #  browser(expr = {any(abs(new_v) > 1000)})
    res <- list(theta = new_theta, v = new_v)
  }
  else if (type %in% 'stiefel') {
    # return:theta: p by k matrix, v: p by k matrix
    res <- compute_stiefel_flow(theta, v, epsilon)
    #   browser(expr={any(v > 1e+10)})
  }

  else {
    stop('Manifold undefined!')
  }

  return(res)
}


#' projection on the tangent of the stiefel manifold
#'
#' @description
#' `projectionTanget_stiefel` projects on the tangent surface of the
#' stiefel manifold
#'
#' @param U R_d by p
#' @param X V_d by p
#'
projectionTanget_stiefel <- function(U, X) {
  # U: R_d by p
  # X: V_d by p
  U - 0.5*X %*% (t(X) %*% U + t(U) %*% X)
}

#' Orthogonal projection
#'
#' @description
#' `orthogonalProj` finds the orthogonal projection of \code{theta} based on
#' the surface type.
#'
#' @param theta A vector, the parameter values
#' @param type A string, either 'euclidean' or 'sphere'
#'
orthogonalProj <- function(theta, type) {
  d <- length(theta)
  res <- diag(d)
  if (type %in% 'euclidean') {
  }
  else if (type %in% 'sphere') {
    res <- res - tcrossprod(theta)
  }
  else {
    stop('Manifold undefined!')
  }
  return(res)
}

#' Projection onto the tangent space
#'
#' @description
#' `projectionTanget` finds the orthogonal projection of \code{tobeprojected} based
#' on the surface type.
#'
#' @param cur_pos Current positions.
#' @param tobeprojected The surface to be projected.
#' @param type A string, either sphere or stiefel.
#'
projectionTanget <- function(cur_pos,tobeprojected, type) {
  if (type %in% 'sphere') {
    ortho_proj <- orthogonalProj(cur_pos, 'sphere')
    res <- c(ortho_proj %*% tobeprojected)
  }
  else if (type %in% 'stiefel') {
    # p by k
    res <- projectionTanget_stiefel(tobeprojected,cur_pos)
  }
  else {
    stop('Manifold undefined!')
  }
  return(res)
}

#' Draw the random drift in Step O of the sampling procedure
#'
#' @description
#' `randomDrift_fast` draws the random drift from a Gaussian distribution in step
#' O of the sampling procedure
#'
#' @param Apara The scalar parameter
#' @param epsilon step size
#' @param para_dim The parameter dimensions
#' @param Bmat V(x^*) in step O
#'
#' @returns A list of vectors, with each of length equal to the parameter
#' dimension.
#' @importFrom stats rnorm
randomDrift_fast <- function(Apara, epsilon, para_dim, Bmat = 0) {
  mean_ <- 0
  res <- vector('list', length = length(para_dim))
  names(res) <- names(para_dim)
  if (is.numeric(Bmat)) {
    cov_ <- (2*Apara - epsilon * Bmat) * epsilon
    # cov_ <- diag(cov_, para_dim)
    std_ <- sqrt(cov_)
    for (i in seq_along(res)) {
      if (length(para_dim[[i]]) == 1) {
        res[[i]] <-  stats::rnorm(para_dim[[i]], mean = mean_, sd = std_)
      }
      else {
        res[[i]] <- stats::rnorm(prod(para_dim[[i]]), mean = mean_, sd = std_)
        res[[i]] <- array(res[[i]], dim = para_dim[[i]])
      }
    }
    #  res[['ltau']] <- rnorm(1, mean = mean_, sd = sqrt(10) * std_)
  }
  else if (is.list(Bmat)) {
    thres <- 2e-5
    cov_ <- lapply(Bmat, function(x) {(2*Apara - epsilon * x) * epsilon})
    cov_ <- lapply(cov_, function(x) {x[x < 0] <- thres; return(x)})
    # print(cov_)
    for (i in seq_along(res)) {
      var_ <- names(res)[i]
      std_ <- sqrt(cov_[[var_]])
      if (length(para_dim[[i]]) == 1) {
        res[[i]] <-  stats::rnorm(para_dim[[i]], mean = mean_, sd = std_)
      }
      else {
        # not checked
        res[[i]] <- stats::rnorm(prod(para_dim[[i]]), mean = mean_, sd = std_)
        res[[i]] <- array(res[[i]], dim = para_dim[[i]])
      }
    }
  }
  else {
    stop('the input type of Bmat is not expected!')
  }# this line has not been checked yet
  #  res[['ltau']] <- rnorm(1, mean = mean_, sd = (2*Apara * 3 - epsilon * Bmat) * epsilon)
  #  term2 <- rnorm(para_dim, mean = mean_, sd = std_)
  return(res)
}


#' Find parameter dimensions
#'
#' @description
#' `findParDims` finds parameter dimensions
#'
#' @param para_names The names of the parameters
#' @param variable_interest A list
#'
#' @returns A list of parameter dimensions
#'
findParDims <- function(para_names, variable_interest) {
  # variable_interest: a list
  res <- lapply(variable_interest,dim)
  discard_names <- setdiff(names(variable_interest),para_names)
  for(i in seq_along(res)) {
    if (is.null(res[[i]])) {
      res[[i]] <- length(variable_interest[[i]])
    }
  }
  if (length(discard_names) != 0) {
    res[which(names(res) %in% discard_names)] <- NULL
    # purrr::discard(res,.p = ~stringr::str_detect(.x,discard_names))
  }
  return(res)
}

#' Format the storage for the gSGNHT output
#'
#' @description
#' `storeUpdates` creates a list by parameter names and initializes
#' a storage container for each parameter.
#'
#' @param para_dims Dimensions of parameters
#' @param niter Total number of iterations
#' @param init_vals A list of initial values by parameter names
#'
#' @returns A list for storing output of gSGNHT sampler
storeUpdates <- function(para_dims, niter, init_vals) {
  # para_dims: a list of par dims, names: par names
  # init_vals: either a scalar or a list
  res <- vector('list', length = length(para_dims))
  names(res) <- names(para_dims)
  if (niter <= 1) {
    for (i in seq_along(para_dims)) {
      if (length(para_dims[[i]]) == 1) {
        res[[i]] <- base::rep(init_vals, length = para_dims[[i]])
      }
      else {
        res[[i]] <- array(init_vals, dim = para_dims[[i]])
      }
    }
  }
  else {
    for (i in seq_along(para_dims)) {
      res[[i]] <- array(NA, dim = c(niter, para_dims[[i]]))
    }
    if (any('axesDir' %in% names(para_dims))) {
      res$lambda <- array(NA, dim = c(niter, para_dims$axesDir))
    }
    if (is.list(init_vals)) {
      for (i in seq_along(para_dims)) {
        temp <- length(dim(res[[names(para_dims)[i]]]))
        if (temp == 2) {
          res[[names(para_dims)[i]]][1,] <- init_vals[[names(para_dims[i])]]
        }
        else if (temp == 3) {
          res[[names(para_dims)[i]]][1,,] <- init_vals[[names(para_dims[i])]]
        }
        else {
          stop('The input cannot be accommodated at the moment.')
        }
      }
      if (any('axesDir' %in% names(para_dims))) {

        res$lambda[1,,] <- init_vals$axesDir %*% diag(init_vals$axesLen)
        # res$lambda[1,,] <- init_vals$axesDir %*% diag(exp(init_vals$axesLen))
      }
    }

  }
  return(res)
}

#' Extract the row indices from each element
#'
#' @description
#' `extractSamples` extracts the \code{rowIdx} from \code{sampleList}.
#'
#' @param sampleList A list of samples
#' @param rowIdx The row index of the elements that we want to extract
#'
#' @returns A list of extract samples
extractSamples <- function(sampleList, rowIdx) {
  # extract the rowIdx from each element in sampleList
  lapply(sampleList, function(x)
    R.utils::extract(x, indices=list(rowIdx), dims=c(1), drop=TRUE))
}

#' Find variable that are not updated in the sampler
#'
#' @description
#' `notUpdatedPar` extracts the variables that are not being updated in the
#'  sampler.
#'
#' @param sampleList A list of samples
#' @param all_vars Variable names
#'
#' @returns A vector, variables that are not being updated in the sampler
#'
notUpdatedPar <- function(sampleList, all_vars) {
  notUpdate <- which(!names(all_vars) %in% names(sampleList))

  return(all_vars[notUpdate])
}

#' Find the total number of parameters that are being updated
#'
#' @description
#' `findParaDims` find the total number of parameters that are being updated in
#' the sampler.
#'
#' @param para_dim_list A list of parameter dimensions
#'
#' @returns The total number of parameters that are being updated
#'
findParaDims <- function(para_dim_list) {
  sum_ <- 0
  for (i in seq_along(para_dim_list)) {
    sum_ <- sum_ + prod(para_dim_list[[i]])
  }
  return(sum_)
}

#' Save the samples
#'
#' @description
#' `collectSamples` transforms the updated parameters to their original scales
#' and saves them in a desired format.
#'
#' @param sampleList A list of the raw samples from the gSGNHT sampling
#' @param notUpdatedList A list of the parameters that are not being updated in the
#' sampling.
#' @param nrow Number of rows
#' @param tauBounds A list of tau bounds, \code{lowerB} and \code{upperB}.
#'
#' @returns A list of organized samples, on the scale of interest.
#'
collectSamples <- function(sampleList, notUpdatedList, nrow, tauBounds) {
  notUpdateName <- names(notUpdatedList)
  for (i in seq_along(notUpdatedList)) { # how about being empty?
    var_ <- notUpdateName[i]
    if (var_ %in% c('lambda','axesDir')) {
      # if (c('axesLen','axesDir') %in% notUpdatedName) {
      dim(notUpdatedList[[var_]]) <- c(1,dim(notUpdatedList[[var_]]))
      #  dim(notUpdatedList$axesDir) <- c(1,dim(notUpdatedList$lambda))
      #  sampleList$lambda <- notUpdatedList$lambda[rep(1, nrow),,]
      sampleList[[var_]] <- notUpdatedList[[var_]][base::rep(1, nrow),,]
    }
    else { #mu, tau, center
      sampleList[[var_]] <- narray::rep(as.numeric(notUpdatedList[[var_]]), nrow, along = 1)
    }
  }
  sampleList$tau <- transformedTau2Tau(sampleList$ltau,tauBounds$lowerB, tauBounds$upperB,
                                       tauBounds$fac)
  sampleList$invsig <- exp(sampleList$linvSig)
  sampleList$ltau <- NULL
  sampleList$linvSig <- NULL
  return(sampleList)
}
#' Assign samples
#'
#' @description
#' `assignSamples` assign the sample from current iteration to a list of samples.
#'
#' @param sampleList A list of samples.
#' @param rowIdx The row index.
#' @param cur_sample A list of samples from the current iteration.
#' @param par_name Parameter names.
#'
#' @returns A list of samples
assignSamples <- function(sampleList, rowIdx, cur_sample, par_name) {
  for (i in seq_along(par_name)) {
    var_ <- par_name[i]
    if (var_ %in% c('lambda','axesDir')) {
      # need to consider changing the code! (CHANGE)
      #  if (any(c('axesDir','axesLen') %in% par_name)) {
      sampleList[[var_]][rowIdx,,] <- cur_sample[[var_]]
      #  }
    }
    else {
      sampleList[[var_]][rowIdx,] <- cur_sample[[var_]]
    }
  }

  return(sampleList)
}




#' Adaptive Metropolis update for the precision and tau.
#'
#' @description
#' `proposeSigTau_ram2` uses adaptive Metropolis (Vihola 2010)
#' to update the noise variances and tau.
#'
#' @param all_theta A list of the parameters
#' @param prior_para A list of prior parameters
#' @param iteration the current iteration index
#' @param SS current covariance matrix
#' @param dat A numeric matrix, the entire data matrix
#' @param case If case is 2, both tau and variances are updated; if case is 1,
#' only variances are updated. Otherwise, neither parameters are updated.
#'
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @references
#' Vihola, M. (2010). Robust adaptive Metropolis algorithm with coerced
#' acceptance rate. Statistics and Computing, 22, 997-1008.
proposeSigTau_ram2 <- function(all_theta, prior_para, iteration, SS, dat,case) {
  # based on the orignal Tau and linvSig -
  # warning: `it risks having tau going off the edge
  # robust adaptive Metropolis (Viola 2012)
  # par, all_cur_vals, prior_para, SS, iteration, n_burnin, dat
  # SS: current covariance mat
  ## iteration: iteration index
  ## nburn_in:
  # likelihood: only the pseudo-normalising-const contains Lambda
  # prior: DE()
  # must be in this order!!
  p <- ncol(dat)
  if (case != 2 && case != 1) {
    return(NULL)
  }
  if (case == 2) {
    namess <- c('ltau','linvSig')
    cur_val <- all_theta[namess]
    cur_val$ltau <- transformedTau2Tau(all_theta$ltau, prior_para$ltau$lowerB, prior_para$ltau$upperB,
                                       prior_para$ltau$fac) # original tau

    # len <- length(cur_val)
    invalid_tau <- 1 # auto reject tau that is negative
    while (invalid_tau) {
      u <- stats::rnorm( p + 1)
      #u <- rt(len, 1)
      change <- split(c(SS %*% u), base::rep(1:2, c(1, p)))
      theta_prop <- Map("+", cur_val, change)
      names(theta_prop) <- namess
      if (prior_para$ltau$lowerB < theta_prop$ltau && prior_para$ltau$upperB > theta_prop$ltau) {
        invalid_tau <- 0
      }
    }



    #this also includes pseudo-minidataset approach

    # log scale
    var_cur <- exp(-all_theta$linvSig)
    var_prop <- exp(-theta_prop$linvSig)

    post_prop <- dellipsoidgaussian_up2const_logged(dat,var_prop,
                                                    all_theta$lambda,
                                                    theta_prop$ltau, all_theta$mu, all_theta$center) +
      # CHANGE!!!
      logDensityTau_prior(theta_prop$ltau, prior_para$ltau) +
      logDensitylinvSig(theta_prop$linvSig, prior_para$linvSig$asig)
    post_cur <- dellipsoidgaussian_up2const_logged(dat,var_cur, all_theta$lambda,
                                                   cur_val$ltau, all_theta$mu,all_theta$center) +
      # CHANGE!!!
      logDensityTau_prior(cur_val$ltau, prior_para$ltau) +
      logDensitylinvSig(all_theta$linvSig, prior_para$linvSig$asig)
  }
  else {
    # case = 1

    namess <- c('linvSig')
    cur_val <- all_theta[namess]
    u <- stats::rnorm(p)
    theta_prop <- Map("+", cur_val, list(c(SS %*% u)))
    names(theta_prop) <- namess
    orig_tau <- transformedTau2Tau(all_theta$ltau, prior_para$ltau$lowerB, prior_para$ltau$upperB,
                                   prior_para$ltau$fac) # original tau

    # log scale
    var_cur <- exp(-all_theta$linvSig)
    var_prop <- exp(-theta_prop$linvSig)

    post_prop <- dellipsoidgaussian_up2const_logged(dat,var_prop,
                                                    all_theta$lambda,
                                                    orig_tau, all_theta$mu, all_theta$center) +
      logDensitylinvSig(theta_prop$linvSig, prior_para$linvSig$asig)
    post_cur <- dellipsoidgaussian_up2const_logged(dat, var_cur,all_theta$lambda,
                                                   orig_tau, all_theta$mu,all_theta$center) +
      logDensitylinvSig(all_theta$linvSig, prior_para$linvSig$asig)
  }

  logratio <- post_prop - post_cur
  # browser(expr={is.na(logratio)})
  acceptance_prob <- min(1, exp(logratio))
  if((logratio > 0) || (log(stats::runif(1)) <= logratio)) {
    accept_status <- 1
    cur_val <- theta_prop

    # print('accept')
  }
  else { accept_status <- 0
  #print('reject')
  }
  # if(iteration <= n_burnin) {
  SS <- ramcmc::adapt_S(SS, u, acceptance_prob, iteration - 1)
  # }
  res <- list(accept_status = accept_status, theta = cur_val, SS = SS)
  return(res)
}


#' Transform tau to the parameter that is being updated
#'
#' @description
#' `tau2transformedTau` transforms tau to the parameter, namely
#'  logit((tau - l ) / (u - l)), where u and l are the tau upper and lower bounds.
#'
#' @param tau The concentration parameter
#' @param tau_bounds A list of tau bounds, named as upperB and lowerB.
#'
tau2transformedTau<- function(tau, tau_bounds) {
  if (is.null(tau_bounds$upperB)) {
    return(log(tau - tau_bounds$lowerB))
  }

  temp <- logit((tau - tau_bounds$lowerB) / (tau_bounds$upperB - tau_bounds$lowerB))

  return(temp * tau_bounds$fac)
}

#' The logit function
#'
#' @description
#' `logit` calculates \eqn{logit(p) = \log(\frac{p}{1 - p})}.
#'
#' @param p The proportion, between 0 and 1.
#'
logit <- function(p) {
  res <- log(p) - log(1 - p)
  return(res)
}
#' Evaluate the log density of the precision parameters
#'
#' @description
#' `logDensitylinvSig` evaluates the log density of the precision parameters.
#'
#' @param tau_prop Proposed tau parameter
#' @param tau_prior A list of tau prior parameters, including cvmf (mean),
#' dvmf (standard deviation), lowerB and upperB.
#' @importFrom stats dgamma
logDensityTau_prior <- function(tau_prop, tau_prior) {
  if (tau_prop < tau_prior$lowerB) {
    return(-Inf)
  }
  if ((!is.null(tau_prior$upperB)) && tau_prop > tau_prior$upperB) {
    return(-Inf)
  }
  res <- stats::dnorm(tau_prop, mean = tau_prior$cvmf, sd = tau_prior$dvmf,
                      log = TRUE)
 # res <- stats::dgamma(tau_prop, tau_prior$cvmf, tau_prior$dvmf, log = T)
  return(res)
}


#' Log-density of log(precision)
#'
#' @description
#' `logDensitylinvSig` calculates the log prior density of the log(precision),
#'  which is the parameter that is being updated in the sampler.
#'
#' @details
#' Each of the variances, 1 / precision, follows a normal distribution.
#'
#' @param linvSig log(precision).
#' @param asig The standard deviation parameter in the prior distribution.
#'
#' @returns The total contribution from the prior on the logarithm scale.
#' @importFrom stats dnorm
logDensitylinvSig <- function(linvSig, asig) {
  sum(stats::dnorm(exp(-linvSig), mean = 0, sd = asig, log = TRUE) - linvSig)
}



#' Missing index of the data entries
#'
#' @description
#' `misIdx` finds the positions for the missing entries in \code{dat}.
#'
#' @param dat a numeric matrix, data.
#'
misIdx <- function(dat) {
  na_yes <- apply(dat, 1, anyNA)
  if(sum(na_yes)) {
    misPos <- which(is.na(dat), arr.ind=TRUE)
    misPos_relative2MisPart <- which(is.na(dat[na_yes,,drop = FALSE]),
                                     arr.ind = TRUE)
    missing_idx <- as.data.frame(misPos_relative2MisPart) # missing idx relative to the subset of the data that contains missingness
    m <- reshape2::melt(missing_idx, id.var = c("row"))
    missing_idx <- reshape2::dcast(m, row ~ value, fill = 0)
    res <- as.list(apply(missing_idx[,-1,drop = FALSE],1,
                         function(x) setdiff(x, 0)))
    # names(res) <- missing_idx$row
    return(list(misPos = misPos, misPos_relative2Mis = misPos_relative2MisPart, misIdx = res, na_yes = na_yes))
  }
  return(NULL)
}

#' Impute missing data
#'
#' @description
#' `simpleImpute` imputes missing data using [mice::mice()].
#'
#' @param dat A numerical matrix of p by k, with missing entries coded as NA.
#'
#' @returns Imputed data, a numerical matrix of p by k.
simpleImpute <- function(dat) {
  imputed <- mice::mice(dat, m = 1,printFlag = FALSE)
  dat_new <- as.matrix(mice::complete(imputed,1))
  return(dat_new)
}
