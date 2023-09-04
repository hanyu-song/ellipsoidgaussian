#' Generate the necessary input to the sampling function
#'
#' @description
#' `gen_input_eg` generates the necessary input to the sampling function
#' [gSGNHT_EG()], which includes the initial values of the Markov chain,
#' prior parameters, algorithm relevant parameters, and whether or not the center
#' should be updated during the sampling.
#'
#' @details
#' `gen_input_eg` initializes the sampler with CTEF
#'  (Melikechi and Dunson 2023). When the sample size
#'  \eqn{n \leq 50}, the minibatch size is 50; else, the minibatch size is \eqn{\max
#'  (4\%n, 50).}  While our sampler allows direct
#'  updates of \eqn{\Lambda} with a Dirichlet-Laplace prior, we find that
#'  updating the identifiable quantifies instead, namely the singular vectors and values
#'  of \eqn{\Lambda}, leads to better mixing. Therefore, we directly place a prior
#'  on the singular vectors and values.
#'
#'  The parameters in the ellipsoid-Gaussian model follow prior distributions
#'  as follows:
#' \deqn{\begin{aligned} {c} &\sim \text{N}({c}_0, \sigma_c^2 {I}_p),\\
#' \tau &\sim \text{N}(\tau_0, \sigma_{\tau}^2) \mathbb{1}(\mathfrak{l}, \mathfrak{u}),\\
#' \mu &\sim \text{vMF}(\mu_0, 7),\\
#' \sigma_{j}^2 &\sim \text{N}(0, \sigma_0^2)\mathbb{1}(0, \infty), \quad j = 1, \ldots, p,\\
#' {U} &\sim \text{Unif}(V_k(\mathbb{R}^p)),\\
#' s_l &\sim \text{N}(s_{l0}, \sigma_s^2)\mathbb{1}(0, \infty), \quad l = 1, \ldots, k,
#' \end{aligned}}
#' where
#' \eqn{V_k(\mathbb{R}^p) =\left\{A\in \mathbb{R} ^{p\times k}:A^{T}A={I}_{k}\right\}},
#' Parameters \eqn{c_0, \tau_0, \mu_0, s_{l0}} are set to be the initial values
#' produced by CTEF. The prior bounds \eqn{\mathfrak{l}} and \eqn{\mathfrak{u}} are set to be
#' \eqn{\tau_0 / 1.5} and \eqn{1.5 \tau_0,} which works well empirically.
#' We place a truncated Gaussian prior on the noise standard deviation
#' \eqn{\sigma_j} because
#' we favor a fit with small noise.
#' Results are generally not sensitive to the choice of the prior parameters.
#'
#' @param dat A numeric matrix of data of dimension \eqn{n \times p.}
#' @param k An integer, the latent dimension
#' @param updateCenter Logical, whether to update the center during the sampling.
#' @param minibatchSize Positive integer. Default is NULL, which means the
#' minibatch size is set to be \eqn{\max
#'  (4\%n, 50)}, where \eqn{n} is the sample size, for \eqn{n > 50}. When
#'  \eqn{n \leq 50}, minibatch size \eqn{= n}.
#' @param epsilon1 the step size in the gSGNHT algorithm. Default is \eqn{{10}^{-4}.}
#' See [ellipsoid_gaussian()] for more details.
#' @param scalar_para the scalar parameter in the gSGNHT algorithm. Default is 0.1.
#' See [ellipsoid_gaussian()] for more details.
#' @param noise_sd Positive scalar, the standard deviation of the noise
#' a-priori.
#' @param tau_sd Positive scalar, the standard deviation of \eqn{\tau} a-priori.
#' @param cen_sd Positive scalar, the standard deviation of the center a-priori.
#'
#' @returns A list of input for the sampler, including initial values (\code{init}),
#' k, prior parameters (\code{prior_par}), tuning parameters for the sampler
#' (\code{alg_par}), whether to update the center (\code{updateCenter}) through the
#' sampling process.
#' @export
#'
#' @importFrom mice complete
#'
#' @references
#' Melikechi, O. and Dunson, D. B. (2023). Ellipsoid fitting with the cayley
#' transform. arXiv Preprint. arXiv:2304.10630.
#'
#' @examplesIf reticulate::py_module_available('ctef')
#' gen_input_eg(shell, k = 3, TRUE)
gen_input_eg <- function(dat,
                         k,
                         updateCenter,
                         minibatchSize = NULL,
                         epsilon1 =1e-4 ,
                         scalar_para = 0.1,
                         noise_sd = 3,
                         tau_sd = 0.1,
                         cen_sd = 1) {
  # lambda_prior: parameter in dirichlet lapalce prior; if null, use uniform prior on the stiefel manifold
  # good default choice for epsilon1: 5 * 1e-4
  # 1e-5 * 5 for air quality data with missing
  if (anyNA(dat)) {
    dat <- mice::complete(mice::mice(dat,printFlag = FALSE),1)
    dat <- as.matrix(dat)
  }
  noise_var <- base::rep(1,ncol(dat)) # default
  # init <- initialization2(dat,init_tau,  noise_var,k, ordered)
  init <- init_CTEF(dat,noise_var,k,ordered = FALSE)

  # init$center <- c(-0.29922293,-0.80661124,  1.07756096 , 0.05516127, -0.20955178)
  # init$tau <- 7.326988472045936
  # init$center <- c(-2.26794017, -0.69567866, 1.34626025, -0.01169344,-0.47735721) #horse mussel
  # init$center <- c(4.9029400, 4.4831723, 2.9917971, 0.4669573, 5.1124659) # air quality
  # init$tau <- 16
  #################
  #init$sph_r <- c(3, 3)
  #init$center <- rep(0, ncol(dat))
  n <- nrow(dat)
  if (is.null(minibatchSize)) {
    minibatchSize <- floor(n / 25) # default is 1%
    if(minibatchSize < 50) { # minibatchsize is at least 50
      minibatchSize <- min(50, n)
    }
  }

  # lambda_prior is only in use if we update Lambda directly in the sampler
  # Currently the sampler updates the axes diretions and lengths, instead of
  # Lambda. Therefore, Lambda_prior is NOT in use.
  parList <- list(init = init, k = k,
                  prior_par = list(lambda_prior = 0.5, noise_sd = noise_sd,
                                   tau_sd = tau_sd,
                                   cen_sd = cen_sd, # experimental: omar-fix center
                                   tau_lowerB = init$tau / 1.5,
                                   #   tau_lowerB = max(min(findTauBound(k, 'min'),init$tau / 1.5 ),0),
                                   tau_upperB = init$tau * 1.5),
                  #      tau_upperB = max(findTauBound(k,'max'), init$tau * 1.5)),
                  alg_par = list(epsilon = epsilon1, scalar_para = scalar_para, minibatchSize = minibatchSize),
                  updateCenter = updateCenter)

  cat(paste0('Center is initialized as '), round(parList$init$center, 3),'\n')
  cat(paste0('The bound of tau is set as (',
             round(parList$prior_par$tau_lowerB, 3), ',',
             round(parList$prior_par$tau_upperB,3), ')'),'\n')



  return(parList)
}
#' Process Lambda
#'
#'@description
#'`preprocess_Lambda_ctef` processes the factor loading matrix and
#'the axes length output from [init_CTEF()] such that the order of the
#'axes lengths match that of the singular values of the factor loading matrix.
#'
#'@param Lambda The factor loading matrix output from CTEF.
#'@param ax_lengths The axes lengths output from CTEF.
#'
#'@returns A list of axes lengths and directions.
#'
preprocess_Lambda_ctef <- function(Lambda, ax_lengths) {
  # both Lambda and ax_lengths are output from ellipsoid_fit_R
  SVD <- svd(Lambda)
  k <- ncol(SVD$u)
  if (identical(diag(k), SVD$v)) {
    newU <- SVD$u
    newD <- SVD$d
  }
  else {
    #ord <- order(params$ax_lengths)
    ord <- match(round(ax_lengths,5), round(SVD$d, 5))
    newD <- SVD$d[ord]
    newU <- SVD$u[,ord]
    newV <- SVD$v[, ord]
    all_sign <- diag(newV)
    newU <- sweep(newU, 2, STATS = all_sign, FUN = '*')
  }

  # print(newU %*% diag(newD))
  stopifnot(all.equal(newU %*% diag(newD), Lambda))
  return(list(sph_r = newD, axesDir = newU))
}

#' Produce initial values to the sampler
#'
#' @description
#' `init_CTEF` provides initial values to the sampler using CTEF.
#'
#' @param dat A numeric matrix of data, N by p.
#' @param noise_var A vector of estimates of the noise variances, length = p.
#' @param k The latent dimension.
#' @param ordered Logical, whether to order the axes lengths in a monotonically
#' increasing manner. Typically we set it to be false.
#'
#' @returns A list of the initial values, ordered by parameter names.
#'
init_CTEF <- function(dat,noise_var,k,ordered) {
  params <- fit_ellipsoid_r(dat, as.integer(k))

  res <- list(center = params$center, mu = params$mu,
              noise_prec = base::rep(10, length(params$center)), tau = params$tau)
  temp <- preprocess_Lambda_ctef(params$Lambda, params$ax_lengths)
  res$axesDir <- temp$axesDir
  res$sph_r <- temp$sph_r
  if (ordered) {
    order_ <- order(res$sph_r)
    res$sph_r <- res$sph_r[order_]
    res$axesDir <- res$axesDir[,order_]
    res$mu <- res$mu[order_]
  }

  return(res)
}
#reticulate::py_run_file(system.file("python/ellipsoid_fit.py", package = "ellipsoidgaussian"))
#reticulate::source_python('inst/python/ellipsoid_fit.py')

#' Fit an ellipsoid to the noisy data using Cayley transform ellipsoid fitting (CTEF)
#'
#'
#' @description
#' `fit_ellipsoid_r` fits an ellipsoid to the noisy data and outputs relevant parameters
#'  using Cayley transform ellipsoid fitting (CTEF) (Melikechi and Dunson 2023).
#' @param X A data set, in the matrix format.
#' @param k The ellipsoid dimension (>1).

#' @returns The parameters of the fitted ellipsoid.
#' @export
#' @references
#' Melikechi, O. and Dunson, D. B. (2023). Ellipsoid fitting with the cayley
#' transform. arXiv Preprint. arXiv:2304.10630.
#'
#' @examplesIf reticulate::py_module_available('ctef')
#' fit_ellipsoid_r(shell, 3)
fit_ellipsoid_r <- function(X, k) {
  # X: data (matrix)
  # k: projection to what dimension of space?
  # output: res(list): X_full (the centered data X-\bar X times the matrix of principal components V. So X = (X-\bar X)V.)
  # inv_len (inverse of axes length), center (center)
  # s(skew-symmetric matrix), R(Caley transform of s, an orthogonal matrix)
#  mod <- reticulate::import('ctef',delay_load = TRUE)
  input <- as.matrix(X)
  res <- ctef$ctef$ctef(input,as.integer(k))
  vmf_pars <- ctef$ctef$vmf_mle(input, res$Lambda_inv, res$center)
  res$ax_lengths <- 1 / res$result$x[1:k]
  res$mu <- vmf_pars[[1]]
  res$tau <- vmf_pars[[2]]
  res$result <- NULL
  return(res)
}


