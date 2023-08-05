#' Random generation from the posterior predictive distribution of an ellipsoid-Gaussian likelihood
#'
#' @description
#' `gen_posterior_predictive` generates from the posterior distribution of an ellipsoid-Gaussian likelihood,
#' using the posterior samples saved in the output format of [ellipsoid_gaussian()].
#'
#' @details
#' The posterior predictive distribution of \eqn{\tilde{x}} given \eqn{\mathbf{X}} is defined by
#' \deqn{{p(\tilde{x}|\mathbf{X}) = \int_{\Theta}p(\tilde{x}|\theta)p(\theta|\mathbf{X})d\theta}.}
#' Therefore, to generate from \eqn{p(\tilde{x}|\mathbf{X}),}
#' we generate from the posterior distribution \eqn{p(\theta|\mathbf{X})} using saved
#' samples, and then using the generated parameters to simulate from the likelihood
#' \eqn{p(\tilde{x}|\theta).}
#'
#' @param samples a list of results, same format as the output of
#' [ellipsoid_gaussian()] or [gSGNHT_EG()].
#' @param n number of samples to be generated from the posterior predictive
#' distribution.
#' @param burnin number of burn-in samples, that is the number of samples to be
#' discarded from \code{samples}.
#' @param include_noise Logical. If true, the noise (which has covariance
#' \code{Sigma}) is included; else, it is excluded and only the resulting data
#' will be on the ellipsoid.
#'
#' @returns A matrix of data, samples from the posterior predictive distribution
#' @export
#' @importFrom Rfast rvmf
#' @examples
#' gen_posterior_predictive(1000, samples, 2500)
gen_posterior_predictive <- function(n,samples,burnin,include_noise = TRUE) {
  niter <- dim(samples$lambda)[1] -1
  stopifnot(burnin < niter)
  dat <- matrix(0, nrow = n, ncol = dim(samples$lambda[1,,])[1])
  for (i in 1:n) {

    iter <- sample((burnin + 1):niter,1)
    #  Lambda <- samples$lambdas[,,iter]
    Lambda <- samples$lambda[iter,,]

    mu <- samples$mu[iter,]
    mu <- mu / sqrt(sum(mu^2))
    center <- samples$center[iter,]
    tau <- samples$tau[iter]
    latFac <- Rfast::rvmf(1, mu,tau) # n by k
    if('error_cov' %in% names(samples)) {
      # cov_ <- diag( samples$error_cov[iter,])
      cov_ <- diag( exp(samples$error_cov[iter,]))
    }
    else {
      cov_ <-  1 / samples$invsig[iter,]
      # cov_ <- diag( 1 / samples$invsig[iter,])
    }


    if (include_noise) {
      noise <- rnorm(length(cov_), mean = 0, sd = sqrt(cov_))
      # noise <- rep(0,length = nrow(samples$lambdas[,,1]))
      dat[i,] <- c(Lambda %*% t(latFac)) +  noise + center
      # dat[i,] <- sweep(latFac %*% t(Lambda) + noise ,2,center, '+')
    }
    else {
      dat[i,] <- c(Lambda %*% t(latFac)) +   center
    }

  }
  return(dat)
}
