#' Fit an ellipsoid-Gaussian distribution
#'
#' @description
#' `ellipsoid_gaussian` fits an ellipsoid-Gaussian distribution (or a von-Mises Fisher linear factor model)
#' to the data using geodesic Stochastic Gradient Nose-Hoover Thermostats (gSGNHT) \insertCite{liu16sgmcmc}{ellipsoidgaussian}.
#'
#' @details
#' The von-Mises Fisher linear factor model is defined by
#' \deqn{\mathbf{x} = \mathbf{c} + \Lambda \eta + \epsilon,\quad \epsilon \sim \text{N}_p({0}, \Sigma), \quad
#' \eta \sim \text{vMF}(\mu, \tau) \quad
#' \text{and} \quad  \Sigma = \text{diag}\left(\sigma_1^2, \ldots, \sigma_p^2\right).}
#'  Marginalizing over the distribution of \eqn{\eta} yields density
#' \deqn{f(\mathbf{x}) = \frac{C_k(\tau)}{(2\pi)^{\frac{p}{2}}\prod_{i = 1}^p \sigma_i}
#' \exp\left\{-\frac{1}{2} (x - {c})^T \Sigma^{-1} (x - {c})\right\}
#' \varsigma\left\{ \tau \mu + \Lambda^T \Sigma^{-1
#' }(x - {c}) , \frac{\Lambda^T \Sigma^{-1}\Lambda}{2}\right\},}
#' where \eqn{C_k(\tau) = ({\tau}/{2})^{k / 2 - 1}\{\Gamma(k / 2)I_{k /2 - 1}(\tau)\}^{-1}} and \eqn{I_{v}(\cdot)}
#' denotes the modified Bessel function of the first kind of order \eqn{v}, with respect to the Lebesgue measure on \eqn{\mathbb{R}^p,}
#' \eqn{\varsigma(\vartheta, A)} is the normalizing constant in a Fisher-Bingham distribution with density
#' \deqn{\frac{1}{\varsigma(\kappa \vartheta, A)} \exp\left(\kappa \vartheta^T{y} - {y}^T A {y}\right).}
#' We call this distribution the Ellipsoid-Gaussian (EG).
#'
#' @param dat A numeric matrix of N by p, where N is number of observations and p is the data dimension.
#' @param k The latent dimension.
#' @param scale_col Logical. If true, [scale()] is applied to the data so that data has mean 0 and standard deviation 1.
#' @param updateCenter Logical. If true, the parameter center in the ellipsoid-Gaussian distribution is updated throughout the sampler.
#' @param niter Total number of iterations of the sampler.
#' @param minibatchSize The minibatch size in the gSGNHT algorithm. The default is NULL, referring to use
#' \eqn{\max(4\%N, 50)} when \eqn{N \geq 50.}
#' @param step_size Step size in the gSGNHT algorithm. The default is \eqn{10^{-4}}; other common options are
#' \eqn{10^{-5}}, \eqn{0.5 \times 10^{-4}}, and \eqn{0.5 \times 10^{-3}}. See additional details in the Details section.
#' @param scalar_para The scalar parameter in the gSGNHT algorithm. The default is 0.1. See additional details in the Details section.
#' @returns All the samples, including the pre and post burn-in ones.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#' @importFrom Rdpack reprompt
#' @examples
#' res <- ellipsoid_gaussian(shell, 3, FALSE, TRUE, 100)
ellipsoid_gaussian <- function(dat,
                               k,
                               scale_col,
                               updateCenter,
                               niter,
                               minibatchSize = NULL,
                               step_size = 1e-4,
                               scalar_para = 0.1) {
  ####################
  #   Input:
  #  data: matrix, N * p
  #  k: number of latent dimensions (i.e. the dimensions of the vMF distribution)
  #  scale_col: True or False, whether to scale the columns to have mean 0 and var 1
  #  updateCenter: True or False, whether to update the center. If not, it is fixed at the
  #               estimate based on the Ellipsoid fitting with the Cayley transform
  #               https://arxiv.org/abs/2304.10630v1
  # niter: total number of iterations (not just post burn-in)
  #  step_size: step size in the geodesic stochastic gradient Nose-Hoover Thermostat. Common choice:
  #            (0.5 * 1e-3, 1e-4, 0.5 * 1e-4, 1e-5). Default: 1e-4
  #  scalar_para: scalar in the geodesic stochastic gradient Nose-Hoover Thermostat.
  #         Common choice: (0.1, 0.01). Default: 0.1
  #        See Alg. 2 in https://ml.cs.tsinghua.edu.cn/~changliu/sggmcmc-sam/sggmc_supp_nips2016.pdf
  #        for more info. of these two parameters
  #  output:  res - the first iteration is just the initialization
  #          axesDir: axes direction. array, (niter + 1) * p * k
  #          mu: parameter for the vMF distribution, matrix, (niter + 1) * k
  #          center: the center of the ellipsoid, matrix, (niter + 1) * p
  #          axesLen: axes length, matrix, (niter + 1) * k
  #          lambda: the factor loading matrix, array, (niter + 1) * p * k
  #          tau: parameter for the vMF distribution, matrix, (niter + 1) * 1
  #          invsig: the precision for the ellispoid-Gaussian distribution, matrix, (niter + 1) * p
  #          grad: gradients of all the parameters, a list
  #          acpt_rate: acceptance rate of tau and precisions when updated
  #                   in an additional step using robust adaptive Metropolis (Viola 2012)
  # example: res <- ellipsoid_Gaussian(dat = dat, k = 3, scale_col = F, updateCenter = T, niter = 5000)
  ####################
  if (scale_col) {
    dat <- scale(dat)
  }
  # generating the initialization
  parList <- gen_input_eg(dat,k, updateCenter,
                                       epsilon1 = step_size , scalar_para = 0.1,
                          minibatchSize = minibatchSize)
  # running geodesic SGNHT algorithm to fit an ellipsoid-Gaussian model
  res <- gSGNHT_EG(dat, parList$k, parList$init, parList$prior_par,
                   parList$alg_par, niter, parList$updateCenter)

  return(res)
}
