% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_gradient_helpers.R
\name{calcHalfNormalPriorGrad}
\alias{calcHalfNormalPriorGrad}
\title{Truncated normal gradient}
\usage{
calcHalfNormalPriorGrad(linvSigma, logInv, asig, mean = 0)
}
\arguments{
\item{linvSigma}{A vector of length p, \eqn{\log(1 / \sigma_j^2), j = 1, \ldots,p}.}

\item{logInv}{Logical, if logInv, compute the gradient w.r.t. \eqn{\log(1 / \text{var})};
else, compute the gradient w.r.t. \eqn{\log(\text{var})}, where \eqn{\text{var}}
follows a truncated normal prior in both cases.}

\item{asig}{The standard deviation parameter in the normal distribution.}

\item{mean}{The mean parameter in the normal distribution.}
}
\value{
A vector of gradient
}
\description{
\code{calcHalfNormalPriorGrad} calculates the gradient w.r.t. \eqn{\log(1 / \sigma^2)},
with \eqn{\sigma^2} following a prior distribution N(\code{mean}, \code{asig}\eqn{^2}).
}
