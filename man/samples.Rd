% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samples.R
\docType{data}
\name{samples}
\alias{samples}
\title{A small data set that is generated fitting an ellipsoid-Gaussian distribution
to the \code{shell} data set to illustrate the output of the function and
the use of the package.}
\format{
An object of class \code{"list"}
\describe{
\item{axesDir}{An array of samples of axes directions, dim = (number of
samples, p, k).}
\item{mu}{A matrix of samples of mu, number of samples by k.}
\item{axesLen}{A matrix of samples of axes lengths, number of samples by k.}
\item{lambda}{An array of samples of factor loadings, dim = (number of samples, p,k).}
\item{center}{A matrix of samples of the center, number of samples by p.}
\item{tau}{A matrix of samples of tau, number of samples by 1.}
\item{invsig}{A matrix of samples of the noise precisions, number of samples by p.}
}
}
\usage{
data(samples)
}
\description{
A small data set that is generated fitting an ellipsoid-Gaussian distribution
to the \code{shell} data set to illustrate the output of the function and
the use of the package.
}
\examples{
data(samples)
ppd <- gen_posterior_predictive(300, samples, 2500)
pairsPlot(shell, ppd)

}
\references{
This data set was created by fitting an ellipsoid-Gaussian model,
specifically via \code{\link[=ellipsoid_gaussian]{ellipsoid_gaussian()}} with \code{updateCenter = FALSE}.
}
\keyword{datasets}
