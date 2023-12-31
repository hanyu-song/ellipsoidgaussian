% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_eg_sampler.R
\name{tau2transformedTau}
\alias{tau2transformedTau}
\title{Transform tau to the parameter that is being updated}
\usage{
tau2transformedTau(tau, tau_bounds)
}
\arguments{
\item{tau}{The concentration parameter}

\item{tau_bounds}{A list of tau bounds, named as upperB and lowerB.}
}
\description{
\code{tau2transformedTau} transforms tau to the parameter, namely
logit((tau - l ) / (u - l)), where u and l are the tau upper and lower bounds.
}
