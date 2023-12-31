% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_eg_sampler.R
\name{notUpdatedPar}
\alias{notUpdatedPar}
\title{Find variable that are not updated in the sampler}
\usage{
notUpdatedPar(sampleList, all_vars)
}
\arguments{
\item{sampleList}{A list of samples}

\item{all_vars}{Variable names}
}
\value{
A vector, variables that are not being updated in the sampler
}
\description{
\code{notUpdatedPar} extracts the variables that are not being updated in the
sampler.
}
