% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hela functions.R
\name{kern_gaussian}
\alias{kern_gaussian}
\title{Calculate the Gaussian kernel.}
\usage{
kern_gaussian(a, b, sigma_sq, lambda)
}
\arguments{
\item{a}{A vector containing the train/test points.}

\item{b}{A vector containing the train/test points.}

\item{sigma_sq}{Signal variance parameter.}

\item{lambda}{Length scale parameter.}
}
\value{
A matrix containing the kernel between each element of \code{a} and \code{b}.
}
\description{
Calculate the Gaussian kernel.
}
