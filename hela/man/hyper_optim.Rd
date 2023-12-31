% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hela functions.R
\name{hyper_optim}
\alias{hyper_optim}
\title{Hyperparameter estimation of the two-level Gaussian process}
\usage{
hyper_optim(
  hypers,
  train_points_low,
  train_points_high,
  l_train_low,
  l_train_high,
  prior_mean_str
)
}
\arguments{
\item{hypers}{A vector containing the initial values of the hyperparameters,
used for optimisation, of the Gaussian process of each level:
signal variance, length scale and autoregressive parameters.}

\item{train_points_low}{A vector containing the training points of the low-level.}

\item{train_points_high}{A vector containing the training points of the high-level.}

\item{l_train_low}{A vector containing the low-level function approximation at the low-level training points.}

\item{l_train_high}{A vector containing the high-level function approximation at the high-level training points.}

\item{prior_mean_str}{A function giving the structure of the prior mean.}
}
\value{
A list containing the estimated hyperparameters, \code{optim_hypers},
and the optimisation result \code{optim_result}.
}
\description{
Compute the estimated hyperparameters for each level: signal variance, length scale parameters
and the autoregressive parameter of the two-level Gaussian process using the Nelder-Mead method of optimisation.
}
