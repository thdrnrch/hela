% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hela functions.R
\name{b_profile}
\alias{b_profile}
\title{Compute the profile log-likelihood of the intercept coefficients.}
\usage{
b_profile(
  hypers,
  train_points_low,
  train_points_high,
  l_train_low,
  l_train_high,
  prior_mean_str
)
}
\arguments{
\item{hypers}{A vector containing the hyperparameters of the Gaussian process of each level:
signal variance, length scale and autoregressive parameters.}

\item{train_points_low}{A vector containing the training points of the low-level.}

\item{train_points_high}{A vector containing the training points of the high-level.}

\item{l_train_low}{A vector containing the low-level function approximation at the low-level training points.}

\item{l_train_high}{A vector containing the high-level function approximation at the high-level training points.}

\item{prior_mean_str}{A function giving the structure of the prior mean.}
}
\value{
A vector \code{b_prof_log_lik_two} containing the profile log-likelihood
at the given values of hyperparameters and training points for the two-level approximation.
}
\description{
Compute the profile log-likelihood of the intercept coefficients.
}
