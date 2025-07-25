% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_test0.R
\name{loglik01_H0_test0_s}
\alias{loglik01_H0_test0_s}
\title{Calculate Log-Likelihood (Test for Zero Method, Sigma2 Parameterization, Null Hypothesis)}
\usage{
loglik01_H0_test0_s(par, y0, K2, i0)
}
\arguments{
\item{par}{A numeric vector of parameters for the null hypothesis. \code{par[1]} is
the variance component for \code{K2}, and \code{par[length(par)]} is \code{s2e}.}

\item{y0}{Numeric vector of observed data for the Y0 subset.}

\item{K2}{A numeric vector representing the diagonal values of the transformed K2 matrix.}

\item{i0}{Indices used to subset the data and K matrices for Y0.}
}
\value{
The log-likelihood value.
}
\description{
This function computes the log-likelihood specifically for the null hypothesis
($h_1^2 = 0$) within the "test for zero" scenario, using sigma2 parameterization.
It's for finding theta_hat under H0, not for conditional evaluation.
}
