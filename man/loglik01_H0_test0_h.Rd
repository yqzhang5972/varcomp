% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_test0.R
\name{loglik01_H0_test0_h}
\alias{loglik01_H0_test0_h}
\title{Calculate Log-Likelihood (Test for Zero Method, H2/S2 Parameterization, Null Hypothesis)}
\usage{
loglik01_H0_test0_h(par, y0, K2, i0)
}
\arguments{
\item{par}{A numeric value representing the h2 for K2.}

\item{y0}{Numeric vector of observed data for the Y0 subset.}

\item{K2}{A numeric vector representing the diagonal values of the transformed K2 matrix.}

\item{i0}{Indices used to subset the data and K matrices for Y0.}
}
\value{
The log-likelihood value.
}
\description{
This function computes the log-likelihood specifically for the null hypothesis
\eqn{h_1^2 = \dots = h_{M-1}^2 = 0} within the "test for zero" scenario, using h2/s2 parameterization.
It's for finding theta_hat under H0 when \code{K1_list} components are zero, leaving only K2 and error.
}
