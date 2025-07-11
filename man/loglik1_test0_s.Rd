% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_test0.R
\name{loglik1_test0_s}
\alias{loglik1_test0_s}
\title{Calculate Log-Likelihood (Test for Zero Method, Sigma2 Parameterization, Full Model)}
\usage{
loglik1_test0_s(par, y1, K1_list, K2, i1)
}
\arguments{
\item{par}{A numeric vector of variance component parameters. The first \code{n_K1}
elements correspond to \code{K1_list}, \code{par[length(par)-1]} corresponds to \code{K2},
and \code{par[length(par)]} is the residual error variance (s2e).}

\item{y1}{Numeric vector of observed data for the Y1 subset.}

\item{K1_list}{A list of full covariance matrices (pre-transformed by V2) corresponding to the first set of parameters.}

\item{K2}{A numeric vector representing the diagonal values of the transformed K2 matrix.}

\item{i1}{Indices used to subset the data and K matrices for Y1.}
}
\value{
The log-likelihood value.
}
\description{
This function computes the log-likelihood for the full model (Y1) in the
"test for zero" scenario, where K2 is assumed to be a single matrix and
its diagonal values are used.
}
