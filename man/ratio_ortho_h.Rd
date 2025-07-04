% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multi_ortho.R
\name{ratio_ortho_h}
\alias{ratio_ortho_h}
\title{Calculate Likelihood Ratio (Ortho Method, H2/S2 Parameterization)}
\usage{
ratio_ortho_h(y, K_list, opt1, s2hat1, opt01)
}
\arguments{
\item{y}{Numeric vector of observed data for the whole Y.}

\item{K_list}{A list of numeric vectors, where each \code{K_list[[i]]} represents the
diagonal values of a transformed covariance matrix.}

\item{opt1}{Optimization result object for the alternative hypothesis (full model).
Must contain \verb{$par} for estimated parameters.}

\item{s2hat1}{Numeric scaler as estimation of s2 based on Y1.}

\item{opt01}{Optimization result object for the null hypothesis. Must contain \verb{$value} for the log-likelihood.}
}
\value{
The likelihood ratio (exp(L_theta1 - L_rho)).
}
\description{
Computes the likelihood ratio for the naive method with h2/s2 parameterization.
}
