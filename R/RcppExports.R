# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

RLRsimCpp <- function(p, m, n, nsim, mu, tauGrid, tau0) {
    .Call(`_varcomp_RLRsimCpp`, p, m, n, nsim, mu, tauGrid, tau0)
}

#' Restricted score test-statistic for a proportion of variation, or heritability
#'
#' @description
#' Compute a restricted score test-statistic for the proportion of variation due to the
#' variance component in a model with one variance component and an error term.
#' The function assumes the covariance matrix corresponding to the variance component is
#' diagonal which in practice usually means the actual covariance matrix
#' has been eigendecomposed and the transformed data are supplied to this
#' function (see details). This proportion is known as heritability in some
#' applications, and therefore denoted `h2` in the code.
#'
#' @param h2 The null hypothesis value, which needs to be in [0,1).
#' @param y A vector of length n of observed responses with diagonal covariance matrix.
#' @param X An n-by-p matrix of predictors, n > p.
#' @param lambda A vector of length n with (non-negative) variances, which are the
#' eigenvalues of the variance component covariance matrix (see details).
#' @param sqRoot If `true`, return statistic for signed square root statistic. Defaults to `false`.
#' @return The test-statistic evaluated at `h2`.
#'
#' @details
#' The function assumes the model
#' \deqn{y \sim N(X \beta, \sigma_g^2 \Lambda + \sigma_e^2 I_n),}
#' where \eqn{\Lambda} is a diagonal matrix with non-negative diagonal elements
#' supplied in the argument vector `lambda`. The parameter of interest is
#' \eqn{h^2=\sigma_g^2/(\sigma^2_g + \sigma^2_e)}.
#'
#' Usually this model results
#' from transformation: If
#' \deqn{\tilde{y} \sim N(\tilde{X} \beta, \sigma_g^2 K + \sigma_e^2 I_n),}
#' where \eqn{K} is a positive semi-definite (covariance) matrix with
#' eigendecomposition \eqn{U\Lambda U^\top}, then the transformed responses
#' \eqn{y = U^\top \tilde{y}} and predictors \eqn{X = U^\top \tilde{X}} satisfy
#' the model the function assumes.
#'
#' A linear mixed model with one random effect,
#' \eqn{\tilde{y} = \tilde{X}\beta + ZU + E}, where \eqn{U\sim N(0, \sigma^2_g I_q)}
#' and \eqn{E \sim N(0, \sigma^2_e I_n)}, is equivalent to the above with
#' \eqn{K = ZZ^\top}.
#'
#' The test-statistic is approximately chi-square with one degree of
#' freedom, even if `h2` is small or equal to zero, that is, near or at the
#' boundary of the parameter set. If `sqRoot = TRUE',
#' then the test-statistic is approximately standard normal.
#'
#' If the parameter of interest is instead \eqn{\tau = \sigma^2_g/\sigma^2_e},
#' note \eqn{h^2 = \tau / (1 + \tau)}, so the function can be evaluated the
#' null hypothesis value for \eqn{\tau}, say `tau`, by calling
#' `varRatioTest1d(h2 = tau / (1 + tau), ...)`.
#'
#' @export
varRatioTest1d <- function(h2, y, X, lambda, sqRoot = FALSE) {
    .Call(`_varcomp_varRatioTest1d`, h2, y, X, lambda, sqRoot)
}

#' Joint restricted score test for a proportion of variation, or heritability, and total variation
#'
#' @description
#' Compute a joint restricted score test-statistic for the total variance and the proportion
#' of variation due to the variance component in a model with one variance component and an error term.
#' See details below, and the documentation of varRatioTest1d.
#'
#' @param h2 The null hypothesis value of \eqn{h^2}, which needs to be in [0,1).
#' @param s2p The null hypothesis value of \eqn{\sigma^2_p = \sigma^2_g + \sigma^2_e}.
#' @param y A vector of length n of observed responses with diagonal covariance matrix.
#' @param X An n-by-p matrix of predictors, n > p.
#' @param lambda A vector of length n with (non-negative) variances, which are the
#' eigenvalues of the variance component covariance matrix (see details).
#' @return The test-statistic evaluated at `h2` and `s2p`.
#'
#' @details
#' The function assumes the model
#' \deqn{y \sim N(X \beta, \sigma_g^2 \Lambda + \sigma_e^2 I_n),}
#' where \eqn{\Lambda} is a diagonal matrix with non-negative diagonal elements
#' supplied in the argument vector `lambda`. See the documentation for
#' `varRatioTest1d` for how this model often results from transforming more
#' common ones using an eigendecomposition.
#'
#' The parameters in the function are \eqn{h^2=\sigma_g^2/\sigma^2_p} and
#' \eqn{\sigma^2_p = \sigma^2_g + \sigma^2_e}.
#'
#' The test-statistic is approximately chi-square with two degrees of
#' freedom, even if `h2` and `s2p` are small.
#' @export
varRatioTest2d <- function(h2, s2p, y, X, lambda) {
    .Call(`_varcomp_varRatioTest2d`, h2, s2p, y, X, lambda)
}

#' Confidence interval for a proportion of variation, or heritability
#'
#' @description
#' Calculate a confidence interval by numerically inverting the test-statistic
#' implemented in the function `varRatioTest1d`. Numerical inversion is done by
#' bisection search for points where the graph of the test-statistic as a function
#' of the null-hypothesis value `h2` crosses the appropriate quantile.
#'
#' @param range_h A vector of length 2 giving the boundaries of the interval
#' within which the bisection search is performed. The endpoints must be in [0,1).
#' @param y A vector of length n of observed responses with diagonal covariance matrix.
#' @param X An n-by-p matrix of predictors, n > p.
#' @param lambda A vector of length n with (non-negative) variances, which are the
#' eigenvalues of the variance component covariance matrix (see details).
#' @param tolerance A positive scalar with the tolerance used in bisection search.
#' @param confLevel A number in (0, 1) with the level of the confidence interval.
#' @param maxiter A positive integer. Stop and warning if number of iterations
#' in search exceeds this value.
#' @param type A string that is either "two-sided", "lower_bd" (lower bound only)
#' or "upper_bd" (upper bound only). Default is "two-sided".
#'
#' @return A vector of length 2 with endpoints of the confidence interval. NA if no root found.
#' @details
#' The function assumes the model
#' \deqn{y \sim N(X \beta, \sigma_g^2 \Lambda + \sigma_e^2 I_n),}
#' where \eqn{\Lambda} is a diagonal matrix with non-negative diagonal elements
#' supplied in the argument vector `lambda`. See the documentation for
#' `varRatioTest1d` for how this model often results from transforming more
#' common ones using an eigendecomposition.
#'
#' The parameter of interest is \eqn{h^2 = \sigma^2_g / (\sigma^2_g + \sigma^2_e).}
#'
#' If the parameter of interest is instead \eqn{\tau = \sigma^2_g/\sigma^2_e},
#' note \eqn{h^2 = \tau / (1 + \tau)}. Therefore, after running the function
#' to compute interval endpoints `a < b` for \eqn{h^2}, a confidence interval for \eqn{\tau}
#' has lower endpoint `a/(1 - a)` and upper endpoint `b / (1 - b)`.
#' @export
confInv <- function(y, X, lambda, range_h = as.numeric( c(0.0, 1.0)), tolerance = 1e-4, confLevel = 0.95, maxiter = 50L, type = "two-sided") {
    .Call(`_varcomp_confInv`, y, X, lambda, range_h, tolerance, confLevel, maxiter, type)
}

#' Joint confidence region for proportion of variation, or heritability, and total variation
#'
#' @description
#' Calculate a confidence region by numerically inverting the test-statistic
#' implemented in the function `varRatioTest2d`. Numerical inversion is done
#' by evaluating the test-statistic on a grid (see details).
#'
#' @param range_h A vector of length 2 with the boundaries for \eqn{h^2}.
#' The endpoints must be in [0,1).
#' @param range_p A vector of length 2 giving the boundaries for \eqn{\sigma^2_p}.
#' The endpoints must be positive.
#' @param y A vector of length n of observed responses with diagonal covariance matrix.
#' @param X An n-by-p matrix of predictors, n > p.
#' @param lambda A vector of length n with (non-negative) variances, which are the
#' eigenvalues of the variance component covariance matrix (see details).
#' @param grid The number of grid points in each interval, meaning the total number
#' of points in the grid is `grid^2`.
#' @return A `grid`-by-`grid` matrix with the test-statistic evaluated at the corresponding
#' grid points. Rows index `h2`, columns index `s2p`.
#' @details
#' The function assumes the model
#' \deqn{y \sim N(X \beta, \sigma_g^2 \Lambda + \sigma_e^2 I_n),}
#' where \eqn{\Lambda} is a diagonal matrix with non-negative diagonal elements
#' supplied in the argument vector `lambda`. See the documentation for
#' `varRatioTest1d` for how this model often results from transforming more
#' common ones using an eigendecomposition.
#'
#' The parameters of interest are \eqn{h^2 = \sigma^2_g / \sigma^2_p}
#' and \eqn{\sigma^2_p = \sigma^2_g + \sigma^2_e}.
#'
#' The function creates a set of feasible values for \eqn{h^2} by taking `grid`
#' evenly spaced points in the interval defined by `range_h`. It creates
#' a set of feasible values for \eqn{\sigma^2_p} by doing the same thing
#' using the interval defined by `range_p`. The Cartesian product of these
#' two sets is the grid on which the test-statistic is evaluated. A point on
#' this grid is in the confidence region if the value of the test-statistic
#' is less than the appropriate quantile of the chi-square distribution with
#' 2 degrees of freedom.
#'
#' @export
confReg <- function(y, X, lambda, range_h = as.numeric( c(0.0, 1.0)), range_p = as.numeric( c(0.0, 1.0)), grid = 200L) {
    .Call(`_varcomp_confReg`, y, X, lambda, range_h, range_p, grid)
}

