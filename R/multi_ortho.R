#' Calculate Log-Likelihood (Orthogonal K, Sigma2 Parameterization)
#'
#' Computes the log-likelihood of Y1 for the orthogonal K matrices case with direct
#' sigma2 parameterization (s21, s22, ..., s2e). K_list contains diagonal values.
#'
#' @param par A numeric vector of variance component parameters (s21, s22, ..., s2e).
#'   The last element is assumed to be the residual error variance (s2e).
#' @param y Numeric vector of observed data of Y1.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix, corresponding to `par[i]`.
#' @return The log-likelihood value.
loglik1_ortho_s <- function(par, y, K_list) {
  # Calculate the denominator: sum(par_i * K_i) + s2e
  deno <- Reduce(`+`, lapply(1:(length(par)-1), function(j) {par[j] * K_list[[j]]}))  + par[length(par)]
  l <- -(sum(log(deno) + y^2 / deno)) / 2
  return(l)
}


#' Calculate Log-Likelihood (Orthogonal K, Sigma2 Parameterization)
#'
#' Computes the log-likelihood of Y0 for the orthogonal K matrices case with direct
#' sigma2 parameterization (s21, s22, ..., s2e). K_list contains diagonal values.
#'
#' @param par A numeric vector of variance component parameters (s21, s22, ..., s2e)
#'   for the parameters being optimized under the null hypothesis (excluding rho).
#'   The last element is assumed to be the residual error variance (s2e).
#' @param rho A numeric vector of fixed variance component parameters (s21, s22, ...)
#'   under the null hypothesis. These correspond to `K1_list`.
#' @param y Numeric vector of observed data of Y0.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix, corresponding to `par[i]`.
#' @return The log-likelihood value.
loglik0_ortho_s <- function(par, rho, y, K_list) { # rho corresponding to K1_list
  # Calculate the denominator: sum(par_i * K_i) + s2e
  deno <- Reduce(`+`, lapply(1:(length(rho)), function(j) {rho[j] * K_list[[j]]}))
  deno <- deno + Reduce(`+`, lapply(1:(length(par)-1), function(j) {par[j] * K_list[[j+length(rho)]]}))  + par[length(par)]

  l <- -(sum(log(deno) + y^2 / deno)) / 2
  return(l)
}


#' Calculate Likelihood Ratio (Orthogonal K, Sigma2 Parameterization)
#'
#' Computes the likelihood ratio for the orthogonal method with sigma2 parameterization.
#'
#' @param y Numeric vector of observed data of Y0.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix, corresponding to `par[i]`.
#' @param opt1 Optimization result object for the alternative hypothesis (full model).
#'   Must contain `$par` for estimated parameters.
#' @param opt01 Optimization result object for the null hypothesis. Must contain `$value` for the log-likelihood.
#' @return The likelihood ratio (exp(L_theta1 - L_rho)).
ratio_ortho_s <- function(y, K_list, opt1, opt01) {
  n_K1 <- length(opt1$par) - length(opt01$par)
  # Compute log-likelihood under theta1
  l_theta1 <- loglik0_ortho_s(par = opt1$par[-(1:n_K1)], rho = opt1$par[1:n_K1],
                              y = y, K_list=K_list)
  # Null model log-likelihood
  l_rho <- opt01$value

  return(exp(l_theta1 - l_rho))
}





#' Calculate Log-Likelihood (Orthogonal K, H2/S2 Parameterization)
#'
#' Computes the log-likelihood for the orthogonal K matrices case with h2/s2 parameterization.
#' `par` is ordered as (h12, h22, ...). s2 has closed form estimation.
#'
#' @param par A numeric vector of parameters. `par[1]` to `par[length(par)-1]`
#'   are `h2_i` (proportions of variance for `K_list[[i]]`). The last parameter,
#'   `par[length(par)]`, is `s2` (total variance).
#' @param y Numeric vector of observed data Y1.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @param return.s2phat Boolean. Determine if estimation of s2 is returned, default at FALSE
#' @return The log-likelihood value.
loglik1_ortho_h <- function(par, y, K_list, return.s2phat= FALSE) {
  n_K <- length(K_list)
  deno <- Reduce(`+`, lapply(1:n_K, function(j) {par[j] * K_list[[j]]})) + 1-sum(par)           # Sigma / s2
  # if (deno < TOL) {
  #   return(-1e10)
  # }
  s2hat <- sum(y^2 / deno) / length(y)
  if(return.s2phat) {
    return(s2hat)
  }
  l <- -(sum(log(deno)) + length(y) * log(s2hat) + length(y)) / 2 # -(sum(log(deno) + y^2 / deno)) / 2
  return(l)
}

#' Calculate Log-Likelihood (Orthogonal K, H2/S2 Parameterization)
#'
#' Computes the log-likelihood for the orthogonal K matrices case with h2/s2 parameterization.
#' `par` is ordered as (h21, h22, ..., h2M).
#'
#' @param par A numeric vector of parameters (h21, h22, ..., h2M). The preceding elements are
#'   proportions of variance (`h2_i`). The total variance `s2` has closed form estimation.
#' @param rho A numeric vector of fixed `h2` values (h2_1, h2_2, ...) under the null
#'   hypothesis. These correspond to `K1_list`.
#' @param y Numeric vector of observed data Y0.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @param s2hat A numeric value of given estimation of s2. Default as NULL.
#' @return The log-likelihood value.
loglik0_ortho_h <- function(par, rho, y, K_list, s2hat = NULL) {
  deno <- Reduce(`+`, lapply(1:length(rho), function(j) {rho[j] * K_list[[j]]}))                # Sigma / s2
  deno <- deno + Reduce(`+`, lapply(1:length(par), function(j) {par[j] * K_list[[j+length(rho)]]})) + 1-sum(par)-sum(rho)

  if (is.null(s2hat)) {
    s2hat <- sum(y^2 / deno) / length(y)
    l <- -(sum(log(deno)) + length(y) * log(s2hat) + length(y)) / 2 # -(sum(log(deno) + y^2 / deno)) / 2
  } else {
    l <- -(sum(log(deno)) + length(y) * log(s2hat) + sum(y^2 / deno) / s2hat) / 2
  }
  return(l)
}


#' Calculate Gradient (Orthogonal K, H2/S2 Parameterization, for h2_i)
#'
#' Computes the gradient of the log-likelihood with respect to a specific `h2_i`
#' parameter in the orthogonal K matrices case.
#'
#' @param h2i_param_index The index of the specific `h2_i` in the `par` vector.
#' @param par The full parameter vector (h12, h22, ..., s2).
#' @param y Numeric vector of observed data.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @param TOL Numerical tolerance to be considered as zero (default: 1e-8).
#' @return The gradient value for `h2_i`.
grad_ortho_h2i <- function(h2i_param_index, par, y, K_list, TOL = 1e-8) {
  n_K <- length(K_list)

  deno <- Reduce(`+`, lapply(1:n_K, function(j) {par[j] * K_list[[j]]})) + 1-sum(par)           # Sigma / s2
  s2hat <- sum(y^2 / deno) / length(y)


  # Derivative of D with respect to h2_i:
  # dD/dh2_i = s2 * K_i - s2
  # Note: K_list[[h2i_param_index]] corresponds to the specific h2_i at that index.
  dD_dh2i <- K_list[[h2i_param_index]] - 1 # grad / s2

  # Gradient of L wrt h2_i
  g <- -sum( (1 / deno) * (1 - y^2 / deno / s2hat) * dD_dh2i ) / 2 # -(sum(dD_dh2i / deno) - length(y)) / 2 # -sum( (1 / deno) * (1 - y^2 / deno) * dD_dh2i ) / 2
  return(g)
}

#' Calculate Gradient (ortho Method, H2/S2 Parameterization, Full Model)
#'
#' This function computes the gradient of the log-likelihood for the full model (Y1)
#' using a "naive" approach with full covariance matrices and an h2/s2 parameterization.
#' `par` is ordered as (h21, h22, ..., h2M).
#'
#' @param par A numeric vector of parameters (h21, h22, ..., h2M). The preceding elements are
#'   proportions of variance (`h2_i`). The total variance `s2` has closed form estimation.
#' @param y Numeric vector of observed data for the Y1 subset.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @param TOL Numerical tolerance to be considered as zero (default: 1e-8).
#' @return A numeric vector of gradients corresponding to each parameter in `par`.
grad1_ortho_h <- function(par, y, K_list, TOL = 1e-8) { # y and K_list are subsets
  n_par <- length(par)
  gradvec <- numeric(n_par)
  for (i in 1:n_par) {
    gradvec[i] <- grad_ortho_h2i(i, par=par, y=y, K_list = K_list)
  }
  return(gradvec)
}

#' Calculate Gradient (Ortho Method, H2/S2 Parameterization, Conditional Model)
#'
#' This function computes the gradient of the log-likelihood for the conditional Model (Y0|Y1)
#' using a "naive" approach with full covariance matrices and an h2/s2 parameterization.
#' `par` is ordered as (h21, h22, ..., h2M).
#'
#' @param par A numeric vector of parameters (h21, h22, ..., h2M) for the
#'   parameters being optimized under the null hypothesis. The preceding elements are
#'   proportions of variance (`h2_i`) for `K2_list`. The the total variance `s2` has closed form estimation.
#' @param rho A numeric vector of fixed `h2` values (h2_1, h2_2, ...) under the null
#'   hypothesis. These correspond to `K1_list`.
#' @param y Numeric vector of observed data for the Y0 subset.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @param TOL Numerical tolerance to be considered as zero (default: 1e-8).
#' @return A numeric vector of gradients corresponding to each parameter in `par`.
grad0_ortho_h <- function(par, rho, y, K_list, TOL = 1e-8) { # par = h12, h22,..., h2M
  n_par <- length(par)
  gradvec <- numeric(n_par)
  for (i in 1:n_par) {
    gradvec[i] <- grad_ortho_h2i(i+length(rho), par=c(rho, par), y=y, K_list = K_list)
  }
  return(gradvec)
}


#' Calculate Likelihood Ratio (Ortho Method, H2/S2 Parameterization)
#'
#' Computes the likelihood ratio for the naive method with h2/s2 parameterization.
#'
#' @param y Numeric vector of observed data for the whole Y.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @param opt1 Optimization result object for the alternative hypothesis (full model).
#'   Must contain `$par` for estimated parameters.
#' @param s2hat1 Numeric scaler as estimation of s2 based on Y1.
#' @param opt01 Optimization result object for the null hypothesis. Must contain `$value` for the log-likelihood.
#' @return The likelihood ratio (exp(L_theta1 - L_rho)).
ratio_ortho_h <- function(y, K_list, opt1, s2hat1, opt01) {
  n_K1 <- length(opt1$par) - length(opt01$par)
  # Compute log-likelihood under theta1
  l_theta1 <- loglik0_ortho_h(par = opt1$par[-(1:n_K1)], rho = opt1$par[1:n_K1], s2hat = s2hat1,
                              y = y, K_list=K_list)
  #print(l_theta1)
  # Null model log-likelihood
  l_rho <- opt01$value
  #print(l_rho)
  return(exp(l_theta1 - l_rho))
}






#' Spatial Likelihood Ratio Test for multiple variance components (KiKj=0 method)
#'
#' Performs a spatial likelihood ratio test optimized for "orthogonal" (diagonal)
#' transformed covariance matrices. This method is highly efficient when the problem
#' can be transformed such that the covariance matrices are diagonal. It supports
#' both `sigma2` and `h2` parameter sets.
#'
#' @param y Numeric vector of observed data.
#' @param K1_list A list of numeric vectors (diagonal matrices) corresponding to parameters
#'   to be tested (fixed at `rho` under the null).
#' @param K2_list A list of numeric vectors (diagonal matrices) corresponding to nuisance parameters.
#' @param i1 Indices for the Y1 subset of data (used in the alternative hypothesis).
#' @param i0 Indices for the Y0 subset of data (used in the null hypothesis and conditional likelihood).
#' @param rho A numeric vector of fixed values for the parameters in `K1_list` under the null hypothesis.
#'   For 'h2' parameter set, these are h2 proportions. For 'sigma2', these are variance components.
#' @param parameter.set A string indicating the parameterization:
#'   'h2' for proportions of variance (h2_i) and total variance (s2), or
#'   'sigma2' for direct variance components (s21, s22, ..., s2e).
#' @param TOL Numerical tolerance for convergence and constraints (default: 1e-8).
#' @return The likelihood ratio test statistic.
#' @export
slrt_ortho <- function(y, K1_list, K2_list, i1, i0, rho, parameter.set = 'h2', TOL = 1e-8) {
  y0 <- y[i0]
  y1 <- y[i1]
  n_K1 <- length(K1_list)
  n_K2 <- length(K2_list)
  n_K <- n_K1 + n_K2

  # Estimate parameters for the full model (alternative hypothesis)
  # K_list_subset_1 will apply index1 to each matrix in K_list
  K1_list_subset_1 <- lapply(K1_list, function(K) K[i1])
  K2_list_subset_1 <- lapply(K2_list, function(K) K[i1])

  # Estimate parameters for the null model (s21 fixed at rho)
  # K_list_subset_0 will apply index0 to each matrix in K_list
  K1_list_subset_0 <- lapply(K1_list, function(K) K[i0])
  K2_list_subset_0 <- lapply(K2_list, function(K) K[i0])

  if (parameter.set == 'sigma2') {
    # Find MLEs
    opt1 <- stats::optim(par = rep(1, n_K+1), loglik1_ortho_s, control = list(fnscale=-1),
                         y=y1, K_list=c(K1_list_subset_1, K2_list_subset_1),
                         method = "L-BFGS-B", lower = rep(TOL, n_K+1))
    opt01 <- stats::optim(par = rep(1, n_K2+1), loglik0_ortho_s, control = list(fnscale=-1),
                          rho = rho, y=y0, K_list=c(K1_list_subset_0, K2_list_subset_0),
                          method = "L-BFGS-B", lower = rep(TOL, n_K2+1))
    teststat <- ratio_ortho_s(y=y0, K_list=c(K1_list_subset_0, K2_list_subset_0), opt1, opt01)

  } else if (parameter.set == 'h2') {
    # Find MLEs
    consMat <- rbind(
      -rep(1, n_K),           # h1 + h2 + ... + hJ ≤ 1 → -(h1 + ... + hJ) ≥ -1
      diag(n_K)               # hi ≥ 0
    )
    consVec <- c(-1+TOL, rep(0, n_K))
    opt1 = stats::constrOptim(theta = rep(1/(n_K+1), n_K), loglik1_ortho_h, control = list(fnscale=-1),
                              grad = grad1_ortho_h, method = "BFGS", ui = consMat, ci=consVec,
                              y=y1, K_list=c(K1_list_subset_1, K2_list_subset_1))
    s2hat1 <- loglik1_ortho_h(par=opt1$par, y=y1, K_list=c(K1_list_subset_1, K2_list_subset_1), return.s2phat= TRUE)

    consMat01 <- rbind(
      -rep(1, n_K2),           # hj + ... + hJ ≤ 1-sum(rho) → -(hj + ... + hJ) ≥ -1+sum(rho)
      diag(n_K2)               # hi ≥ 0
    )
    consVec01 <- c(-1+sum(rho)+TOL, rep(0, n_K2))
    opt01 = stats::constrOptim(theta = rep((1-sum(rho))/(n_K2+1), n_K2), loglik0_ortho_h, control = list(fnscale=-1),
                               grad = grad0_ortho_h, method = "BFGS", ui = consMat01, ci=consVec01,
                               rho=rho, y=y0, K_list=c(K1_list_subset_0, K2_list_subset_0))
    #print(opt1)
    #print(opt01)
    teststat <- ratio_ortho_h(y=y0, K_list=c(K1_list_subset_0, K2_list_subset_0), opt1, s2hat1=s2hat1, opt01)
  } else {
    stop("parameter.set must be 'h2' or 'sigma2'.")
  }

  return(teststat)
}
