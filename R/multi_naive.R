
library(matrixcalc)

#' Calculate Log-Likelihood (Naive Method, Sigma2 Parameterization, Full Model)
#'
#' This function computes the log-likelihood for the full model (Y1) using a
#' "naive" approach with full covariance matrices and a direct sigma2 parameterization
#' (s21, s22, ..., s2e).
#'
#' @param par A numeric vector of variance component parameters (s21, s22, ..., s2e).
#'   The last element is assumed to be the residual error variance (s2e).
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices corresponding to the first set of parameters.
#' @param K2_list A list of full covariance matrices corresponding to the second set of parameters.
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @return The log-likelihood value.
loglik1_naive_s <- function(par, y1, K1_list, K2_list, i1) {
  K_list <- c(K1_list, K2_list)
  n_K <- length(K_list)

  # Construct S = sum_j par[j] * K_list[[j]][i1, i1]
  S <- Reduce(`+`, lapply(1:n_K, function(j) par[j] * K_list[[j]][i1, i1]))

  # Add identity term (residual error variance)
  S <- S + par[length(par)] * diag(length(y1))  # Last element is residual variance

  # Find upper triangular Cholesky factor of S
  RS <- chol(S)

  # Compute log-likelihood
  l <- -(2*sum(log(diag(RS))) + crossprod(solve(t(RS), y1))) / 2  # -(determinant(S)$modulus[1] + crossprod(y1, chol2inv(chol(S)) %*% y1)) / 2
  return(l)
}

#' Calculate Log-Likelihood (Naive Method, Sigma2 Parameterization, Conditional Model)
#'
#' This function computes the log-likelihood for the conditional model (Y0|Y1) using a
#' "naive" approach with full covariance matrices and a direct sigma2 parameterization.
#'
#' @param par A numeric vector of variance component parameters (s21, s22, ..., s2e)
#'   for the parameters being optimized under the null hypothesis (excluding rho).
#'   The last element is assumed to be the residual error variance (s2e).
#' @param rho A numeric vector of fixed variance component parameters (s21, s22, ...)
#'   under the null hypothesis. These correspond to `K1_list`.
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices corresponding to the `rho` parameters.
#' @param K2_list A list of full covariance matrices corresponding to the `par` parameters (excluding s2e).
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @param i0 Indices used to subset the full covariance matrices for Y0.
#' @return The log-likelihood value.
loglik01_naive_s <- function(par, rho, y0, y1, K1_list, K2_list, i1, i0) { # K1_list <--> rho, K2_list <--> par to be optimized
  # Number of kernel matrices
  n_K1 <- length(K1_list)
  n_K2 <- length(K2_list)

  # Construct S00 = rho * K1[i0, i0] + sum_j par[j] * K_list[[j]][i0, i0]
  S00 <- Reduce(`+`, lapply(1:(n_K1), function(j) rho[j] * K1_list[[j]][i0, i0])) +
    Reduce(`+`, lapply(1:(n_K2), function(j) par[j] * K2_list[[j]][i0, i0])) +
    par[length(par)] * diag(length(y0))

  # Construct S01 = rho * K1[i0, i1] + sum_j par[j] * K_list[[j]][i0, i1]
  S01 <- Reduce(`+`, lapply(1:(n_K1), function(j) rho[j] * K1_list[[j]][i0, i1])) +
    Reduce(`+`, lapply(1:(n_K2), function(j) par[j] * K2_list[[j]][i0, i1]))

  # Construct S11 = rho * K1[i1, i1] + sum_j par[j] * K_list[[j]][i1, i1]
  S11 <- Reduce(`+`, lapply(1:(n_K1), function(j) rho[j] * K1_list[[j]][i1, i1])) +
    Reduce(`+`, lapply(1:(n_K2), function(j) par[j] * K2_list[[j]][i1, i1])) +
    par[length(par)] * diag(length(y1))

  # Compute conditional covariance
  # S01.11_inv <- S01 %*% chol2inv(chol(S11))
  # S <- S00 - tcrossprod(S01.11_inv, S01)

  S11inv.01t <- solve(S11, t(S01))
  S <- S00 - S01 %*% S11inv.01t

  # Symmetrize
  S <- Matrix::forceSymmetric(S)  # (S + t(S)) / 2 # as.matrix

  # Find upper triangular Cholesky factor of S
  RS <- try(chol(S), silent = T)

  # check if S is pd, matrixcalc::is.positive.definite(S)
  if(inherits(RS, "try-error")) {
    return(-1e10)
  }

  # Compute log-likelihood
  res <- y0 - crossprod(S11inv.01t, y1)  # S01.11_inv %*% y1
  l <- -(2*sum(log(diag(RS))) + crossprod(solve(t(RS), res))) / 2  # 2*sum(log(diag(chol(S)))) = determinant(S)$modulus[1], crossprod(res, chol2inv(chol(S)) %*% res)

  #print(paste("par, rho =", paste(par, collapse = ", "), rho))
  return(l)
}

#' Calculate Likelihood Ratio (Naive Method, Sigma2 Parameterization)
#'
#' Computes the likelihood ratio for the naive method with sigma2 parameterization.
#'
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices corresponding to the parameters that are fixed under the null.
#' @param K2_list A list of full covariance matrices corresponding to the parameters that are optimized.
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @param i0 Indices used to subset the full covariance matrices for Y0.
#' @param opt1 Optimization result object for the alternative hypothesis (full model).
#'   Must contain `$par` for estimated parameters.
#' @param opt01 Optimization result object for the null hypothesis. Must contain `$value` for the log-likelihood.
#' @return The likelihood ratio (exp(L_theta1 - L_rho)).
ratio_naive_s <- function(y0, y1, K1_list, K2_list, i1, i0, opt1, opt01) {
  # Compute log-likelihood under theta1
  l_theta1 <- loglik01_naive_s(par = opt1$par[-(1:length(K1_list))], rho = opt1$par[1:length(K1_list)],
                               y0 = y0, y1 = y1,
                               K1_list=K1_list, K2_list=K2_list, i1 = i1, i0 = i0)
  # Null model log-likelihood
  l_rho <- opt01$value

  return(exp(l_theta1 - l_rho))
}



#' Calculate Log-Likelihood (Naive Method, H2/S2 Parameterization, Full Model)
#'
#' This function computes the log-likelihood for the full model (Y1) using a
#' "naive" approach with full covariance matrices and an h2/s2 parameterization.
#' `par` is ordered as (h12, h22, ..., s2).
#'
#' @param par A numeric vector of parameters (h12, h22, ..., s2). The last element
#'   `par[length(par)]` is the total variance `s2`. The preceding elements are
#'   proportions of variance (`h2_i`).
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices corresponding to the first set of `h2_i` parameters.
#' @param K2_list A list of full covariance matrices corresponding to the second set of `h2_i` parameters.
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @return The log-likelihood value.
loglik1_naive_h <- function(par, y1, K1_list, K2_list, i1) {
  K_list <- c(K1_list, K2_list)
  n_K <- length(K_list)
  S <- Reduce(`+`, lapply(1:n_K, function(j) {par[j] * K_list[[j]][i1, i1]}))
  S <- par[length(par)] * (S + (1 - sum(par[-length(par)])) * diag(length(y1)))

  # Find upper triangular Cholesky factor of S
  RS <- chol(S)

  l <- -(2*sum(log(diag(RS))) + crossprod(solve(t(RS), y1))) / 2 # -(determinant(S)$modulus[1] + crossprod(y1, chol2inv(chol(S)) %*% y1)) / 2
  return(l)
}

#' Calculate Gradient (Naive Method, H2/S2 Parameterization, Full Model)
#'
#' This function computes the gradient of the log-likelihood for the full model (Y1)
#' using a "naive" approach with full covariance matrices and an h2/s2 parameterization.
#' `par` is ordered as (h12, h22, ..., s2).
#'
#' @param par A numeric vector of parameters (h12, h22, ..., s2). The last element
#'   `par[length(par)]` is the total variance `s2`. The preceding elements are
#'   proportions of variance (`h2_i`).
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices corresponding to the first set of `h2_i` parameters.
#' @param K2_list A list of full covariance matrices corresponding to the second set of `h2_i` parameters.
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @return A numeric vector of gradients corresponding to each parameter in `par`.
grad1_naive_h <- function(par, y1, K1_list, K2_list, i1) { # par = h12, h22,..., sigma2
  K_list <- c(K1_list, K2_list)
  n_K <- length(K_list)

  S <- Reduce(`+`, lapply(1:n_K, function(j) {par[j] * K_list[[j]][i1, i1]}))
  S <- par[length(par)] * (S + (1 - sum(par[-length(par)])) * diag(length(y1)))

  # Find upper triangular Cholesky factor of S
  RS <- chol(S)
  # res <- Sinv%*%y1
  res <- solve(RS, solve(t(RS), y1))

  u_vec <- sapply(1:n_K, function(j) {
    Dj <- K_list[[j]][i1, i1] - diag(length(y1))
    par[length(par)]/2 * (crossprod(res, Dj %*% res) - sum(diag(solve(RS, solve(t(RS), Dj))))) # sum(diag(Sinv%*%Dj))
  })

  usigma2 <- 1/2/par[length(par)] * (crossprod(y1, res) - length(y1))
  return(c(u_vec, usigma2))
}


#' Calculate Log-Likelihood (Naive Method, H2/S2 Parameterization, Conditional Model)
#'
#' This function computes the log-likelihood for the conditional model (Y0|Y1) using a
#' "naive" approach with full covariance matrices and an h2/s2 parameterization.
#'
#' @param par A numeric vector of parameters (h2_k, ..., h2_total, s2) for the
#'   parameters being optimized under the null hypothesis. The last element
#'   `par[length(par)]` is the total variance `s2`. The preceding elements are
#'   proportions of variance (`h2_i`) for `K2_list`.
#' @param rho A numeric vector of fixed `h2` values (h2_1, h2_2, ...) under the null
#'   hypothesis. These correspond to `K1_list`.
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices corresponding to the `rho` parameters.
#' @param K2_list A list of full covariance matrices corresponding to the `par` parameters (excluding s2).
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @param i0 Indices used to subset the full covariance matrices for Y0.
#' @return The log-likelihood value.
loglik01_naive_h <- function(par, rho, y0, y1, K1_list, K2_list, i1, i0) { # K1_list <--> rho, K2_list <--> par to be optimized
  n_K1 <- length(K1_list)
  n_K2 <- length(K2_list)

  S00 <- par[length(par)]*(Reduce(`+`, lapply(1:n_K1, function(j) {rho[j] * K1_list[[j]][i0, i0]})) +
                             (1-rho-sum(par[-length(par)]))*diag(length(y0)) +
                             Reduce(`+`, lapply(1:n_K2, function(j) {par[j] * K2_list[[j]][i0, i0]})))

  S01 <- par[length(par)]*(Reduce(`+`, lapply(1:n_K1, function(j) {rho[j] * K1_list[[j]][i0, i1]})) +
                             Reduce(`+`, lapply(1:n_K2, function(j) {par[j] * K2_list[[j]][i0, i1]})))
  S11 <- par[length(par)]*(Reduce(`+`, lapply(1:n_K1, function(j) {rho[j] * K1_list[[j]][i1, i1]})) +
                             (1-rho-sum(par[-length(par)]))*diag(length(y1)) +
                             Reduce(`+`, lapply(1:n_K2, function(j) {par[j] * K2_list[[j]][i1, i1]})))

  # S01.11_inv <- S01 %*% chol2inv(chol(S11))
  # S <- S00 - tcrossprod(S01.11_inv, S01)
  # # Symmetrize
  # S <- (S + t(S)) / 2
  # if (!matrixcalc::is.positive.definite(S)) {
  #   return(-1e10)
  # }
  # res <- y0-S01.11_inv%*%y1
  # l <- -(determinant(S)$modulus[1] +
  #          crossprod(res, chol2inv(chol(S))%*%res))/2

  S11inv.01t <- solve(S11, t(S01))
  S <- S00 - S01 %*% S11inv.01t

  # Symmetrize
  S <- Matrix::forceSymmetric(S)  # (S + t(S)) / 2 # as.matrix

  # Find upper triangular Cholesky factor of S
  RS <- try(chol(S), silent = T)

  # check if S is pd, matrixcalc::is.positive.definite(S)
  if(inherits(RS, "try-error")) {
    return(-1e10)
  }

  # Compute log-likelihood
  res <- y0 - crossprod(S11inv.01t, y1)  # S01.11_inv %*% y1
  l <- -(2*sum(log(diag(RS))) + crossprod(solve(t(RS), res))) / 2
  return(l)
}

#' Calculate Likelihood Ratio (Naive Method, H2/S2 Parameterization)
#'
#' Computes the likelihood ratio for the naive method with h2/s2 parameterization.
#'
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices corresponding to the parameters that are fixed under the null.
#' @param K2_list A list of full covariance matrices corresponding to the parameters that are optimized.
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @param i0 Indices used to subset the full covariance matrices for Y0.
#' @param opt1 Optimization result object for the alternative hypothesis (full model).
#'   Must contain `$par` for estimated parameters.
#' @param opt01 Optimization result object for the null hypothesis. Must contain `$value` for the log-likelihood.
#' @return The likelihood ratio (exp(L_theta1 - L_rho)).
ratio_naive_h <- function(y0, y1, K1_list, K2_list, i1, i0, opt1, opt01) { # y0
  l_theta1 <- loglik01_naive_h(opt1$par[-(1:length(K1_list))], opt1$par[1:length(K1_list)],
                               y0=y0, y1=y1, K1_list=K1_list, K2_list=K2_list, i1=i1, i0=i0)
  l_rho <- opt01$value
  return(exp(l_theta1 - l_rho))
}



#' Split Likelihood Ratio Test for multiple variance components (Naive Method)
#'
#' Performs a spatial likelihood ratio test using a direct, "naive" approach
#' that works with the full covariance matrices. This method is generally applicable
#' but can be computationally intensive for large datasets.
#'
#' @param y Numeric vector of observed data.
#' @param K1_list A list of full covariance matrices corresponding to parameters
#'   to be tested (fixed at `rho` under the null).
#' @param K2_list A list of full covariance matrices corresponding to nuisance parameters.
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
slrt_naive <- function(y, K1_list, K2_list, i1, i0, rho = 0, parameter.set = 'h2', TOL = 1e-8) {
  n_K1 <- length(K1_list)
  n_K2 <- length(K2_list)
  n_K <- n_K1 + n_K2
  y0 <- y[i0]
  y1 <- y[i1]

  if (parameter.set == 'h2') {
    # find consMat and consVec
    consMat <- rbind(
      c(-rep(1, n_K),0),           # h1 + h2 + ... + hJ ≤ 1 → -(h1 + ... + hJ) ≥ -1
      diag(n_K+1)               # hi ≥ 0 , s2≥ 0
    )
    consVec <- c(-1+TOL, rep(0, n_K+1))

    # find mles
    opt1 = stats::constrOptim(theta = c(rep(1/(n_K+1), n_K), 1), loglik1_naive_h, control = list(fnscale=-1),
                              grad = grad1_naive_h, method = "BFGS", ui = consMat, ci=consVec,
                              y1=y1, K1_list=K1_list, K2_list=K2_list, i1=i1)

    consMat01 <- rbind(
      c(-rep(1, n_K2),0),           # hj + ... + hJ ≤ 1-sum(rho) → -(hj + ... + hJ) ≥ -1+sum(rho)
      diag(n_K2+1)               # hi ≥ 0 , s2≥ 0
    )
    consVec01 <- c(-1+sum(rho)+TOL, rep(0, n_K2+1))
    opt01 = stats::constrOptim(theta = c(rep((1-sum(rho))/(n_K2+1), n_K2), 1), loglik01_naive_h, control = list(fnscale=-1),
                               grad = NULL, method = "Nelder-Mead", ui = consMat01, ci=consVec01,
                               rho=rho, y0=y0, y1=y1, K1_list=K1_list, K2_list=K2_list, i1=i1, i0=i0)

    teststat <- ratio_naive_h(y0=y0, y1=y1, K1_list=K1_list, K2_list=K2_list, i1=i1, i0=i0, opt1=opt1, opt01=opt01)
    return(teststat)

  } else if (parameter.set == 'sigma2') { # num of pars is n_K+1
    opt1 = stats::optim(par = rep(1, n_K+1), loglik1_naive_s, control = list(fnscale=-1),
                        y1=y1, K1_list=K1_list, K2_list=K2_list, i1=i1, method = "L-BFGS-B", lower = rep(TOL, n_K+1))

    opt01 = stats::optim(par=rep(1, n_K2+1), fn = loglik01_naive_s, control = list(fnscale=-1),
                         rho=rho, y0=y0, y1=y1, K1_list=K1_list, K2_list=K2_list,
                         i1=i1, i0=i0,
                         method = "L-BFGS-B", lower = rep(TOL, n_K2+1))
    teststat <- ratio_naive_s(y0, y1, K1_list, K2_list, i1, i0, opt1, opt01)
    return(teststat)
  }
  else {
    stop("parameter.set must be 'h2' or 'sigma2'.")
  }
}
