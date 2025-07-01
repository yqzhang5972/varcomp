

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

  # # Find upper triangular Cholesky factor of S
  # RS <- chol(S)
  # Find upper triangular Cholesky factor of S
  RS <- try(chol(S), silent = T)

  # check if S is pd, matrixcalc::is.positive.definite(S)
  if(inherits(RS, "try-error")) {
    return(-1e10)
  }

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
  S <- as.matrix(Matrix::forceSymmetric(S))  # (S + t(S)) / 2 # as.matrix

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
#' `par` is ordered as (h21, h22, ..., h2M). estimation of s2 has closed form
#'
#' @param par A numeric vector of parameters (h21, h22, ..., h2M). The elements are
#'   proportions of variance (`h2_i`). The total variance `s2` has closed form estimation.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices corresponding to the first set of `h2_i` parameters.
#' @param K2_list A list of full covariance matrices corresponding to the second set of `h2_i` parameters.
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @param return.s2phat Boolean. Determine if estimation of s2 is returned, default at FALSE
#' @return The log-likelihood value.
loglik1_naive_h <- function(par, y1, K1_list, K2_list, i1, return.s2phat= FALSE) {
  n1 <- length(y1)
  K_list <- c(K1_list, K2_list)
  n_K <- length(K_list)
  S <- Reduce(`+`, lapply(1:n_K, function(j) {par[j] * K_list[[j]][i1, i1]}))
  S <- S + (1 - sum(par)) * diag(n1) # Sigma / s2

  # # Find upper triangular Cholesky factor of S
  # RS <- chol(S)
  # Find upper triangular Cholesky factor of S
  RS <- try(chol(S), silent = T)

  # check if S is pd, matrixcalc::is.positive.definite(S)
  if(inherits(RS, "try-error")) {
    return(-1e10)
  }

  # Find estimation of s2
  s2hat <- (crossprod(solve(t(RS), y1))) / n1

  if(return.s2phat) {
    return(s2hat)
  }

  l <- -(2*sum(log(diag(RS))) + n1 * log(s2hat) + n1) / 2 # -(determinant(S)$modulus[1] + crossprod(y1, chol2inv(chol(S)) %*% y1)) / 2
  return(l)
}

#' Calculate Gradient (Naive Method, H2/S2 Parameterization, Full Model)
#'
#' This function computes the gradient of the log-likelihood for the full model (Y1)
#' using a "naive" approach with full covariance matrices and an h2/s2 parameterization.
#' `par` is ordered as (h21, h22, ..., h2M).
#'
#' @param par A numeric vector of parameters (h21, h22, ..., h2M). The preceding elements are
#'   proportions of variance (`h2_i`). The total variance `s2` has closed form estimation.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices corresponding to the first set of `h2_i` parameters.
#' @param K2_list A list of full covariance matrices corresponding to the second set of `h2_i` parameters.
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @return A numeric vector of gradients corresponding to each parameter in `par`.
grad1_naive_h <- function(par, y1, K1_list, K2_list, i1) { # par = h12, h22,..., h2M
  K_list <- c(K1_list, K2_list)
  n_K <- length(K_list)
  n1 <- length(y1)

  S <- Reduce(`+`, lapply(1:n_K, function(j) {par[j] * K_list[[j]][i1, i1]}))
  S <- S + (1 - sum(par)) * diag(n1) # Sigma / s2

  # Find upper triangular Cholesky factor of S
  RS <- chol(S)

  # Find estimation of s2
  temp <- solve(t(RS), y1)
  s2hat <- (crossprod(temp)) / n1

  # res <- Sigma_inv%*%y1 * s2
  res <- solve(RS, temp)

  u_vec <- sapply(1:n_K, function(j) {
    Dj <- K_list[[j]][i1, i1] - diag(n1)
    (crossprod(res, Dj %*% res) / s2hat - sum(diag(solve(RS, solve(t(RS), Dj))))) / 2 # sum(diag(Sinv%*%Dj)) # par[length(par)]*
  })

  # usigma2 <- 1/2/par[length(par)] * (crossprod(y1, res) - length(y1))
  return(u_vec) # (c(u_vec, usigma2))
}


#' Calculate Log-Likelihood (Naive Method, H2/S2 Parameterization, Conditional Model)
#'
#' This function computes the log-likelihood for the conditional model (Y0|Y1) using a
#' "naive" approach with full covariance matrices and an h2/s2 parameterization.
#'
#' @param par A numeric vector of parameters (h21, h22, ..., h2M) for the
#'   parameters being optimized under the null hypothesis. The preceding elements are
#'   proportions of variance (`h2_i`) for `K2_list`. The the total variance `s2` has closed form estimation.
#' @param rho A numeric vector of fixed `h2` values (h2_1, h2_2, ...) under the null
#'   hypothesis. These correspond to `K1_list`.
#' @param s2hat A numeric value of given estimation of s2. Default as NULL.
#' @param y Numeric vector of whole observed data for the Y.
#' @param K1_list A list of full covariance matrices corresponding to the `rho` parameters.
#' @param K2_list A list of full covariance matrices corresponding to the `par` parameters (excluding s2).
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @param i0 Indices used to subset the full covariance matrices for Y0.
#' @return The log-likelihood value.
loglik01_naive_h <- function(par, rho, s2hat = NULL, y, K1_list, K2_list, i1, i0) { # K1_list <--> rho, K2_list <--> par to be optimized
  n_K1 <- length(K1_list)
  n_K2 <- length(K2_list)
  y0 <- y[i0]
  y1 <- y[i1]
  n0 <- length(y0)
  n1 <- length(y1)

  S00 <- Reduce(`+`, lapply(1:n_K1, function(j) {rho[j] * K1_list[[j]][i0, i0]})) +
                             (1-rho-sum(par))*diag(n0) +
                             Reduce(`+`, lapply(1:n_K2, function(j) {par[j] * K2_list[[j]][i0, i0]})) # Sigma00 / s2, sum(par[-length(par)]) -> sum(par)

  S01 <- Reduce(`+`, lapply(1:n_K1, function(j) {rho[j] * K1_list[[j]][i0, i1]})) +
                             Reduce(`+`, lapply(1:n_K2, function(j) {par[j] * K2_list[[j]][i0, i1]}))
  S11 <- Reduce(`+`, lapply(1:n_K1, function(j) {rho[j] * K1_list[[j]][i1, i1]})) +
                             (1-rho-sum(par))*diag(n1) +
                             Reduce(`+`, lapply(1:n_K2, function(j) {par[j] * K2_list[[j]][i1, i1]}))

  S11inv.01t <- solve(S11, t(S01))
  S <- S00 - S01 %*% S11inv.01t

  # Symmetrize
  S <- as.matrix(Matrix::forceSymmetric(S))  # (S + t(S)) / 2 # as.matrix

  # Find upper triangular Cholesky factor of S
  RS <- try(chol(S), silent = T)

  # check if S is pd, matrixcalc::is.positive.definite(S)
  if(inherits(RS, "try-error")) {
    return(-1e10)
  }

  res <- y0 - crossprod(S11inv.01t, y1)  # S01.11_inv %*% y1

  if (is.null(s2hat)) {                    # optimizing step
    # Find estimation of s2
    s2hat <- crossprod(solve(t(RS), res)) / n0
    # Compute log-likelihood
    l <- -(2*sum(log(diag(RS))) + n0 * log(s2hat) + n0) / 2
  } else {                                 # evaluation of theta1 step
    l <- -(2*sum(log(diag(RS))) + n0 * log(s2hat) + crossprod(solve(t(RS), res)) / s2hat) / 2
  }

  return(l)
}

#' Calculate Gradient (Naive Method, H2/S2 Parameterization, Conditional Model)
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
#' @param y Numeric vector of observed data for the whole subset.
#' @param K1_list A list of full covariance matrices corresponding to the first set of `h2_i` parameters.
#' @param K2_list A list of full covariance matrices corresponding to the second set of `h2_i` parameters.
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @param i0 Indices used to subset the full covariance matrices for Y0.
#' @return A numeric vector of gradients corresponding to each parameter in `par`.
grad01_naive_h <- function(par, rho, y, K1_list, K2_list, i1, i0) { # par = h12, h22,..., h2M

  g_all <- grad1_naive_h(par = c(rho, par), y1 = y, K1_list = K1_list, K2_list = K2_list, i1 = 1:length(y))
  g_1 <- grad1_naive_h(par = c(rho, par), y1 = y[i1], K1_list = K1_list, K2_list = K2_list, i1 = i1)

  g_01 <- g_all - g_1
  return(g_01[-(1:length(rho))]) # (c(u_vec, usigma2))
}

#' Calculate Likelihood Ratio (Naive Method, H2/S2 Parameterization)
#'
#' Computes the likelihood ratio for the naive method with h2/s2 parameterization.
#'
#' @param y Numeric vector of observed data for the whole Y.
#' @param K1_list A list of full covariance matrices corresponding to the parameters that are fixed under the null.
#' @param K2_list A list of full covariance matrices corresponding to the parameters that are optimized.
#' @param i1 Indices used to subset the full covariance matrices for Y1.
#' @param i0 Indices used to subset the full covariance matrices for Y0.
#' @param opt1 Optimization result object for the alternative hypothesis (full model).
#'   Must contain `$par` for estimated parameters.
#' @param s2hat1 Numeric scaler as estimation of s2 based on Y1.
#' @param opt01 Optimization result object for the null hypothesis. Must contain `$value` for the log-likelihood.
#' @return The likelihood ratio (exp(L_theta1 - L_rho)).
ratio_naive_h <- function(y, K1_list, K2_list, i1, i0, opt1, s2hat1, opt01) { # y0
  l_theta1 <- loglik01_naive_h(par=opt1$par[-(1:length(K1_list))], rho=opt1$par[1:length(K1_list)], s2hat = s2hat1,
                               y=y, K1_list=K1_list, K2_list=K2_list, i1=i1, i0=i0)
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
      -rep(1, n_K),           # h1 + h2 + ... + hJ ≤ 1 → -(h1 + ... + hJ) ≥ -1
      diag(n_K)               # hi ≥ 0
    )
    consVec <- c(-1+TOL, rep(0, n_K))

    # find mles
    opt1 = stats::constrOptim(theta = rep(1/(n_K+1), n_K), loglik1_naive_h, control = list(fnscale=-1),
                              grad = grad1_naive_h, method = "BFGS", ui = consMat, ci=consVec,
                              y1=y1, K1_list=K1_list, K2_list=K2_list, i1=i1)
    s2hat1 <- loglik1_naive_h(par=opt1$par, y1=y1, K1_list=K1_list, K2_list=K2_list, i1=i1, return.s2phat= TRUE)

    consMat01 <- rbind(
      -rep(1, n_K2),           # hj + ... + hJ ≤ 1-sum(rho) → -(hj + ... + hJ) ≥ -1+sum(rho)
      diag(n_K2)               # hi ≥ 0
    )
    consVec01 <- c(-1+sum(rho)+TOL, rep(0, n_K2))
    opt01 = stats::constrOptim(theta = rep((1-sum(rho))/(n_K2+1), n_K2), loglik01_naive_h, control = list(fnscale=-1),
                               grad = grad01_naive_h, method = "BFGS", ui = consMat01, ci=consVec01,
                               rho=rho, y=y, K1_list=K1_list, K2_list=K2_list, i1=i1, i0=i0)

    teststat <- ratio_naive_h(y=y, K1_list=K1_list, K2_list=K2_list, i1=i1, i0=i0, opt1=opt1, s2hat1=s2hat1, opt01=opt01)
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
