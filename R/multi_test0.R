#' Calculate Log-Likelihood (Test for Zero Method, Sigma2 Parameterization, Full Model)
#'
#' This function computes the log-likelihood for the full model (Y1) in the
#' "test for zero" scenario, where K2 is assumed to be a single matrix and
#' its diagonal values are used.
#'
#' @param par A numeric vector of variance component parameters. The first `n_K1`
#'   elements correspond to `K1_list`, `par[length(par)-1]` corresponds to `K2`,
#'   and `par[length(par)]` is the residual error variance (s2e).
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices (pre-transformed by V2) corresponding to the first set of parameters.
#' @param K2 A numeric vector representing the diagonal values of the transformed K2 matrix.
#' @param i1 Indices used to subset the data and K matrices for Y1.
#' @return The log-likelihood value.
loglik1_test0_s <- function(par, y1, K1_list, K2, i1) { # K2 is a vector of diagonal value
  n_K1 <- length(K1_list)

  # Construct S = sum_j par[j] * K_list[[j]][i1, i1]
  S <- Reduce(`+`, lapply(1:n_K1, function(j) par[j] * K1_list[[j]][i1, i1]))

  # Add identity term (residual error variance)
  S <- S + diag(par[length(par)-1] * K2[i1]  + par[length(par)] )  # Last element is residual variance

  # Find upper triangular Cholesky factor of S
  RS <- try(chol(S), silent = T)

  # check if S is pd, matrixcalc::is.positive.definite(S)
  if(inherits(RS, "try-error")) {
    return(-1e10)
  }

  # Compute log-likelihood
  l <- -(2*sum(log(diag(RS))) + crossprod(solve(t(RS), y1))) / 2 # -(determinant(S)$modulus[1] + crossprod(y1, chol2inv(chol(S)) %*% y1)) / 2
  return(l)
}



#' Calculate Log-Likelihood (Test for Zero Method, Sigma2 Parameterization, Conditional Model)
#'
#' This function computes the log-likelihood for the conditional model (Y0|Y1) in the
#' "test for zero" scenario, where K2 is assumed to be a single matrix and
#' its diagonal values are used.
#'
#' @param par A numeric vector of variance component parameters for the `K2` matrix
#'   and the residual error (`s2e`). `par[1]` is the variance component for `K2`,
#'   and `par[length(par)]` is `s2e`.
#' @param rho A numeric vector of fixed variance component parameters (s21, s22, ...)
#'   under the null hypothesis. These correspond to `K1_list`.
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices (pre-transformed by V2) corresponding to `rho`.
#' @param K2 A numeric vector representing the diagonal values of the transformed K2 matrix.
#' @param i1 Indices used to subset the data and K matrices for Y1.
#' @param i0 Indices used to subset the data and K matrices for Y0.
#' @return The log-likelihood value.
loglik01_test0_s <- function(par, rho, y0, y1, K1_list, K2, i1, i0) { # K1_list <--> rho, K2 <--> par:sigma2M, sigmae2 to be optimized
  # Number of kernel matrices
  n_K1 <- length(K1_list)

  # Construct S00 = rho * K1[i0, i0] + sum_j par[j] * K_list[[j]][i0, i0]
  S00 <- Reduce(`+`, lapply(1:n_K1, function(j) rho[j] * K1_list[[j]][i0, i0])) +
    diag(par[1] * K2[i0] + par[length(par)])

  # Construct S01 = rho * K1[i0, i1] + sum_j par[j] * K_list[[j]][i0, i1]
  S01 <- Reduce(`+`, lapply(1:n_K1, function(j) rho[j] * K1_list[[j]][i0, i1]))

  # Construct S11 = rho * K1[i1, i1] + sum_j par[j] * K_list[[j]][i1, i1]
  S11 <- Reduce(`+`, lapply(1:n_K1, function(j) rho[j] * K1_list[[j]][i1, i1])) +
    diag(par[1] * K2[i1] + par[length(par)])

  # # Compute conditional covariance
  # S01.11_inv <- S01 %*% chol2inv(chol(S11))
  # S <- S00 - tcrossprod(S01.11_inv, S01)
  # # Symmetrize
  # S <- (S + t(S)) / 2
  # if (!matrixcalc::is.positive.definite(S)) {
  #   return(-1e10)
  # }
  # # Compute log-likelihood
  # res <- y0 - S01.11_inv %*% y1
  # l <- -(determinant(S)$modulus[1] + crossprod(res, chol2inv(chol(S)) %*% res)) / 2

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
  l <- -(2*sum(log(diag(RS))) + crossprod(solve(t(RS), res))) / 2

  return(l)
}


#' Calculate Log-Likelihood (Test for Zero Method, Sigma2 Parameterization, Null Hypothesis)
#'
#' This function computes the log-likelihood specifically for the null hypothesis
#' ($h_1^2 = 0$) within the "test for zero" scenario, using sigma2 parameterization.
#' It's for finding theta_hat under H0, not for conditional evaluation.
#'
#' @param par A numeric vector of parameters for the null hypothesis. `par[1]` is
#'   the variance component for `K2`, and `par[length(par)]` is `s2e`.
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param K2 A numeric vector representing the diagonal values of the transformed K2 matrix.
#' @param i0 Indices used to subset the data and K matrices for Y0.
#' @return The log-likelihood value.
loglik01_H0_test0_s <- function(par, y0, K2, i0) { # to find theta_hat under h_1^2 = 0. NOT for evaluation of L_Y0|Y1(theta_hat1)

  # Compute conditional covariance
  S <- par[1] * K2[i0] + par[length(par)]

  # Compute log-likelihood
  l <- -(sum(log(S)) + sum(y0^2 / S)) / 2

  return(l)
}

#' Calculate Likelihood Ratio (Test for Zero Method, Sigma2 Parameterization)
#'
#' Computes the likelihood ratio for the "test for zero" method with sigma2 parameterization.
#'
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices (pre-transformed by V2) corresponding to parameters fixed under the null.
#' @param K2 A numeric vector representing the diagonal values of the transformed K2 matrix.
#' @param i1 Indices used to subset the data and K matrices for Y1.
#' @param i0 Indices used to subset the data and K matrices for Y0.
#' @param opt1 Optimization result object for the alternative hypothesis (full model).
#'   Must contain `$par` for estimated parameters.
#' @param opt01 Optimization result object for the null hypothesis. Must contain `$value` for the log-likelihood.
#' @return The likelihood ratio (exp(L_theta1 - L_rho)).
ratio_test0_s <- function(y0, y1, K1_list, K2, i1, i0, opt1, opt01) {
  # Compute log-likelihood under theta1
  l_theta1 <- loglik01_test0_s(par=opt1$par[(length(K1_list)+1):length(opt1$par)], rho=opt1$par[1:length(K1_list)],
                               y0 = y0, y1 = y1,
                               K1_list=K1_list, K2=K2, i1 = i1, i0 = i0)
  # Null model log-likelihood
  l_rho <- opt01$value

  return(exp(l_theta1 - l_rho))
}


#' Calculate Log-Likelihood (Test for Zero Method, H2/S2 Parameterization, Full Model)
#'
#' This function computes the log-likelihood for the full model (Y1) in the
#' "test for zero" scenario, using h2/s2 parameterization. K2 is assumed to be a
#' single matrix whose diagonal values are used.
#'
#' @param par A numeric vector of parameters (h21, h22, ..., h2M). The elements are
#'   proportions of variance (`h2_i`). The total variance `s2` has closed form estimation.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices (pre-transformed by V2) corresponding to the first set of `h2_i` parameters.
#' @param K2 A numeric vector representing the diagonal values of the transformed K2 matrix.
#' @param i1 Indices used to subset the data and K matrices for Y1.
#' @param return.s2phat Boolean. Determine if estimation of s2 is returned, default at FALSE
#' @return The log-likelihood value.
loglik1_test0_h <- function(par, y1, K1_list, K2, i1, return.s2phat = FALSE) { # par = h12, h22,..., h2M, K2 is a vector of diagonal values
  n1 <- length(y1)
  n_K1 <- length(K1_list)

  S <- Reduce(`+`, lapply(1:n_K1, function(j) {par[j] * K1_list[[j]][i1, i1]})) +
                             diag(par[length(par)]*K2[i1] + (1 - sum(par))) # par[length(par)] * ()

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


#' Calculate Gradient (Test for Zero Method, H2/S2 Parameterization, Full Model)
#'
#' This function computes the gradient of the log-likelihood for the full model (Y1)
#' in the "test for zero" scenario, using h2/s2 parameterization.
#'
#' @param par A numeric vector of parameters (h21, h22, ..., h2M). The preceding elements are
#'   proportions of variance (`h2_i`). The total variance `s2` has closed form estimation.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices (pre-transformed by V2) corresponding to the first set of `h2_i` parameters.
#' @param K2 A numeric vector representing the diagonal values of the transformed K2 matrix.
#' @param i1 Indices used to subset the data and K matrices for Y1.
#' @return A numeric vector of gradients corresponding to each parameter in `par`.
grad1_test0_h <- function(par, y1, K1_list, K2, i1) { # par = h12, h22,...h2M, sigma2, K2 is vector
  n_K1 <- length(K1_list)
  n1 <- length(y1)

  S <- Reduce(`+`, lapply(1:n_K1, function(j) {par[j] * K1_list[[j]][i1, i1]})) +
                             diag(par[length(par)]*K2[i1] + (1 - sum(par))) # par[length(par)] * () # Sigma / s2

  # Find upper triangular Cholesky factor of S
  RS <- chol(S)

  # Find estimation of s2
  temp <- solve(t(RS), y1)
  s2hat <- (crossprod(temp)) / n1

  # res <- Sigma_inv%*%y1 * s2
  res <- solve(RS, temp) # / s2hat

  u1_vec <- sapply(1:n_K1, function(j) {
    Dj <- K1_list[[j]][i1, i1] - diag(length(y1))
    (crossprod(res, Dj %*% res) / s2hat - sum(diag(solve(RS, solve(t(RS), Dj))))) / 2 # sum(diag(Sinv %*% Dj)))  # par[length(par)] * ()
  })

  D2 <- K2[i1] - 1
  u2 <- (crossprod(res, D2 * res) / s2hat - sum(diag(solve(RS, solve(t(RS), D2))))) / 2 # sum(diag(Sinv * D2)))   # par[length(par)]  * ()
  # usigma2 <- 1/2/par[length(par)] * (crossprod(y1, res) - length(y1))

  return(c(u1_vec, u2)) # , usigma2
}


#' Calculate Log-Likelihood (Test for Zero Method, H2/S2 Parameterization, Conditional Model)
#'
#' This function computes the log-likelihood for the conditional model (Y0|Y1) in the
#' "test for zero" scenario, using h2/s2 parameterization.
#'
#' @param par A numeric vector of parameters for (h21, h22, ..., h2M, s2). `par[1,...,M]` is
#'   the `h2` for `K2`, and `par[length(par)]` is the total variance `s2`.
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices (pre-transformed by V2) corresponding to `rho`.
#' @param K2 A numeric vector representing the diagonal values of the transformed K2 matrix.
#' @param i1 Indices used to subset the data and K matrices for Y1.
#' @param i0 Indices used to subset the data and K matrices for Y0.
#' @return The log-likelihood value.
loglik01_test0_h <- function(par, y0, y1, K1_list, K2, i1, i0) { # length of par is M+1: h21,...,h2M, s2     # previous:[par:h2M, sigma2, rho:h21,...,h2_{M-1}]
  n_K1 <- length(K1_list)

  S00 <- par[length(par)]*(Reduce(`+`, lapply(1:n_K1, function(j) {par[j] * K1_list[[j]][i0, i0]})) +
                             diag(par[n_K1+1] * K2[i0] + (1 - sum(par[-length(par)]))))
  S01 <- par[length(par)]*(Reduce(`+`, lapply(1:n_K1, function(j) {par[j] * K1_list[[j]][i0, i1]})))       # + par[1]*K2[i0, i1])
  S11 <- par[length(par)]*(Reduce(`+`, lapply(1:n_K1, function(j) {par[j] * K1_list[[j]][i1, i1]})) +
                             diag(par[n_K1+1] * K2[i1] + (1 - sum(par[-length(par)]))))
  # # print(paste("last eigen:", eigen(S11)$values[dim(S11)[1]]))
  # S01.11_inv <- S01 %*% chol2inv(chol(S11))
  # S <- S00 - tcrossprod(S01.11_inv, S01)
  # l <- -(determinant(S)$modulus[1] +
  #          crossprod(y0-S01.11_inv%*%y1, chol2inv(chol(S))%*%(y0-S01.11_inv%*%y1))) / 2

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
  l <- -(2*sum(log(diag(RS))) + crossprod(solve(t(RS), res))) / 2
  return(l)
}

#' Calculate Log-Likelihood (Test for Zero Method, H2/S2 Parameterization, Null Hypothesis)
#'
#' This function computes the log-likelihood specifically for the null hypothesis
#' ($h_1^2 =...=h_{M-1}^2= 0$) within the "test for zero" scenario, using h2/s2 parameterization.
#' It's for finding theta_hat under H0 when `K1_list` components are zero, leaving only K2 and error.
#'
#' @param par A numeric value representing the h2 for K2.
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param K2 A numeric vector representing the diagonal values of the transformed K2 matrix.
#' @param i0 Indices used to subset the data and K matrices for Y0.
#' @return The log-likelihood value.
loglik01_H0_test0_h <- function(par, y0, K2, i0) { # to find theta_hat under h_1^2 = 0. NOT for evaluation of L_Y0|Y1(theta_hat1)
  dseq <- par * K2[i0] + 1-par
  s2hat <- sum(y0^2 / dseq) / length(y0)
  l <- -(length(y0)*log(s2hat) + sum(log(dseq)) + length(y0)) / 2  #sum(log(dseq) + y0^2/s2/dseq)) / 2
  return(l)
}

#' Calculate Likelihood Ratio (Test for Zero Method, H2/S2 Parameterization)
#'
#' Computes the likelihood ratio for the "test for zero" method with h2/s2 parameterization.
#'
#' @param y0 Numeric vector of observed data for the Y0 subset.
#' @param y1 Numeric vector of observed data for the Y1 subset.
#' @param K1_list A list of full covariance matrices (pre-transformed by V2) corresponding to parameters fixed under the null.
#' @param K2 A numeric vector representing the diagonal values of the transformed K2 matrix.
#' @param i1 Indices used to subset the data and K matrices for Y1.
#' @param i0 Indices used to subset the data and K matrices for Y0.
#' @param opt1 Optimization result object for the alternative hypothesis (full model).
#'   Must contain `$par` for estimated parameters.
#' @param s2hat1 Numeric scaler as estimation of s2 based on Y1.
#' @param opt01 Optimization result object for the null hypothesis. Must contain `$value` for the log-likelihood.
#' @return The likelihood ratio (exp(L_theta1 - L_rho)).
ratio_test0_h <- function(y0, y1, K1_list, K2, i1, i0, opt1, s2hat1, opt01) { # K2 is vector
  l_theta1 <- loglik01_test0_h(c(opt1$par, s2hat1),  # opt1$par[(length(K1_list)+1):length(opt1$par)], opt1$par[1:length(K1_list)],
                               y0=y0, y1=y1, K1_list=K1_list, K2=K2, i1=i1, i0=i0)
  l_rho <- opt01$value
  return(exp(l_theta1 - l_rho))
}



#' Split Likelihood Ratio Test for multiple variance components (Test for Zero Method)
#'
#' Conducts a specialized likelihood ratio test to test for a specific genetic component
#' being zero (e.g., h12 = 0). This function is tailored for scenarios where a single
#' "nuisance" matrix (`K2_list[[1]]`) is provided, and it leverages an eigen-decomposition
#' of this matrix for potentially faster computations. It supports both `sigma2` and `h2` parameter sets.
#'
#' @param y Numeric vector of observed data.
#' @param K1_list A list of full covariance matrices corresponding to parameters
#'   to be tested (fixed at 0 under the null). These will be transformed by V2.
#' @param K2_list A list containing a single full covariance matrix (K2) which is the nuisance parameter.
#'   This matrix will be eigen-decomposed.
#' @param i1 Indices for the Y1 subset of data (used in the alternative hypothesis).
#' @param i0 Indices for the Y0 subset of data (used in the null hypothesis and conditional likelihood).
#' @param parameter.set A string indicating the parameterization:
#'   'h2' for proportions of variance (`h2_i`) and total variance (s2), or
#'   'sigma2' for direct variance components (s21, s22, ..., s2e).
#' @param TOL Numerical tolerance for convergence and constraints (default: 1e-8).
#' @return The likelihood ratio test statistic.
#' @export
slrt_test0 <- function(y, K1_list, K2_list, i1, i0, parameter.set = 'h2', TOL = 1e-8) { # rho is 0
  n_K1 <- length(K1_list)
  if (length(K2_list) > 1) {
    stop("This method can only deal with one nuisance matrix.")
  }

  eoK2 <- eigen(K2_list[[1]])
  V2 <- eoK2$vec
  K1_list <- lapply(1:n_K1, function(j) {crossprod(V2, K1_list[[j]]) %*% V2})
  K2 <- eoK2$val

  y <- crossprod(V2, y)
  y0 <- y[i0]
  y1 <- y[i1]

  if (parameter.set == 'h2') {
    # find consMat and consVec
    consMat <- rbind(
      -rep(1, n_K1+1),           # h1 + h2 + ... + hJ ≤ 1 → -(h1 + ... + hJ) ≥ -1
      diag(n_K1+1)               # hi ≥ 0
    )
    consVec <- c(-1+TOL, rep(0, n_K1+1))

    # find mles
    opt1 = stats::constrOptim(theta = rep(1/(n_K1+2), n_K1+1), loglik1_test0_h, control = list(fnscale=-1),
                              grad = grad1_test0_h, method = "BFGS", ui = consMat, ci=consVec,
                              y1=y1, K1_list=K1_list, K2=K2, i1=i1)
    s2hat1 <- loglik1_test0_h(par=opt1$par, y1=y1, K1_list=K1_list, K2=K2, i1=i1, return.s2phat = TRUE)

    opt01 = stats::optim(par = 0, loglik01_H0_test0_h, control = list(fnscale=-1), # par is h2_M
                         method = "L-BFGS-B", lower = TOL, upper = 1-TOL,
                         y0=y0, K2=K2, i0=i0)
    teststat <- ratio_test0_h(y0=y0, y1=y1, K1_list=K1_list, K2=K2, i1=i1, i0=i0, opt1=opt1, s2hat1=s2hat1, opt01=opt01)
    return(teststat)

  } else if (parameter.set == 'sigma2') { # num of pars is n_K+1
    opt1 = stats::optim(par = rep(1, n_K1+2), loglik1_test0_s, control = list(fnscale=-1),
                        y1=y1, K1_list=K1_list, K2=K2, i1=i1, method = "L-BFGS-B", lower = rep(TOL, n_K1+2))

    opt01 = stats::optim(par=rep(1, 2), fn = loglik01_H0_test0_s, control = list(fnscale=-1),
                         y0=y0, K2=K2, i0=i0, method = "L-BFGS-B", lower = rep(TOL, 2))
    teststat <- ratio_test0_s(y0, y1, K1_list, K2, i1, i0, opt1, opt01)
    return(teststat)
  }
  else {
    stop("parameter.set must be 'h2' or 'sigma2'.")
  }
}
