#' Calculate Log-Likelihood (Orthogonal K, Sigma2 Parameterization)
#'
#' Computes the log-likelihood for the orthogonal K matrices case with direct
#' sigma2 parameterization (s21, s22, ..., s2e). K_list contains diagonal values.
#'
#' @param par A numeric vector of variance component parameters (s21, s22, ..., s2e).
#'   The last element is assumed to be the residual error variance (s2e).
#' @param y Numeric vector of observed data.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix, corresponding to `par[i]`.
#' @return The log-likelihood value.
loglik_ortho_s <- function(par, y, K_list) {
  # Calculate the denominator: sum(par_i * K_i) + s2e
  deno <- Reduce(`+`, lapply(1:(length(par)-1), function(j) {par[j] * K_list[[j]]}))  + par[length(par)]
  l <- -(sum(log(deno) + y^2 / deno)) / 2
  return(l)
}

#' Calculate Log-Likelihood (Orthogonal K, H2/S2 Parameterization)
#'
#' Computes the log-likelihood for the orthogonal K matrices case with h2/s2 parameterization.
#' `par` is ordered as (h12, h22, ..., s2).
#'
#' @param par A numeric vector of parameters. `par[1]` to `par[length(par)-1]`
#'   are `h2_i` (proportions of variance for `K_list[[i]]`). The last parameter,
#'   `par[length(par)]`, is `s2` (total variance).
#' @param y Numeric vector of observed data.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @return The log-likelihood value.
loglik_ortho_h <- function(par, y, K_list) {
  n_K <- length(K_list)

  deno <- par[length(par)] * (Reduce(`+`, lapply(1:n_K, function(j) {par[j] * K_list[[j]]})) +
                                1-sum(par[1:n_K]))
  # # change parameter set to s2i and s2e
  # s2_i_list <- lapply(par[1:n_K], function(h2) h2 * par[length(par)])
  # s2e <- par[length(par)] * (1 - sum(par[1:n_K]))
  #
  # deno <- Reduce(`+`, lapply(1:n_K, function(j) {s2_i_list[[j]] * K_list[[j]]}))  + s2e
  l <- -(sum(log(deno) + y^2 / deno)) / 2
  return(l)
}


#' Calculate Gradient (Orthogonal K, Sigma2 Parameterization, for s2_i)
#'
#' Computes the gradient of the log-likelihood with respect to a specific `s2_i`
#' parameter in the orthogonal K matrices case.
#'
#' @param s2i_param_index The index of the specific `s2_i` in the `par` vector.
#' @param par The full parameter vector (s21, s22, ..., s2e).
#' @param y Numeric vector of observed data.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @param TOL Numerical tolerance to be considered as zero (default: 1e-8).
#' @return The gradient value for `s2_i`.
grad_ortho_s2i <- function(s2i_param_index, par, y, K_list, TOL = 1e-8) {
  s2i <- par[s2i_param_index]
  s2e <- par[length(par)]
  ki <- K_list[[s2i_param_index]]
  # Find non-zero entries in Ki to optimize calculations
  nonzero <- which(abs(ki) > TOL)
  ki <- ki[nonzero]
  y <- y[nonzero]

  # Calculate denominator: s2i_val * Ki_mat + sum(other_s2_j * K_j) + s2e_val
  deno <- s2i * ki + s2e
  g <- sum(ki / deno * (1 - y^2 / deno)) / 2
  return(g)
}

#' Calculate Gradient (Orthogonal K, Sigma2 Parameterization, for s2e)
#'
#' Computes the gradient of the log-likelihood with respect to `s2e` (residual error variance)
#' in the orthogonal K matrices case.
#'
#' @param par The full parameter vector (s21, s22, ..., s2e).
#' @param y Numeric vector of observed data.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @return The gradient value for `s2e`.
grad_ortho_s2e <- function(par, y, K_list) {
  # Calculate denominator: sum(par_s2_i * K_i) + s2e_val
  deno <- Reduce(`+`, lapply(1:(length(par)-1), function(j) {par[j] * K_list[[j]]}))  + par[length(par)]
  g <- sum(1 / deno * (1 - y^2 / deno)) / 2
  return(g)
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

  deno <- par[length(par)] * (Reduce(`+`, lapply(1:n_K, function(j) {par[j] * K_list[[j]]})) +
                                1-sum(par[1:n_K]))

  # Derivative of D with respect to h2_i:
  # dD/dh2_i = s2 * K_i - s2
  # Note: K_list[[h2i_param_index]] corresponds to the specific h2_i at that index.
  dD_dh2i <- par[length(par)] * (K_list[[h2i_param_index]] - 1)

  # Gradient of L wrt h2_i
  g <- sum( (1 / deno) * (1 - y^2 / deno) * dD_dh2i ) / 2
  return(g)
}




#' Coordinate Descent Optimization (Orthogonal K)
#'
#' This function performs coordinate descent optimization to find maximum likelihood
#' estimates for parameters in the orthogonal K matrices setting. It supports both
#' sigma2 and h2 parameterizations.
#'
#' @param y Numeric vector of observed data.
#' @param K_list A list of numeric vectors, where each `K_list[[i]]` represents the
#'   diagonal values of a transformed covariance matrix.
#' @param grad_list A list of gradient functions. For 'sigma2' it's `list(grad_ortho_s2i, grad_ortho_s2e)`.
#'   For 'h2', it's `list(grad_ortho_h2i)` (as s2 has a closed form).
#' @param init Initial guess for parameters. If `NULL`, initializes based on `parameter.set`.
#' @param lr Learning rate for the optimization.
#' @param TOL Numerical tolerance for convergence and constraints (default: 1e-8).
#' @param max_iter Maximum number of iterations for the coordinate descent.
#' @param fix_s2_indices A numeric vector of indices of parameters to fix during optimization.
#'   (e.g., `c(1)` to fix the first parameter).
#' @param fixed_s2_values A numeric vector of values for the fixed parameters,
#'   corresponding to `fix_s2_indices`.
#' @param parameter.set A string indicating the parameterization:
#'   'h2' for proportions of variance (h2_i) and total variance (s2), or
#'   'sigma2' for direct variance components (s21, s22, ..., s2e).
#' @return A numeric vector of the estimated parameters.
coordinate_descent_ortho <- function(y, K_list, grad_list, init = NULL, lr = 0.1, TOL = 1e-8,
                                     max_iter = 1000, fix_s2_indices = NULL, fixed_s2_values = NULL, parameter.set) {
  n_K <- length(K_list)
  n_par <- n_K + 1 # s21, s22, ..., s2k, s2e

  # Initialize parameters
  if (is.null(init)) {
    if (parameter.set == 'sigma2') {
      params <- rep(1, n_K+1)
    } else {
      params <- c(rep(0.1, n_K), 1)
    }
  } else { params <- init }

  # Ensure fixed_s2_indices and fixed_s2_values are consistent
  if (!is.null(fix_s2_indices)) {
    if (length(fix_s2_indices) != length(fixed_s2_values)) {
      stop("Length of fix_s2_indices must match length of fixed_s2_values.")
    }
    params[fix_s2_indices] <- fixed_s2_values
  }
  grad_sq_sum <- rep(0, n_par)

  for (iter in 1:max_iter) {
    params_old <- params
    gvec <- numeric(n_par)

    # Calculate gradients for s2_i parameters
    for (i in 1:n_K) {
      if (is.null(fix_s2_indices) || !(i %in% fix_s2_indices)) {
        gvec[i] <- grad_list[[1]](i, params, y, K_list, TOL=TOL)
      } else {
        gvec[i] <- 0 # Gradient is zero for fixed parameters
      }
    }
    # Calculate gradient for s2e
    if (parameter.set == 'h2') {
      gvec[n_par] <- 0 # mle of s2 has closed form
    } else {
      gvec[n_par] <- grad_list[[2]](params, y, K_list)
    }

    # update new gradient sum of square
    grad_sq_sum <- grad_sq_sum + gvec^2
    # Compute adaptive learning rate for each parameter
    adaptive_lr <- lr / (sqrt(grad_sq_sum) + TOL)

    # Update each variable using its gradient
    for (i in 1:n_par) {
      if (is.null(fix_s2_indices) || !(i %in% fix_s2_indices)) {
        params[i] <- params[i] - adaptive_lr[i] * gvec[i]
        # Apply non-negativity constraint
        params[i] <- max(params[i], TOL)
      }
    }

    # apply sum of h2 less than 1 constraint by scaling params proportionally so their sum is less than 1
    if (parameter.set == 'h2') {
      # sum of h2
      sum_h2 <- sum(params[-length(params)])
      if (sum_h2 > 1-TOL) {
        # Scale s2_params proportionally so their sum is 1
        scale_factor <- (1-TOL) / sum_h2
        for (i in 1:n_K) {
          if (is.null(fix_s2_indices) || !(i %in% fix_s2_indices)) {
            params[i] <- params[i] * scale_factor
            # Ensure non-negativity after scaling
            params[i] <- max(params[i], TOL)
          }
        }
      }
      # find estimator of s2
      deno <- Reduce(`+`, lapply(1:n_K, function(j) {params[j] * K_list[[j]]})) + 1-sum(params[1:n_K])
      params[length(params)] <- sum(y^2 / deno) / length(y)
    }
    # Convergence check: Stop if changes are small
    if (max(abs(params - params_old)) < TOL) {
      break
    }
  }
  if (iter == max_iter) { # reach maximum iteration
    print("DIDN'T FOUND")
  }
  return(params)
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
  Y0 <- y[i0]
  Y1 <- y[i1]

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
    est1 <- coordinate_descent_ortho(y = Y1, K_list = c(K1_list_subset_1, K2_list_subset_1),
                                     grad_list = list(grad_ortho_s2i, grad_ortho_s2e), parameter.set=parameter.set, TOL=TOL)
    # print(est1)
    est0 <- coordinate_descent_ortho(y = Y0, K_list = c(K1_list_subset_0, K2_list_subset_0),
                                     grad_list = list(grad_ortho_s2i, grad_ortho_s2e),
                                     fix_s2_indices = 1:length(K1_list), fixed_s2_values = rho, parameter.set=parameter.set, TOL=TOL)
    #print(est0)
  } else if (parameter.set == 'h2') {
    # Find MLEs
    est1 <- coordinate_descent_ortho(y = Y1, K_list = c(K1_list_subset_1, K2_list_subset_1),
                                     grad_list = list(grad_ortho_h2i), parameter.set=parameter.set, TOL=TOL)
    est0 <- coordinate_descent_ortho(y = Y0, K_list = c(K1_list_subset_0, K2_list_subset_0),
                                     grad_list = list(grad_ortho_h2i),
                                     fix_s2_indices = 1:length(K1_list), fixed_s2_values = rho, parameter.set=parameter.set, TOL=TOL)
  } else {
    stop("parameter.set must be 'h2' or 'sigma2'.")
  }

  # Calculate log-likelihood for the alternative hypothesis using Y0 and K_list[index0]
  l_theta1 <- loglik_ortho_s(par = est1, y = Y0, K_list = c(K1_list_subset_0, K2_list_subset_0))
  # Calculate log-likelihood for the null hypothesis using Y0 and K_list[index0]
  l_rho <- loglik_ortho_s(par = est0, y = Y0, K_list = c(K1_list_subset_0, K2_list_subset_0))

  return(exp(l_theta1 - l_rho))
}
