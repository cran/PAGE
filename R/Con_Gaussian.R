#' Estimation of network structure and variable selection in the linear model via the conditional likelihood function.
#'
#' This function focuses on multivariate linear regression models Y = XB + \eqn{\epsilon} subject to measurement error in responses and covariates, where with B is a matrix of parameters and \eqn{\epsilon} is a noise term with zero expectation. We aim to detect the network structure of responses and select informative covaraites. The estimation method is motivated by the conditional likelihood function and uses the conditional expectation to correct for measurement error.
#'
#' @param W A n × m response matrix, the variables can be error-prone or precisely measured.
#' @param Z A n × p matrix of continuous covariates. The variables can be error-prone or precisely measured.
#' @param sigma_eta A p × p covariance matrix of the noise term \eqn{\eta} in the classical measurement error model Z = X + \eqn{\eta}, where X is the unobserved version of Z.
#' @param sigma_delta A m × m covariance matrix of the noise term \eqn{\delta} in the classical measurement error model W = Y + \eqn{\delta}, where Y is the unobserved version of W.
#' @param alpha_1 A tuning parameter associated with parameter B.
#' @param alpha_2 A tuning parameter associated with parameter, denoted as \eqn{\Gamma}, that reflects the network in Y.
#' @param alpha_1_list A list of tuning parameters for the model averaging estimator of B. The default value is NULL.
#' @param alpha_2_list A list of tuning parameters for the model averaging estimator of \eqn{\Gamma}. The default value is NULL.
#' @param max_iter A maximum number for iterations for updated values of B and \eqn{\Gamma}. The default value is 30.
#' @param tol A prespecified tolerance \eqn{\zeta} for iterations for updated values of B and \eqn{\Gamma}. The default value is \eqn{10^{-6}}.
#' @param label_name The name of the response variables. The default value is TRUE, which reflects the labels from the input data. Else, users can input the required labels manually.
#'
#' @return
#'   \item{Beta}{An estimator of B.}
#'   \item{gamma}{An estimator of the network in Y.}
#'   \item{graph}{A visualization of the estimated network structure by gamma.}
#'   \item{Beta_BICs}{A vector of Bayesian Information Criterion (BIC) weights for the model averaging estimator of B under candidate models alpha_1_list.}
#'   \item{Gamma_BICs}{A vector of Bayesian Information Criterion (BIC) weights for the model averaging estimator of \eqn{\Gamma} under candidate models alpha_2_list.}
#'
#'
#' @author
#' Wan-Yi Chang and Li-Pang Chen \cr
#' Maintainer: Wan-Yi Chang \email{jessica306a@gmail.com}
#'
#' @examples
#' n <- 100
#' Z <- matrix(rnorm(n * 5), n, 5)
#' W <- matrix(rnorm(n * 5), n, 5)
#' sigma_eta <- diag(0.15, ncol(Z))
#' sigma_delta <- diag(0.3, ncol(W))
#'
#' Cond_Gaussian(W, Z, sigma_eta, sigma_delta,
#'                         alpha_1 = 0.1, alpha_2 = 0.1,
#'                         alpha_1_list = NULL,
#'                         alpha_2_list = NULL,
#'                         max_iter = 1, tol = 1e-6, label_name = TRUE)
#'
#'
#'
#' @export
Cond_Gaussian = function(W, Z, sigma_eta, sigma_delta, alpha_1, alpha_2,
                         alpha_1_list = NULL, alpha_2_list = NULL,
                         max_iter = 30, tol = 1e-6,label_name = TRUE) {

  #optimize Beta and Gamma
  compute_Beta_gamma <- function(alpha_1_val, alpha_2_val) {
    gamma <- gamma_transform_matrix(diag(ncol(W)))
    beta <- optimize_beta(Z, W, alpha_1_val, sigma_delta, sigma_eta)
    loss_history <- c()

    for (i in 1:max_iter) {
      gamma <- optimize_gamma(Z, W, beta, alpha_2_val, sigma_delta, sigma_eta)
      beta <- optimize_beta_new(Z, W, gamma, alpha_1_val, sigma_delta, sigma_eta)

      total_loss <- objective_function_1(as.vector(beta), Z, W, alpha_1_val, sigma_delta, sigma_eta)
      loss_history <- c(loss_history, total_loss)
    }

    return(list(beta = beta, gamma = gamma, loss_history = loss_history))
  }

  #computation with single alpha_1, alpha_2
  base_result <- compute_Beta_gamma(alpha_1, alpha_2)
  final_Beta <- base_result$beta
  final_gamma <- base_result$gamma
  Gamma_BICs <- compute_BIC_gamma(W, Z, final_gamma, final_Beta, sigma_delta, sigma_eta)
  Beta_BICs <- compute_BIC_beta(W, Z, final_gamma, final_Beta, sigma_delta, sigma_eta)

  # MA average Beta
  if (!is.null(alpha_1_list)) {
    beta_list <- list()
    BICs <- c()

    for (lam in alpha_1_list) {
      res <- compute_Beta_gamma(lam, alpha_2)
      beta_list[[length(beta_list) + 1]] <- res$beta
      BICs <- c(BICs, compute_BIC_beta(W, Z, final_gamma, res$beta, sigma_delta, sigma_eta))
    }

    if (length(BICs) > 0) {
      weights <- BICs / sum(BICs)
      final_Beta <- Reduce(`+`, Map(function(w, b) w * b, weights, beta_list))
      Beta_BICs <- BICs
    }
  }

  # MA average Gamma
  if (!is.null(alpha_2_list)) {
    gamma_list <- list()
    BICs <- c()

    for (r in alpha_2_list) {
      res <- compute_Beta_gamma(alpha_1, r)
      gamma_list[[length(gamma_list) + 1]] <- res$gamma
      BICs <- c(BICs, compute_BIC_gamma(W, Z, res$gamma, final_Beta, sigma_delta, sigma_eta))
    }

    if (length(BICs) > 0) {
      weights <- BICs / sum(BICs)
      final_gamma <- Reduce(`+`, Map(function(w, g) w * g, weights, gamma_list))
      Gamma_BICs <- BICs
    }
  }

  # Force gamma to be symmetric
  final_gamma <- transform_gamma_symmetric(final_gamma)

  # Graph the network structure
  net <- final_gamma
  net <- network::network(net, directed = FALSE)
  network::network.vertex.names(net) <- paste0("X", network::network.vertex.names(net))
  graph <- GGally::ggnet2(net, size = 10, node.color = "lightgray", label = label_name, label.size = 3, mode = "circle")

  return(list(Beta = final_Beta, gamma = final_gamma, graph = graph,Beta_BICs = Beta_BICs, Gamma_BICs = Gamma_BICs))
}

# initial of beta in equation
objective_function_1 <- function(beta_vec, Z, W, alpha_1, sigma_delta, sigma_eta) {
  g <- ncol(W)
  beta <- matrix(beta_vec, nrow = ncol(Z), ncol = g)

  total_loss <- 0

  for (k in 1:g) {
    residual <- W[, k] - Z %*% beta[, k]
    total_loss <- total_loss + sum(residual^2)
  }

  var_term <- 0
  for (k in 1:g) {
    var_term <- var_term + nrow(Z) * (sum(diag(sigma_delta))
                                      + t(beta[, k]) %*% sigma_eta %*% beta[, k])
  }

  reg_term <- alpha_1 * sum(abs(beta))

  total_loss <- total_loss - var_term + reg_term
  return(total_loss)
}

optimize_beta <- function(Z, W, alpha_1, sigma_delta, sigma_eta) {
  beta_init <- rep(0, ncol(Z) * ncol(W))
  result <- stats::optim(par = beta_init,
                  fn = objective_function_1,
                  Z = Z, W = W, alpha_1 = alpha_1, sigma_delta = sigma_delta, sigma_eta = sigma_eta,
                  method = "BFGS")
  return(matrix(result$par, nrow = ncol(Z), ncol = ncol(W)))
}

# objective gamma
objective_function_2 <- function(gamma_vec, Z, W, beta, alpha_2, sigma_delta, sigma_eta) {
  g <- ncol(W)
  gamma <- matrix(gamma_vec, nrow = g - 1, ncol = g)
  total_loss <- 0

  for (k in 1:g) {
    residual <- (W[, k] - Z %*% beta[, k]) - ((W[, -k] - Z %*% beta[, -k]) %*% gamma[, k])
    total_loss <- total_loss + sum(residual^2)
  }

  total_var_term <- 0
  for (k in 1:g) {
    gamma_k <- as.matrix(gamma[, k])
    total_var_term <- total_var_term +
      nrow(W) * sigma_delta[k, k] +
      t(gamma_k) %*% sigma_delta[-k, -k] %*% gamma_k -
      2 * t(gamma_k) %*% (nrow(W) * sigma_delta[-k, k]) +
      nrow(W) * (t(beta[, k]) %*% sigma_eta %*%  beta[, k]) +
      t(gamma_k) %*% (t(beta[, -k]) %*% sigma_eta %*%  beta[, -k]) %*% gamma_k -
      2 * t(gamma_k) %*% ((nrow(W) * t(beta[, -k]) %*% sigma_eta %*%  beta[, k]))
  }

  reg_term <- alpha_2 * sum(abs(gamma))

  total_loss <- total_loss - total_var_term + reg_term
  return(total_loss)
}

optimize_gamma <- function(Z, W, beta, alpha_2, sigma_delta, sigma_eta) {
  gamma_init <- rep(0, (ncol(W) - 1) * ncol(W))
  result <- stats::optim(par = gamma_init,
                  fn = objective_function_2,
                  Z = Z, W = W, beta = beta, alpha_2 = alpha_2, sigma_delta = sigma_delta, sigma_eta = sigma_eta,
                  method = "BFGS")
  return(matrix(result$par, nrow = ncol(W) - 1, ncol = ncol(W)))
}

# optimize beta
objective_function_3 <- function(beta_vec, Z, W, gamma, alpha_1, sigma_delta, sigma_eta) {
  g <- ncol(W)
  beta <- matrix(beta_vec, nrow = ncol(Z), ncol = g)
  total_loss <- 0

  # 第一項: sum_k ||W^k - X * beta_k||_2^2
  for (k in 1:g) {
    residual <- (W[, k] - Z %*% beta[, k]) - ((W[, -k] - Z %*% beta[, -k]) %*% gamma[, k])
    total_loss <- total_loss + sum(residual^2)
  }

  # 第二項: g(k)
  total_var_term <- 0
  for (k in 1:g) {
    gamma_k <- as.matrix(gamma[, k])
    total_var_term <- total_var_term +
      nrow(W) * sigma_delta[k, k] +
      t(gamma_k) %*% sigma_delta[-k, -k] %*% gamma_k -
      2 * t(gamma_k) %*% (nrow(W) * sigma_delta[-k, k])+
      nrow(W) * (t(beta[, k]) %*% sigma_eta %*%  beta[, k]) +
      t(gamma_k) %*% (t(beta[, -k]) %*% sigma_eta %*%  beta[, -k]) %*% gamma_k -
      2 * t(gamma_k) %*% ((nrow(W) * t(beta[, -k]) %*% sigma_eta %*%  beta[, k]))
  }

  # 第三項: Lasso 正則化
  reg_term <- alpha_1 * sum(abs(beta))

  total_loss <- total_loss - total_var_term + reg_term
  return(total_loss)
}
optimize_beta_new <- function(Z, W, gamma, alpha_1, sigma_delta, sigma_eta) {
  beta_init <- rep(0, ncol(Z) * ncol(W))
  result <- stats::optim(par = beta_init,
                  fn = objective_function_3,
                  Z = Z, W = W, gamma = gamma, alpha_1 = alpha_1, sigma_delta = sigma_delta, sigma_eta = sigma_eta,
                  method = "BFGS")
  return(matrix(result$par, nrow = ncol(Z), ncol = ncol(W)))
}

#compute BIC of gamma in differnet penalty
compute_BIC_gamma <- function(W, Z, gamma, B, sigma_delta, sigma_eta) {
  n <- nrow(W)
  g <- ncol(W)
  if (g < 2) stop("Number of variables g must be >= 2")

  total_likelihood <- 0
  total_var_term <- 0
  total_nonzero <- sum(gamma != 0)

  for (k in 1:g) {
    gamma_k <- as.matrix(gamma[, k])
    if (g == 1) {
      residual <- (W[, k] - Z %*% B[, k])
    } else {
      residual <- (W[, k] - Z %*% B[, k]) - ((W[, -k, drop = FALSE] - Z %*% B[, -k, drop = FALSE]) %*% gamma_k)
    }
    total_likelihood <- total_likelihood + sum(residual^2)
  }

  for (k in 1:g) {
    gamma_k <- as.matrix(gamma[, k])

    if (g == 1) {
      var_term <- as.numeric(
        n * sigma_delta[k, k] + n * t(B[, k]) %*% sigma_eta %*% B[, k]
      )
    } else {
      B_minus_k <- B[, -k, drop = FALSE]
      sigma_d_minus_k <- sigma_delta[-k, -k, drop = FALSE]
      sigma_d_k <- sigma_delta[-k, k, drop = FALSE]

      var_term <- as.numeric(
        n * sigma_delta[k, k] +
          t(gamma_k) %*% sigma_d_minus_k %*% gamma_k -
          2 * t(gamma_k) %*% (n * sigma_d_k) +
          n * t(B[, k]) %*% sigma_eta %*% B[, k] +
          t(gamma_k) %*% t(B_minus_k) %*% sigma_eta %*% B_minus_k %*% gamma_k -
          2 * t(gamma_k) %*% (n * t(B_minus_k) %*% sigma_eta %*% B[, k])
      )
    }
    total_var_term <- total_var_term + var_term
  }

  BIC <- -2 * (total_likelihood + total_var_term) + log(n) * total_nonzero
  return(BIC)
}

compute_BIC_beta <- function(W, Z, gamma, B, sigma_delta, sigma_eta) {
  n <- nrow(W)
  g <- ncol(W)
  if (g < 2) stop("Number of variables g must be >= 2")

  total_likelihood <- 0
  total_var_term <- 0
  total_nonzero <- sum(B != 0)

  for (k in 1:g) {
    gamma_k <- as.matrix(gamma[, k])
    if (g == 1) {
      residual <- (W[, k] - Z %*% B[, k])
    } else {
      residual <- (W[, k] - Z %*% B[, k]) - ((W[, -k, drop = FALSE] - Z %*% B[, -k, drop = FALSE]) %*% gamma_k)
    }
    total_likelihood <- total_likelihood + sum(residual^2)
  }

  for (k in 1:g) {
    gamma_k <- as.matrix(gamma[, k])

    if (g == 1) {
      var_term <- as.numeric(
        n * sigma_delta[k, k] + n * t(B[, k]) %*% sigma_eta %*% B[, k]
      )
    } else {
      B_minus_k <- B[, -k, drop = FALSE]
      sigma_d_minus_k <- sigma_delta[-k, -k, drop = FALSE]
      sigma_d_k <- sigma_delta[-k, k, drop = FALSE]

      var_term <- as.numeric(
        n * sigma_delta[k, k] +
          t(gamma_k) %*% sigma_d_minus_k %*% gamma_k -
          2 * t(gamma_k) %*% (n * sigma_d_k) +
          n * t(B[, k]) %*% sigma_eta %*% B[, k] +
          t(gamma_k) %*% t(B_minus_k) %*% sigma_eta %*% B_minus_k %*% gamma_k -
          2 * t(gamma_k) %*% (n * t(B_minus_k) %*% sigma_eta %*% B[, k])
      )
    }
    total_var_term <- total_var_term + var_term
  }

  BIC <- -2 * (total_likelihood + total_var_term) + log(n) * total_nonzero
  return(BIC)
}

# symmetric matrix transform to condition gamma_k
gamma_transform_matrix = function(mat) {
  n = nrow(mat)
  new_mat = matrix(NA, n-1, n)
  for (i in 1:n) {
    diag_element = mat[i, i]
    new_mat[, i] = -mat[-i, i] / diag_element
  }
  return(new_mat)
}

#gamma_k transform symmetric matrix
transform_gamma_symmetric = function(original_matrix) {
  # Initialize the new matrix with NA
  new_matrix = matrix(NA, nrow = ncol(original_matrix), ncol = ncol(original_matrix))

  # Fill the new matrix based on the original matrix
  for (i in 1:nrow(original_matrix)) {
    new_matrix[(i+1):ncol(original_matrix), i] = original_matrix[i:nrow(original_matrix), i]
    new_matrix[1:(i-1), i] = original_matrix[1:(i-1), i]
    new_matrix[i, i] = 1
  }

  # Set the last element in the matrix to 1
  new_matrix[ncol(original_matrix), ncol(original_matrix)] = 1

  # Fill the last column of the new matrix with the last column of the original matrix
  new_matrix[1:nrow(original_matrix), ncol(original_matrix)] = original_matrix[, ncol(original_matrix)]

  return(new_matrix)
}



