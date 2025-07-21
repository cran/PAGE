#' Estimation of network structure and variable selection in the linear model via the Gaussian maximum likelihood.
#'
#' This function focuses on multivariate linear regression models Y = XB + \eqn{\epsilon} subject to measurement error in the responses and covariates, where with B is a matrix of parameters and \eqn{\epsilon} is a noise term with zero expectation. We aim to detect the network structure of responses and select informative covaraites. The estimation method is motivated by the Gaussian maximum likelihood function and uses the conditional expectation to correct for measurement error.
#'
#' @param W A n × m response matrix, the variables can be error-prone or precisely measured.
#' @param Z A n × p matrix of continuous covariates. The variables can be error-prone or precisely measured.
#' @param sigma_eta A p × p covariance matrix of the noise term \eqn{\eta} in the classical measurement error model Z = X + \eqn{\eta}, where X is the unobserved version of Z.
#' @param sigma_delta A m × m covariance matrix of the noise term \eqn{\delta} in the classical measurement error model W = Y + \eqn{\delta}, where Y is the unobserved version of W.
#' @param alpha_1 A tuning parameter associated with the parameter B.
#' @param alpha_2 A tuning parameter associated with precision matrix C, which is the inverse of the covariance matrix of \eqn{\epsilon}.
#' @param alpha_1_list A list of tuning parameters for the model averaging estimator of B. The default value is NULL.
#' @param alpha_2_list A list of tuning parameters for the model averaging estimator of C. The default value is NULL.
#' @param label_name The name of the response variable. The default value is TRUE, which reflects the labels from the input data. Else, users can input the required labels manually.
#'
#' @return
#'   \item{B}{An estimator of B.}
#'   \item{C}{An estimator of C.}
#'   \item{graph}{A visualization of the estimated network structure by C.}
#'   \item{Beta_BICs}{A vector of Bayesian Information Criterion (BIC) weights for the model averaging estimator of B under candidate models alpha_1_list.}
#'   \item{Gamma_BICs}{A vector of Bayesian Information Criterion (BIC) weights for the model averaging estimator of C under candidate models alpha_2_list.}
#'
#'
#' @author
#' Wan-Yi Chang and Li-Pang Chen \cr
#' Maintainer: Wan-Yi Chang \email{jessica306a@gmail.com}
#'
#' @examples
#' n <- 100
#' Z <- matrix(rnorm(n * 10), n, 10)
#' W <- matrix(rnorm(n * 5), n, 5)
#' sigma_eta <- diag(0.15, ncol(Z))
#' sigma_delta <- diag(0.3, ncol(W))
#'
#'Joint_Gaussian(W, Z, sigma_eta, sigma_delta,
#'                          alpha_1 = 0.1, alpha_2 = 0.1,
#'                          alpha_1_list = c(0.1, 0.3),
#'                          alpha_2_list = c(0.1, 0.3),
#'                          label_name = TRUE)
#'
#' @useDynLib PAGE, .registration = TRUE
#' @export
Joint_Gaussian = function(W, Z, sigma_eta, sigma_delta, alpha_1, alpha_2,
                          alpha_1_list = NULL, alpha_2_list = NULL,label_name = TRUE) {
  #compute every B and C
  compute_Beta_gamma <- function(alpha_1_val, alpha_2_val) {
    Estimates <- JointSol(W, Z, alpha_1 = alpha_1_val, alpha_2 = alpha_2_val)
    Estimates_beta <- Estimates$B

    Beta <- MASS::ginv(stats::var(Z) - sigma_eta) %*% stats::var(Z) %*% Estimates_beta
    A <- t(W - Z %*% Beta) %*% (W - Z %*% Beta) / nrow(W) -
      (t(Beta) %*% sigma_eta %*% Beta) - sigma_delta + diag(ncol(sigma_delta),ncol(sigma_delta),ncol(sigma_delta))

    # covariance matrix is positive definite
    if (all(eigen(A)$values > 0)) {
      gamma <- glasso::glasso(A, rho = alpha_2_val)$wi
      return(list(Beta = Beta, gamma = gamma))
    } else {
      return(NULL)
    }
  }

  # initial estimator
  base_result <- compute_Beta_gamma(alpha_1, alpha_2)
  # covariance matrix is not positive definite
  if (is.null(base_result)) stop("Base model (alpha_1, alpha_2) is not valid.")

  final_Beta <- base_result$Beta
  final_gamma <- base_result$gamma
  Beta_BICs <- compute_BIC_beta_joint(W, Z, final_Beta, sigma_eta, sigma_delta)
  Gamma_BICs <- compute_BIC_gamma_joint(W, Z, final_gamma, final_Beta, sigma_eta, sigma_delta)


  # average B
  if (!is.null(alpha_1_list)) {
    beta_list <- list()
    BICs <- c()

    for (lam in alpha_1_list) {
      res <- compute_Beta_gamma(lam, alpha_2)
      if (!is.null(res)) {
        beta_list[[length(beta_list) + 1]] <- res$Beta
        BICs <- c(BICs, compute_BIC_beta_joint(W, Z, res$Beta, sigma_eta, sigma_delta))
      }
    }

    if (length(BICs) > 0) {
      weights <- BICs / sum(BICs)
      final_Beta <- Reduce(`+`, Map(function(w, b) w * b, weights, beta_list))
      Beta_BICs <- BICs
    }
  }

  # average C
  if (!is.null(alpha_2_list)) {
    gamma_list <- list()
    BICs <- c()

    for (r in alpha_2_list) {
      res <- compute_Beta_gamma(alpha_1, r)
      if (!is.null(res)) {
        gamma_list[[length(gamma_list) + 1]] <- res$gamma
        BICs <- c(BICs, compute_BIC_gamma_joint(W, Z, res$gamma, final_Beta, sigma_eta,sigma_delta))
      }
    }

    if (length(BICs) > 0) {
      weights <- BICs / sum(BICs)
      final_gamma <- Reduce(`+`, Map(function(w, g) w * g, weights, gamma_list))
      Gamma_BICs <- BICs
    }
  }

  #graph the network structure
  net = final_gamma
  net = network::network(net, directed = FALSE)
  network::network.vertex.names(net)=paste0("X",network::network.vertex.names(net))
  graph = GGally::ggnet2(net,size=10,node.color = "lightgray",label=label_name,label.size = 3,mode = "circle")

  return(list(B = final_Beta, C = final_gamma,graph = graph,Beta_BICs = Beta_BICs, Gamma_BICs = Gamma_BICs))
}

JointSol = function(W, Z, alpha_1, alpha_2){

  N = dim(W)[1]
  M = dim(W)[2]
  P = dim(Z)[2]

  ## Centered Data
  cent.W = W - matrix(rep(colMeans(W), each=N), N, M)
  cent.Z = Z - matrix(rep(colMeans(Z), each=N), N, P)

  ## Estimating B and C
  update_BC = function(alpha_1, alpha_2){

    ## Initial B
    init.B = matrix(0, P, M)

    for(i in 1:M){
      lasso = lars::lars(cent.Z, cent.W[, i], use.Gram=FALSE, intercept=FALSE)
      init.B[, i] = lars::predict.lars(lasso, cent.Z, alpha_1/2, type="coefficients", mode="lambda")$coefficients
    }

    ## Initial C
    init.A = t(cent.W - cent.Z %*% init.B) %*% (cent.W - cent.Z %*% init.B) / N
    init.C = glasso::glasso(init.A, rho = alpha_2, penalize.diagonal=FALSE)$wi
    mat = matrix(1, P, M)

    ## Updating_BC
    output = JointLasso(cent.W, cent.Z, init.B, init.C, mat, alpha_1, alpha_2)

    return(output)
  }

  result = update_BC(alpha_1, alpha_2)
  resultB = result[["B"]]
  resultC = result[["C"]]

  return(list(B=resultB, C=resultC))
}

JointLasso = function(cent.W, cent.Z, init.B, init.C, V, alpha_1, alpha_2){

  cur.B = init.B
  cur.C = init.C

  update = Onetime_update(cent.W, cent.Z, cur.B, cur.C, V, alpha_1, alpha_2)
  updateB = update[["updateB"]]
  updateC = update[["updateC"]]

  maxDiff = max(max(abs(cur.B - updateB)), max(abs(cur.C - updateC)))
  iter_n = 1

  while(maxDiff > 1.0e-3 && iter_n < 10)
  {
    cur.B = updateB
    cur.C = updateC

    update = Onetime_update(cent.W, cent.Z, cur.B, cur.C, V, alpha_1, alpha_2)
    updateB = update[["updateB"]]
    updateC = update[["updateC"]]

    maxDiff = max(max(abs(cur.B - updateB)), max(abs(cur.C - updateC)))
    iter_n = iter_n + 1
  }

  return(list(C = updateC, B = updateB, iteration = iter_n))
}

Onetime_update = function(cent.W, cent.Z, cur.B, cur.C, V, alpha_1, alpha_2){

  N = dim(cent.W)[1]
  M = dim(cent.W)[2]
  P = dim(cent.Z)[2]
  dims = c(N, M, P)

  ## update C
  A = t(cent.W - cent.Z %*% cur.B) %*% (cent.W - cent.Z %*% cur.B) / N
  Cest = glasso::glasso(A, rho = alpha_2, penalize.diagonal=FALSE, wi.init=cur.C)
  updateC = Cest$wi

  ## Update B
  updateB = numeric(P * M)
  iteration = numeric(1)

  #dyn.load("C:/Users/user/Desktop/jessica/PAGE/R/UpdateBeta.dll")
  output = .C("UpdateBeta",
              as.double(t(cent.W)),
              as.double(t(cent.Z)),
              as.double(t(cur.B)),
              as.integer(dims),
              as.double(t(cur.C)),
              as.double(t(V)),
              as.double(alpha_1),
              updateB = as.double(updateB),
              iteration = as.integer(iteration), PACKAGE = "PAGE")
  #dyn.unload("C:/Users/user/Desktop/jessica/PAGE/R/UpdateBeta.dll")

  updateB = matrix(output[["updateB"]], P, M, byrow = TRUE)

  return(list(updateC = updateC, updateB = updateB))
}

#compute BIC beta for penalty
compute_BIC_beta_joint = function(W, Z, B, sigma_eta,sigma_delta) {
  n = nrow(W)
  likelihood = sum(diag(t(W - Z %*% B)%*% (W - Z %*% B) - sigma_delta - t(B)%*%sigma_eta%*%B))
  nonzero_B =  sum(B != 0)
  BIC = -2 * likelihood + log(n) * nonzero_B
  return(BIC)
}

#compute BIC gamma for penalty
compute_BIC_gamma_joint = function(W, Z, C, B, sigma_eta,sigma_delta) {
  n = nrow(W)
  trace_A = sum(diag(t(W - Z %*% B)%*% (W - Z %*% B) - sigma_delta - t(B)%*%sigma_eta%*%B))
  likelihood = n*log(det(C)) - trace_A
  nonzero_C =  sum(C[upper.tri(C)])
  BIC = -2 * likelihood + log(n) * nonzero_C
  return(BIC)
}
