#' Estimation of network structure and variable selection in the nonlinear model with measurement errors in responses and covariates.
#'
#' This function characterizes Y and X by nonlinear models and is designed for detecting network structure and variable selection with measurement error in responses and covariates. Here the components of Y can be continuous, binary, or count. The estimation strategy in this function includes the regression calibration for correcting error-prone responses and covariates, the random forest method for marginally characterizing the response and covariates, and the distance correlation and graphical lasso for detecting the network structure among the responses.
#'
#' @param W A n × m response matrix. The variables can be error-prone or precisely measured, and can include continuous, binary, or count random variables.
#' @param Z A n × p matrix of continuous covariates. The variables can be error-prone or precisely measured.
#' @param sigma_eta A p × p covariance matrix of the noise term \eqn{\eta} in the classical measurement error model Z = X + \eqn{\eta}, where X is the unobserved version of Z.
#' @param rho A tuning parameter for the graphical lasso.
#' @param sigma_delta The common value in the diagonal covariance matrix of the noise term \eqn{\delta} in the classical measurement error model for continuous components in W. The default value is 0.5.
#' @param r A probability r for misclassification when components in W are binary. The default value is 0.8.
#' @param lambda A parameter \eqn{\lambda} in the Poisson distribution that provides the increasing measurement error effects when components in W are count. The default value is 1.
#' @param pi A parameter \eqn{\pi} in [0,1] for the Binomial distribution that characterizes the decreasing measurement error effects when components in W are count. The default value is 0.8.
#' @param label_name The name of the response variable. The default value is TRUE, which reflects the labels from the input data. Else, users can input the required labels manually.
#' @param var_thred A positive value used to retain important covariates. That is, covariates will be selected when refitting the model if their importance scores are greater than var_thred. The default value is 5.
#'
#' @return
#'   \item{W_hat}{The n × m matrix of corrected responses determined by regression calibration.}
#'   \item{Z_hat}{The n × p matrix of corrected covariates determined by regression calibration..}
#'   \item{PSE}{The Frobenius norm of the residual corresponding to W_hat.}
#'   \item{importance_score}{A matrix containing importance scores for the covariates.}
#'   \item{precision_matrix}{An estimated matrix reflecting the network structure of the responses.}
#'   \item{graph}{An visualization of the estimated network structure by \code{precision_matrix}.}
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
#'
#' NP_Graph(W, Z, sigma_eta, rho = 0.2,
#'                    sigma_delta = 0.5, r = 0.8,
#'                    lambda = 1, pi = 0.8,
#'                    label_name = TRUE, var_thred = 3)
#'
#' @export
NP_Graph = function(W, Z, sigma_eta, rho, sigma_delta = 0.5, r = 0.8, lambda = 1, pi = 0.8, label_name, var_thred  = 5) {
  #correcting W
  correct_W = function(W, sigma_delta, r, lambda, pi) {
    m = ncol(W)
    n = nrow(W)

    model = c()
    for (j in 1:m) {
      if (sum(W[, j] %% 1 == 0) == n & sum(W[, j] != 0 & W[, j] != 1) > 0) {
        model = c(model, "counts")
      } else if (sum(W[, j] %% 1 == 0) == n & sum(W[, j] != 0 & W[, j] != 1) == 0) {
        model = c(model, "binary")
      } else {
        model = c(model, "continuous")
      }
    }

    Y = data.frame(matrix(ncol = m, nrow = n))
    colnames(Y) = colnames(W)
    for (j in 1:m) {
      if (model[j] == "continuous") {
        Y[, j] = mean(W[, j]) + (stats::var(W[, j]) - sigma_delta) / stats::var(W[, j]) * (W[, j] - mean(W[, j]))
      } else if (model[j] == "binary") {
        S = stats::rbinom(n, 1, r)
        Y[, j] = (W[, j] + S - 1) / (2 * S - 1)
      } else {
        Y[, j] = ((mean(W[, j]) - lambda) / (1 - pi)) +
          ((mean(W[, j]) - lambda) / (1 - pi)) *
          (mean(W[, j]) - lambda) / (lambda + (3 * pi + 1) / (1 - pi) * (mean(W[, j]) - lambda)) * (W[, j] - mean(W[, j]))
      }
    }
    return(Y)
  }
  #correcting Z
  correct_Z = function(Z, sigma_eta) {
    covariance_matrix_Z = as.matrix(stats::cor(Z))
    n = nrow(Z)
    p = ncol(Z)
    X = matrix(0, nrow = n, ncol = p)
    for (i in 1:n) {
      X[i, ] = colMeans(Z) + t(covariance_matrix_Z - sigma_eta) %*% solve(covariance_matrix_Z + diag(0.2, dim(Z)[2])) %*% (Z[i, ] - colMeans(Z))
    }
    return(as.data.frame(X))
  }

  Y = as.matrix(correct_W(W, sigma_delta, r, lambda, pi))
  X = as.matrix(correct_Z(Z, sigma_eta))
  #random forest and importance score
  models = list()
  importance = list()
  for (i in 1:ncol(Y)) {
    df_Y = Y[, i]
    df = cbind(X, Y = df_Y)
    model = randomForest::randomForest(Y ~ ., data = df, importance = TRUE)
    models[[i]] = model
    importance[[i]] = caret::varImp(model)
  }
  importance_matrix = do.call(cbind, importance)
  #refit
  rf_model_refit = list()
  for (i in 1:ncol(Y)) {
    df_Y = Y[, i]
    df = cbind(X, Y = df_Y)
    imp_sel = importance[[i]][, "Overall"] > var_thred
    important_vars = rownames(importance[[i]][imp_sel, , drop = FALSE])

    if (length(important_vars) > 0) {
      rf_model_refit[[i]] = randomForest::randomForest(Y ~ ., data = df[, c("Y", important_vars)])
    } else {
      max_var = rownames(importance[[i]])[which.max(importance[[i]][, "Overall"])]
      rf_model_refit[[i]] = randomForest::randomForest(Y ~ ., data = df[, c("Y", max_var)])
    }
  }

  #residual
  r_list = list()
  for (i in 1:ncol(Y)) {
    r_i = Y[, i] - stats::predict(rf_model_refit[[i]], X)
    r_list[[i]] = r_i
  }
  r_matrix = do.call(cbind, r_list)
  PSE = mean(sapply(r_list, function(r_matrix) sum(r_matrix^2)))
  y_matrix = as.matrix(Y)

  #distance_correlation
  n_resp = ncol(y_matrix)
  correlation = matrix(1, nrow = n_resp, ncol = n_resp)
  for (i in 1:n_resp) {
    for (j in 1:n_resp) {
      if (i != j) {
        correlation[i, j] = metrica::dcorr(obs = y_matrix[, i], pred = r_matrix[, j])
      }
    }
  }
  for (i in 1:(n_resp - 1)) {
    for (j in (i + 1):n_resp) {
      max_value = max(correlation[i, j], correlation[j, i])
      correlation[i, j] = max_value
      correlation[j, i] = max_value
    }
  }
  #graphical_lasso
  glasso_result = glasso::glasso(correlation, rho = rho)
  precision_matrix = glasso_result$wi

  #graph the network structure
  net = precision_matrix
  net = network::network(net, directed = FALSE)
  network::network.vertex.names(net)=paste0("Y",network::network.vertex.names(net))
  graph = GGally::ggnet2(net,size=10,node.color = "lightgray",label=label_name,label.size = 3,mode = "circle")

  return(list(
    W_hat = Y,
    Z_hat = X,
    PSE = PSE,
    importance_matrix = importance_matrix,
    precision_matrix = precision_matrix,
    graph = graph
  ))
}
