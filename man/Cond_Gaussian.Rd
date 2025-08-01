% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Con_Gaussian.R
\name{Cond_Gaussian}
\alias{Cond_Gaussian}
\title{Estimation of network structure and variable selection in the linear model via the conditional likelihood function.}
\usage{
Cond_Gaussian(
  W,
  Z,
  sigma_eta,
  sigma_delta,
  alpha_1,
  alpha_2,
  alpha_1_list = NULL,
  alpha_2_list = NULL,
  max_iter = 30,
  tol = 1e-06,
  label_name = TRUE
)
}
\arguments{
\item{W}{A n × m response matrix, the variables can be error-prone or precisely measured.}

\item{Z}{A n × p matrix of continuous covariates. The variables can be error-prone or precisely measured.}

\item{sigma_eta}{A p × p covariance matrix of the noise term \eqn{\eta} in the classical measurement error model Z = X + \eqn{\eta}, where X is the unobserved version of Z.}

\item{sigma_delta}{A m × m covariance matrix of the noise term \eqn{\delta} in the classical measurement error model W = Y + \eqn{\delta}, where Y is the unobserved version of W.}

\item{alpha_1}{A tuning parameter associated with parameter B.}

\item{alpha_2}{A tuning parameter associated with parameter, denoted as \eqn{\Gamma}, that reflects the network in Y.}

\item{alpha_1_list}{A list of tuning parameters for the model averaging estimator of B. The default value is NULL.}

\item{alpha_2_list}{A list of tuning parameters for the model averaging estimator of \eqn{\Gamma}. The default value is NULL.}

\item{max_iter}{A maximum number for iterations for updated values of B and \eqn{\Gamma}. The default value is 30.}

\item{tol}{A prespecified tolerance \eqn{\zeta} for iterations for updated values of B and \eqn{\Gamma}. The default value is \eqn{10^{-6}}.}

\item{label_name}{The name of the response variables. The default value is TRUE, which reflects the labels from the input data. Else, users can input the required labels manually.}
}
\value{
\item{B}{An estimator of B.}
  \item{gamma}{An estimator of the network in Y.}
  \item{graph}{A visualization of the estimated network structure by gamma.}
  \item{Beta_BICs}{A vector of Bayesian Information Criterion (BIC) weights for the model averaging estimator of B under candidate models alpha_1_list.}
  \item{Gamma_BICs}{A vector of Bayesian Information Criterion (BIC) weights for the model averaging estimator of \eqn{\Gamma} under candidate models alpha_2_list.}
}
\description{
This function focuses on multivariate linear regression models Y = XB + \eqn{\epsilon} subject to measurement error in responses and covariates, where with B is a matrix of parameters and \eqn{\epsilon} is a noise term with zero expectation. We aim to detect the network structure of responses and select informative covaraites. The estimation method is motivated by the conditional likelihood function and uses the conditional expectation to correct for measurement error.
}
\examples{
n <- 100
Z <- matrix(rnorm(n * 5), n, 5)
W <- matrix(rnorm(n * 5), n, 5)
sigma_eta <- diag(0.15, ncol(Z))
sigma_delta <- diag(0.3, ncol(W))

Cond_Gaussian(W, Z, sigma_eta, sigma_delta,
                        alpha_1 = 0.1, alpha_2 = 0.1,
                        alpha_1_list = NULL,
                        alpha_2_list = NULL,
                        max_iter = 1, tol = 1e-6, label_name = TRUE)



}
\author{
Wan-Yi Chang and Li-Pang Chen \cr
Maintainer: Wan-Yi Chang \email{jessica306a@gmail.com}
}
