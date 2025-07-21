#' @title Predictor-Assisted Graphical Models under Error-in-Variables
#'
#' @description
#' This package has three functions that characterize the multivariate responses and covariates under linear or nonlinear structures. All functions are valid to handle measurement error, detection of network in responses, and selection of informative covariates.
#'
#' @details
#' Given the responses and covariates that can be error-prone or precisely measured, NP_graph is a function used to detect the network structure of responses and select important covariates under nonlinear structures. Under linear structures, this package provides two different estimation methods: the function Joint_Gaussian implements the error-corrected  Gaussian maximum likelihood method, and the function Cond_Gaussian extends the neighborhood selection strategy to construct the corrected conditional likelihood function. All functions are able to estimate the network strucure of responses and perform variable selection to identify important covariates with measurement error taken into account.
#' @return
#' PAGE_package
#' @name PAGE_package
#'
NULL
