#' Sensitivity analysis on binary regression coefficients
#'
#' Sensitivity analysis on regression coefficients to assess how measurement error
#' (misclassification) in a binary variable affects the estimated treatment effect.
#' It calculates extreme bounds on the coefficient by flipping binary classifications
#' in the data.
#'
#' @docType package
#' @author Leonardo Dantas \email{leodantas@princeton.edu}
#' @name futureproofR
#'
#' @keywords package
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom tidyr pivot_longer spread drop_na
#' @importFrom purrr map
#' @importFrom stringr str_extract
#' @importFrom stats coef lm update formula quantile rbinom rnorm cor complete.cases variable.names resid model.frame na.omit confint
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom mvtnorm rmvnorm
#' @importFrom condMVNorm rcmvnorm
NULL
