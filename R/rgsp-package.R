#' rgsp: Graph Signal Processing in R (Port of PyGSP)
#'
#' Core graph signal processing algorithms with C++ acceleration via RcppArmadillo.
#'
#' @keywords internal
#' @useDynLib rgsp, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import Matrix
#' @import methods
#' @importFrom stats runif rnorm rbinom dist setNames
#' @importFrom utils combn
#' @importFrom grDevices hcl.colors
#' @importFrom graphics axis barplot hist par
#' @importFrom rlang .data
"_PACKAGE"

# Avoid R CMD check note for .data pronoun used in ggplot2 aes()
utils::globalVariables(".data")
