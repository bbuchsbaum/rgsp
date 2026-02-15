# High-level learning wrappers

#' Graph Tikhonov regression
#'
#' Solves min_x 0.5||y - x||^2 + tau x^T L x
#'
#' @param g graph
#' @param y observed labels/signals (n or n x m)
#' @param tau regularization strength
#' @param tol solver tolerance
#' @param maxit max iterations
#' @param threads optional threads
#' @export
regression_tikhonov <- function(g, y, tau, tol = 1e-6, maxit = 1000, threads = NULL) {
  tikhonov_smooth(g, y, tau = tau, tol = tol, maxit = maxit, threads = threads)
}

#' Graph Tikhonov classification (one-vs-rest smoothing)
#'
#' Smooths one-hot label columns over the graph Laplacian.
#'
#' @param g graph
#' @param labels factor or integer vector length n
#' @param tau regularization strength
#' @param tol solver tolerance
#' @param maxit max iterations
#' @param threads optional threads
#' @return matrix of smoothed class scores (n x K)
#' @export
classification_tikhonov <- function(g, labels, tau, tol = 1e-6, maxit = 1000, threads = NULL) {
  lab <- labels
  if (is.factor(lab)) lab <- as.integer(lab)
  classes <- sort(unique(lab))
  n <- g$n
  Y <- matrix(0, nrow = n, ncol = length(classes))
  for (i in seq_along(classes)) {
    Y[, i] <- as.numeric(lab == classes[i])
  }
  scores <- tikhonov_smooth(g, Y, tau = tau, tol = tol, maxit = maxit, threads = threads)
  colnames(scores) <- paste0("class_", classes)
  scores
}
