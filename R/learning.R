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

#' Graph Tikhonov classification with simplex constraints
#'
#' Semi-supervised classification with missing labels (mask). Solves
#'
#' \deqn{\arg\min_X \|M X - Y\|_F^2 + \tau \, \mathrm{tr}(X^\top L X) \quad
#' \text{s.t. each row of } X \text{ lies on the simplex}}
#'
#' where `Y` is a one-hot encoding of the observed labels and `M` is a diagonal
#' masking operator (observed rows).
#'
#' This is a port of `pygsp.learning.classification_tikhonov_simplex` using a
#' projected (proximal) gradient method.
#'
#' @param g graph
#' @param y integer/factor labels of length `g$n`; missing labels may be `NA`
#' @param mask logical vector of observed nodes (length `g$n`); if `NULL`,
#'   defaults to `!is.na(y)`
#' @param tau regularization strength (must be > 0)
#' @param tol stopping tolerance on relative iterate change
#' @param maxit maximum iterations
#' @param lmax optional spectral radius for step size (defaults to `lambda_max(g)`)
#' @param step optional step size; if `NULL`, uses `0.5/(1 + tau*lmax)`
#'
#' @return matrix of class logits (n x K) with rows summing to 1 and non-negative
#' @export
classification_tikhonov_simplex <- function(g,
                                            y,
                                            mask = NULL,
                                            tau = 0.1,
                                            tol = 1e-6,
                                            maxit = 200,
                                            lmax = NULL,
                                            step = NULL) {
  if (!inherits(g, "gsp_graph")) stop("'g' must be a gsp_graph", call. = FALSE)
  if (!is.finite(tau) || tau <= 0) stop("tau must be > 0", call. = FALSE)

  n <- g$n
  y0 <- y
  if (is.factor(y0)) y0 <- as.integer(y0)
  y0 <- as.integer(y0)
  if (length(y0) != n) stop("y must have length g$n", call. = FALSE)

  if (is.null(mask)) {
    mask <- !is.na(y0)
  }
  mask <- as.logical(mask)
  if (length(mask) != n) stop("mask must have length g$n", call. = FALSE)

  obs <- which(mask)
  if (length(obs) == 0L) stop("mask selects no observed nodes", call. = FALSE)

  classes <- sort(unique(y0[obs]))
  if (any(is.na(classes))) stop("observed labels must be non-missing", call. = FALSE)
  k <- length(classes)
  if (k < 2L) stop("need at least 2 classes in observed labels", call. = FALSE)

  y_idx <- match(y0, classes)

  Y <- matrix(0, nrow = n, ncol = k)
  for (ii in obs) {
    Y[ii, y_idx[ii]] <- 1
  }

  L <- graph_laplacian(g, normalized = g$normalized)

  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)
  if (!is.finite(lmax) || lmax <= 0) stop("failed to estimate lmax", call. = FALSE)

  step <- step %||% (0.5 / (1 + tau * lmax))
  step <- as.numeric(step)[1]
  if (!is.finite(step) || step <= 0) stop("step must be > 0", call. = FALSE)

  mask_num <- as.numeric(mask)

  X <- Y
  # Ensure all rows start on simplex (unobserved rows in particular).
  X <- proj_simplex_rows_cpp(X)

  for (iter in seq_len(as.integer(maxit))) {
    X_old <- X

    LX <- as.matrix(L %*% X)
    grad <- (X - Y) * mask_num + tau * LX
    X <- proj_simplex_rows_cpp(X - (2 * step) * grad)

    diff <- sqrt(sum((X - X_old)^2))
    base <- sqrt(sum(X_old^2))
    if (diff < tol * max(1, base)) break
  }

  colnames(X) <- paste0("class_", classes)
  X
}
