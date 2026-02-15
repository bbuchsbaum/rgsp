# Learning / optimization primitives

#' L1 soft-threshold (proximal)
#' @param x vector or matrix
#' @param tau threshold parameter
#' @export
prox_l1 <- function(x, tau) {
  prox_l1_cpp(as.matrix(x), tau)
}

#' 1D total variation proximal
#' @param y numeric vector
#' @param lambda weight
#' @export
prox_tv1d <- function(y, lambda) {
  prox_tv1d_cpp(as.numeric(y), lambda)
}

#' Box projection proximal (component-wise clamp)
#' @param y numeric vector
#' @param lower lower bound
#' @param upper upper bound
#' @return numeric vector with each element clamped to the interval
#' @export
prox_box <- function(y, lower, upper) {
  out <- prox_box_cpp(as.numeric(y), lower, upper)
  drop(out)
}

#' Simplex projection
#' @param y numeric vector
#' @export
proj_simplex <- function(y) {
  drop(proj_simplex_cpp(as.numeric(y)))
}

#' l2-ball projection
#' @param y numeric vector
#' @param r radius
#' @export
proj_l2_ball <- function(y, r) {
  drop(proj_l2_ball_cpp(as.numeric(y), r))
}

#' Tikhonov smoothing: solve (I + tau L) x = b
#' @param g graph
#' @param b vector or matrix of observations
#' @param tau smoothing strength
#' @param tol CG tolerance
#' @param maxit CG max iterations
#' @param threads optional number of threads (reserved for future parallel solver)
#' @export
tikhonov_smooth <- function(g, b, tau, tol = 1e-6, maxit = 1000, threads = NULL) {
  if (is.null(dim(b))) b <- matrix(b, ncol = 1)
  L <- graph_laplacian(g, normalized = g$normalized)
  x <- tikhonov_solve_cpp(L, b, tau = tau, tol = tol, maxit = maxit)
  if (ncol(x) == 1) drop(x) else x
}

#' Interpolation / inpainting of missing nodes via Laplacian smoothing
#' @param g graph
#' @param y vector of length n with NA for missing entries
#' @param alpha smoothing parameter
#' @param tol CG tolerance
#' @param maxit CG max iterations
#' @param threads optional number of threads (reserved)
#' @export
interpolate_laplacian <- function(g, y, alpha = 1e-2, tol = 1e-6, maxit = 1000, threads = NULL) {
  y_vec <- as.numeric(y)
  mask <- !is.na(y_vec)
  b <- y_vec
  b[!mask] <- 0
  L <- graph_laplacian(g, normalized = g$normalized)
  # M = diag(mask)
  Mdiag <- Matrix::Diagonal(x = as.numeric(mask))
  A <- Mdiag + alpha * L
  rhs <- b
  # convert to dense matrix for cg helper
  x <- tikhonov_solve_cpp(A, matrix(rhs, ncol = 1), tau = 0, tol = tol, maxit = maxit)
  drop(x)
}

#' TV-based inpainting (missing entries) using primal-dual prox
#' @param g graph
#' @param y length-n vector with NA for missing entries
#' @param lambda TV weight
#' @param maxit iterations
#' @param tau primal step
#' @param sigma dual step
#' @param tol stopping tolerance
#' @export
tv_inpaint <- function(g, y, lambda, maxit = 300, tau = 0.5, sigma = 0.5, tol = 1e-6) {
  y_vec <- as.numeric(y)
  mask_obs <- !is.na(y_vec)
  y_filled <- y_vec
  y_filled[!mask_obs] <- 0

  inc <- graph_incidence(g)
  rows <- inc$rows; cols <- inc$cols; w <- inc$w

  x <- y_filled
  # projected gradient: data term only on observed entries via mask
  for (iter in seq_len(maxit)) {
    x_old <- x
    x <- prox_graph_tv_cpp(x, rows, cols, w, lambda, maxit = 50, tau = tau, sigma = sigma, tol = tol)
    # enforce data fidelity on observed
    x[mask_obs] <- y_vec[mask_obs]
    if (sqrt(sum((x - x_old)^2)) < tol * max(1, sqrt(sum(x_old^2)))) break
  }
  x
}
#' Graph TV denoising via primal-dual (Chambolle-Pock)
#' @param g graph
#' @param y numeric vector length n
#' @param lambda TV weight
#' @param maxit iterations
#' @param tau primal step
#' @param sigma dual step
#' @param tol stopping tolerance
#' @export
tv_denoise <- function(g, y, lambda, maxit = 300, tau = 0.5, sigma = 0.5, tol = 1e-6) {
  y <- as.numeric(y)
  # build edge list from adjacency
  adj <- g$adjacency
  rows <- adj@i
  cols <- rep(seq_along(adj@p[-length(adj@p)]) - 1L, diff(adj@p))
  w <- adj@x
  # use only upper-tri to avoid double edges
  keep <- rows < cols
  rows <- rows[keep]; cols <- cols[keep]; w <- w[keep]
  prox_graph_tv_cpp(y, rows, cols, w, lambda, maxit = maxit, tau = tau, sigma = sigma, tol = tol)
}
