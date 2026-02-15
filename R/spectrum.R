# Spectral tools: eigenpairs, GFT/IGFT, Chebyshev filtering

#' Compute (and cache) Laplacian eigenpairs
#' @param g gsp_graph
#' @param k number of eigenpairs; NULL for full
#' @param which passed to RSpectra::eigs_sym when k < n (default "SM")
#' @param normalized use normalized Laplacian
#' @export
graph_eigenpairs <- function(g, k = NULL, which = "SM", normalized = g$normalized) {
  L <- graph_laplacian(g, normalized = normalized)
  n <- nrow(L)
  key <- paste0("eig_", normalized, "_", if (is.null(k)) "full" else k, "_", which)
  if (!is.null(g$cache[[key]])) return(g$cache[[key]])

  if (is.null(k) || k >= n) {
    eig <- eigen(as.matrix(L), symmetric = TRUE)
    res <- list(values = eig$values, vectors = eig$vectors)
  } else {
    eig <- RSpectra::eigs_sym(L, k = k, which = which)
    res <- list(values = eig$values, vectors = eig$vectors)
  }
  g$cache[[key]] <- res
  res
}

#' Graph Fourier Transform
#' @param g graph
#' @param signal numeric vector (length n) or matrix (n x m)
#' @param k eigenpairs to use; NULL for full basis
#' @export
gft <- function(g, signal, k = NULL) {
  eig <- graph_eigenpairs(g, k = k)
  U <- eig$vectors
  if (is.null(dim(signal))) {
    signal <- matrix(signal, ncol = 1)
  }
  if (nrow(signal) != nrow(U)) stop("signal length must match number of nodes")
  t(U) %*% signal
}

#' Inverse Graph Fourier Transform
#' @param g graph
#' @param coeffs spectral coefficients (k x m)
#' @param k eigenpairs to use; NULL for full basis
#' @export
igft <- function(g, coeffs, k = NULL) {
  eig <- graph_eigenpairs(g, k = k)
  U <- eig$vectors
  if (is.null(dim(coeffs))) {
    coeffs <- matrix(coeffs, ncol = 1)
  }
  U %*% coeffs
}

#' Estimate largest eigenvalue (spectral radius)
#' @param g a gsp_graph object
#' @param normalized logical; use normalized Laplacian
#' @return the largest eigenvalue
#' @export
lambda_max <- function(g, normalized = g$normalized) {
  key <- paste0("lmax_", normalized)
  if (!is.null(g$cache[[key]])) return(g$cache[[key]])
  L <- graph_laplacian(g, normalized = normalized)
  val <- RSpectra::eigs_sym(L, k = 1, which = "LA")$values[1]
  g$cache[[key]] <- val
  val
}

#' Chebyshev polynomial filtering
#' @param g graph
#' @param signal numeric vector/matrix (n x m)
#' @param coeffs Chebyshev coefficients
#' @param lambda_max_opt optional precomputed lambda_max
#' @export
chebyshev_filter <- function(g, signal, coeffs, lambda_max_opt = NULL) {
  if (is.null(dim(signal))) {
    signal <- matrix(signal, ncol = 1)
  }
  L <- graph_laplacian(g, normalized = g$normalized)
  lmax <- lambda_max_opt %||% lambda_max(g, normalized = g$normalized)
  out <- chebyshev_filter_cpp(L, signal, coeffs, lmax)
  if (ncol(out) == 1) drop(out) else out
}
