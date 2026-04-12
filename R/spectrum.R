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
  key_t <- paste0("eig_", g$normalized, "_", if (is.null(k)) "full" else k, "_", "SM", "_t")
  Ut <- g$cache[[key_t]]
  if (is.null(Ut)) {
    Ut <- t(U)
    g$cache[[key_t]] <- Ut
  }
  if (is.null(dim(signal))) {
    signal <- matrix(signal, ncol = 1)
  }
  if (nrow(signal) != nrow(U)) stop("signal length must match number of nodes")
  Ut_ptr <- g$cache[[paste0(key_t, "_xptr")]]
  if (is.null(Ut_ptr)) {
    Ut_ptr <- dense_mat_xptr_cpp(Ut)
    g$cache[[paste0(key_t, "_xptr")]] <- Ut_ptr
  }
  dense_matmul_xptr_cpp(Ut_ptr, signal)
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
  key <- paste0("eig_", g$normalized, "_", if (is.null(k)) "full" else k, "_", "SM")
  U_ptr <- g$cache[[paste0(key, "_xptr")]]
  if (is.null(U_ptr)) {
    U_ptr <- dense_mat_xptr_cpp(U)
    g$cache[[paste0(key, "_xptr")]] <- U_ptr
  }
  dense_matmul_xptr_cpp(U_ptr, coeffs)
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
  val <- tryCatch({
    out <- suppressWarnings(RSpectra::eigs_sym(L, k = 1, which = "LA"))
    if (length(out$values) < 1) {
      NA_real_
    } else {
      as.numeric(out$values[1])
    }
  }, error = function(e) NA_real_)

  if (!is.finite(val) || is.na(val)) {
    # Robust upper bounds (Gershgorin) as fallback.
    if (isTRUE(normalized)) {
      val <- 2
    } else {
      deg <- graph_degree(g)
      val <- 2 * max(as.numeric(deg))
    }
  }
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
