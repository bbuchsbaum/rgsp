# Graph-based randomized low-pass bases and compression

#' Default probe block size for Chebyshev filtering
#'
#' Chooses a column block size so that filtering a probe block uses roughly
#' `target_mb` of RAM (in doubles), capped by the number of probe vectors.
#'
#' @param n_nodes number of graph nodes (rows)
#' @param m_probes number of probe vectors (columns)
#' @param target_mb target memory in megabytes (default 300)
#' @return integer block size (columns)
#' @export
default_probe_block_cols <- function(n_nodes, m_probes, target_mb = 300) {
  stopifnot(n_nodes >= 1, m_probes >= 1, target_mb > 0)
  target_bytes <- target_mb * 1e6
  block <- max(1L, min(m_probes, ceiling(target_bytes / (8 * n_nodes))))
  as.integer(block)
}


#' Randomized low-pass basis via Chebyshev filtering (no eigenpairs)
#'
#' Builds a k-dimensional, graph-smooth orthonormal basis by filtering random
#' probes with a low-pass kernel approximated by Chebyshev polynomials. No
#' eigen-decomposition of the Laplacian is required.
#'
#' @param g gsp_graph object
#' @param k target basis size
#' @param kernel low-pass kernel; one of "heat", "step", or a function
#'   accepting lambda. For "heat", provide `tau`; for "step", provide `cutoff`.
#' @param tau diffusion time for heat kernel (ignored unless kernel = "heat").
#'   If NULL, chosen so that lambda_c = 0.05 * lmax (i.e., tau = 1 / lambda_c).
#' @param cutoff step cutoff (keeps lambda <= cutoff); if NULL and kernel =
#'   "step", defaults to 0.05 * lmax.
#' @param K Chebyshev order. If missing, defaults to 30 (or 50 when
#'   `jackson = FALSE`).
#' @param oversample oversampling dimension; if missing, defaults to
#'   `max(16, ceiling(0.1 * k))`.
#' @param jackson logical; apply Jackson damping to Chebyshev coefficients
#' @param lmax optional spectral radius (else computed)
#' @param block optional number of probe columns to filter at once (for memory).
#'   If NULL, chosen to target ~300 MB per filtered probe block (keeps memory in
#'   the 200-500 MB range), capped at the total number of probes.
#' @param seed optional RNG seed for reproducibility
#' @param verbose logical; print progress messages (default FALSE)
#' @return list with elements:
#'   \describe{
#'     \item{loadings}{n x k orthonormal basis (columns)}
#'     \item{lambdas}{length-k Ritz eigenvalue estimates}
#'     \item{lmax}{spectral radius used}
#'     \item{kernel}{kernel identifier}
#'     \item{K}{Chebyshev order}
#'     \item{block}{filtered probe block size used}
#'   }
#' @export
random_lowpass_basis <- function(g, k,
                                 kernel = c("heat", "step"),
                                 tau = NULL,
                                 cutoff = NULL,
                                 K = 30,
                                 oversample = 16,
                                 jackson = TRUE,
                                 lmax = NULL,
                                 block = NULL,
                                 seed = NULL,
                                 verbose = FALSE) {
  kernel <- if (is.function(kernel)) kernel else match.arg(kernel)

  if (missing(oversample)) {
    oversample <- max(16L, ceiling(0.1 * k))
  }
  if (missing(K)) {
    K <- if (isFALSE(jackson)) 50L else 30L
  }
  stopifnot(k >= 1, oversample >= 0, K >= 1)

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit({
      if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv)
      else assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
  }

  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)

  kernel_fun <- switch(
    if (is.function(kernel)) "custom" else kernel,
    heat = {
      if (is.null(tau)) {
        lambda_c <- 0.05 * lmax
        tau <- 1 / lambda_c
      }
      kernel_heat(tau)
    },
    step = {
      if (is.null(cutoff)) cutoff <- 0.05 * lmax
      kernel_rectangle(0, cutoff)
    },
    custom = kernel
  )

  m <- k + oversample
  n <- g$n

  # Random probes (Rademacher to reduce variance)
  Z <- matrix(sample(c(-1, 1), size = n * m, replace = TRUE), nrow = n, ncol = m)

  # Chebyshev coefficients for low-pass kernel
  cheb <- cheby_coeffs(kernel_fun, K = K, lmax = lmax, jackson = jackson)

  L <- graph_laplacian(g, g$normalized)

  # Filter in blocks if requested
  if (is.null(block)) {
    block <- default_probe_block_cols(n, m, target_mb = 300)
  }
  if (isTRUE(verbose)) {
    est_mb <- round(8 * n * block / 1e6, 1)
    cli::cli_inform("random_lowpass_basis: using block={block} (~{est_mb} MB) with m={m}, n={n}")
  }
  Y_parts <- vector("list", ceiling(m / block))
  start <- 1
  idx <- 1
  while (start <= m) {
    end <- min(start + block - 1, m)
    Y_parts[[idx]] <- chebyshev_filter_cpp(L, Z[, start:end, drop = FALSE], cheb, lmax)
    start <- end + 1
    idx <- idx + 1
  }
  Y <- do.call(cbind, Y_parts)

  # Orthonormalize
  Q <- qr.Q(qr(Y))
  if (ncol(Q) < k) {
    stop("Obtained basis has fewer columns than requested k; increase oversample or K.")
  }

  # Rayleigh-Ritz refinement
  B <- Matrix::crossprod(Q, L %*% Q)
  eig_small <- eigen(B, symmetric = TRUE)
  order_idx <- order(eig_small$values)
  Vk <- eig_small$vectors[, order_idx[seq_len(k)], drop = FALSE]
  lambdas <- eig_small$values[order_idx[seq_len(k)]]

  U <- Q %*% Vk

  list(
    loadings = U,
    lambdas = lambdas,
    lmax = lmax,
    kernel = if (is.function(kernel)) "custom" else kernel,
    K = K,
    block = block
  )
}


#' Compress signals into a graph-smooth basis (no eigendecomposition)
#'
#' Convenience wrapper that builds a randomized low-pass basis and projects
#' node signals into that basis.
#'
#' @param g gsp_graph object
#' @param X matrix with columns = nodes (ncols must equal g$n); rows are samples
#'   (e.g., time points). A numeric vector is treated as a single sample.
#' @param k target basis size
#' @param block_rows optional row-block size when projecting X (to limit memory)
#' @param ... passed to [random_lowpass_basis()]
#' @return list with elements:
#'   \describe{
#'     \item{basis}{projection of X into k-dim space (n_samples x k)}
#'     \item{loadings}{n_nodes x k basis matrix}
#'     \item{lambdas}{Ritz eigenvalue estimates}
#'   }
#' @export
graph_compress <- function(g, X, k, block_rows = NULL, ...) {
  if (is.null(dim(X))) {
    X <- matrix(X, nrow = 1)
  }
  if (ncol(X) != g$n) {
    stop("ncol(X) must equal number of graph nodes (", g$n, ")")
  }

  basis_res <- random_lowpass_basis(g, k = k, ...)
  U <- basis_res$loadings

  # Project X into the basis, optionally in row blocks
  n_rows <- nrow(X)
  block_rows <- block_rows %||% n_rows
  B <- matrix(0, nrow = n_rows, ncol = k)

  start <- 1
  while (start <= n_rows) {
    end <- min(start + block_rows - 1, n_rows)
    B[start:end, ] <- X[start:end, , drop = FALSE] %*% U
    start <- end + 1
  }

  basis_res$basis <- B
  basis_res
}
