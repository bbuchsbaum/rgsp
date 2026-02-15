# Spectral sparsification utilities (Spielman-Srivastava style)

.graph_is_connected <- function(adjacency) {
  if (!inherits(adjacency, "dgCMatrix")) {
    adjacency <- methods::as(adjacency, "dgCMatrix")
  }
  n <- nrow(adjacency)
  if (n <= 1) return(TRUE)

  p <- adjacency@p
  i <- adjacency@i

  seen <- rep(FALSE, n)
  stack <- integer(n)

  # Start from first node with any incident edge.
  start <- which(diff(p) > 0L)[1]
  if (is.na(start)) return(FALSE)

  top <- 1L
  stack[top] <- start
  seen[start] <- TRUE

  while (top > 0L) {
    v <- stack[top]
    top <- top - 1L

    if (p[v] < p[v + 1L]) {
      nbr <- i[(p[v] + 1L):p[v + 1L]] + 1L
      for (u in nbr) {
        if (!seen[u]) {
          seen[u] <- TRUE
          top <- top + 1L
          stack[top] <- u
        }
      }
    }
  }

  all(seen)
}

.edge_list_undirected <- function(adjacency, min_weight = 0) {
  adjacency <- Matrix::drop0(adjacency)
  if (!inherits(adjacency, "dgCMatrix")) {
    adjacency <- methods::as(adjacency, "dgCMatrix")
  }
  s <- Matrix::summary(Matrix::triu(adjacency, 1))
  i <- as.integer(s$i)
  j <- as.integer(s$j)
  w <- as.numeric(s$x)
  keep <- w > min_weight
  list(i = i[keep], j = j[keep], w = w[keep])
}

.approx_effective_resistance_edges <- function(L, edges, n_probes, reg) {
  stopifnot(inherits(L, "dgCMatrix"))
  n <- nrow(L)
  stopifnot(ncol(L) == n)

  if (n_probes < 1) stop("n_probes must be >= 1")
  if (!is.finite(reg) || reg <= 0) stop("reg must be > 0")

  L_reg <- L + Matrix::Diagonal(n, x = rep(reg, n))
  chol <- Matrix::Cholesky(L_reg, LDL = TRUE, perm = TRUE)

  G <- matrix(sample(c(-1, 1), size = n * n_probes, replace = TRUE), nrow = n, ncol = n_probes)
  X <- Matrix::solve(chol, G)

  i <- edges$i
  j <- edges$j
  block <- 10000L
  R <- numeric(length(i))

  start <- 1L
  while (start <= length(i)) {
    end <- min(start + block - 1L, length(i))
    Xi <- X[i[start:end], , drop = FALSE]
    Xj <- X[j[start:end], , drop = FALSE]
    R[start:end] <- rowMeans((Xi - Xj)^2)
    start <- end + 1L
  }

  pmax(R, 0)
}

#' Spectral graph sparsification (Spielman–Srivastava edge sampling)
#'
#' Sparsifies an undirected graph by sampling edges proportionally to
#' approximate effective resistance, yielding a reweighted sparse graph whose
#' Laplacian approximately preserves the quadratic form.
#'
#' This is a practical, R-native approximation intended to reduce `nnz` and
#' speed up downstream routines (Chebyshev filtering, TV, diffusion, etc.).
#'
#' @param g A `gsp_graph` (undirected).
#' @param epsilon Sparsification parameter in `[1/sqrt(n), 1)`. Smaller values
#'   sample more edges (closer to original).
#' @param maxiter Maximum number of attempts. If a sparsified graph is
#'   disconnected, `epsilon` is reduced toward `1/sqrt(n)` and retried.
#' @param seed Optional RNG seed for reproducibility.
#' @param n_probes Number of random probe vectors used to approximate effective
#'   resistances (default scales with `log(n)` and is capped).
#' @param reg Diagonal regularization added to the combinatorial Laplacian for
#'   linear solves (default scaled to graph degree).
#' @param q Number of edge samples (with replacement). If `NULL`, uses a
#'   Spielman–Srivastava-inspired heuristic `q ~ 0.16 * n * log(n) / epsilon^2`.
#'
#' @return A new `gsp_graph` with a sparser adjacency.
#' @export
graph_sparsify_spectral <- function(g,
                                    epsilon,
                                    maxiter = 10,
                                    seed = NULL,
                                    n_probes = NULL,
                                    reg = NULL,
                                    q = NULL) {
  if (!inherits(g, "gsp_graph")) stop("'g' must be a gsp_graph")
  if (isTRUE(g$directed)) stop("graph_sparsify_spectral() requires an undirected graph")

  n <- g$n
  if (n < 2) return(g)

  if (!is.finite(epsilon) || epsilon <= 0) stop("epsilon must be > 0")
  eps_min <- 1 / sqrt(n)
  if (epsilon < eps_min || epsilon >= 1) {
    stop("epsilon must be in [1/sqrt(n), 1); for n=", n, ", min epsilon=", signif(eps_min, 4))
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }

  adjacency <- Matrix::drop0(g$adjacency)
  edges <- .edge_list_undirected(adjacency, min_weight = 1e-10)
  if (length(edges$w) == 0) stop("graph has no edges")

  L <- graph_laplacian(g, normalized = FALSE)

  if (is.null(n_probes)) {
    n_probes <- max(16L, min(64L, as.integer(ceiling(8 * log(n)))))
  } else {
    n_probes <- as.integer(n_probes)
  }

  if (is.null(reg)) {
    deg_scale <- max(1, stats::median(as.numeric(Matrix::diag(L)), na.rm = TRUE))
    reg <- 1e-6 * deg_scale
  }

  attempt_eps <- epsilon

  for (attempt in seq_len(maxiter)) {
    # approximate effective resistances and sample edges
    Re <- .approx_effective_resistance_edges(L, edges, n_probes = n_probes, reg = reg)
    Pe <- pmax(0, edges$w) * pmax(0, Re)
    sPe <- sum(Pe)
    if (!is.finite(sPe) || sPe <= 0) stop("failed to form sampling distribution (sum=0)")
    Pe <- Pe / sPe

    q_eff <- q
    if (is.null(q_eff)) {
      q_eff <- ceiling(0.16 * n * log(n) / (attempt_eps^2))
    }
    q_eff <- max(as.integer(q_eff), n - 1L)

    draws <- sample.int(length(Pe), size = q_eff, replace = TRUE, prob = Pe)
    counts <- tabulate(draws, nbins = length(Pe))
    keep <- which(counts > 0L)

    new_w <- counts[keep] * edges$w[keep] / (q_eff * Pe[keep])
    new_i <- edges$i[keep]
    new_j <- edges$j[keep]

    W_sparse <- Matrix::sparseMatrix(
      i = c(new_i, new_j),
      j = c(new_j, new_i),
      x = c(new_w, new_w),
      dims = c(n, n)
    )
    Matrix::diag(W_sparse) <- 0
    W_sparse <- Matrix::drop0(W_sparse)
    W_sparse <- Matrix::forceSymmetric(W_sparse, uplo = "U")

    if (.graph_is_connected(W_sparse)) {
      out <- new_graph(
        methods::as(W_sparse, "dgCMatrix"),
        coords = g$coords,
        normalized = g$normalized,
        directed = FALSE
      )
      out$sparsify <- list(
        method = "spectral",
        epsilon = attempt_eps,
        n_probes = n_probes,
        reg = reg,
        q = q_eff,
        edges_in = length(edges$w),
        edges_out = length(.edge_list_undirected(out$adjacency)$w)
      )
      return(out)
    }

    if (attempt == maxiter) {
      warning("graph_sparsify_spectral(): sparsified graph disconnected after maxiter attempts")
      break
    }

    # Move epsilon toward eps_min to sample more edges.
    attempt_eps <- attempt_eps - (attempt_eps - eps_min) / 2
  }

  # Fallback: return original graph unchanged.
  g
}

