# Graph coarsening and pyramid transforms

#' Kron reduction (Schur complement) of a graph Laplacian
#'
#' Reduces a graph by eliminating a subset of nodes via Schur complement.
#'
#' @param g gsp_graph
#' @param keep integer vector of node indices to keep (1-based)
#' @return reduced gsp_graph on |keep| nodes
#' @export
kron_reduction <- function(g, keep) {
  keep <- as.integer(keep)
  n <- g$n
  if (any(keep < 1 | keep > n)) stop("keep indices must be within [1, n]")
  keep <- sort(unique(keep))
  elim <- setdiff(seq_len(n), keep)
  if (length(elim) == 0) return(g)

  L <- graph_laplacian(g, normalized = FALSE)
  Lkk <- L[keep, keep, drop = FALSE]
  Lkc <- L[keep, elim, drop = FALSE]
  Lcc <- L[elim, elim, drop = FALSE]

  # Solve Schur complement: Lred = Lkk - Lkc %*% solve(Lcc, Lck)
  Lcc_inv_Lck <- Matrix::solve(Lcc, Matrix::t(Lkc))
  Lred <- Lkk - Lkc %*% Lcc_inv_Lck

  # Convert back to adjacency: A = D - L (combinatorial)
  deg <- Matrix::diag(Lred)
  Ared <- Matrix::Diagonal(x = deg) - Lred
  Matrix::diag(Ared) <- 0
  Ared <- Matrix::forceSymmetric(Ared, uplo = "U")
  Ared <- Matrix::drop0(Ared)

  coords <- if (!is.null(g$coords)) g$coords[keep, , drop = FALSE] else NULL
  new_graph(Ared, coords = coords, normalized = g$normalized, directed = FALSE)
}

#' Multi-level graph coarsening
#'
#' Builds a hierarchy by repeated random-node Kron reduction.
#'
#' @param g gsp_graph
#' @param levels number of coarsening levels
#' @param keep_fraction fraction of nodes to keep each level (0,1)
#' @param seed optional seed for reproducibility
#' @return list with elements `graphs` (level 0..L) and `keeps` (indices per level)
#' @export
graph_multiresolution <- function(g, levels = 1, keep_fraction = 0.5, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  graphs <- list(g)
  keeps <- list(seq_len(g$n))
  gcur <- g
  for (lvl in seq_len(levels)) {
    n <- gcur$n
    k <- max(1, ceiling(n * keep_fraction))
    keep <- sort(sample.int(n, k, replace = FALSE))
    gcur <- kron_reduction(gcur, keep)
    graphs[[lvl + 1]] <- gcur
    keeps[[lvl + 1]] <- keep
    if (gcur$n <= 2) break
  }
  list(graphs = graphs, keeps = keeps)
}

#' Pyramid analysis (coarse + detail) via simple decimation
#'
#' @param g gsp_graph
#' @param signal numeric vector length n
#' @param keep indices kept (if NULL, keep every other node)
#' @return list with coarse, detail, keep, elim
#' @export
pyramid_analysis <- function(g, signal, keep = NULL) {
  n <- g$n
  if (is.null(keep)) {
    keep <- seq(1, n, by = 2)
  }
  keep <- sort(unique(as.integer(keep)))
  elim <- setdiff(seq_len(n), keep)
  coarse <- signal[keep]
  detail <- signal[elim]
  list(coarse = coarse, detail = detail, keep = keep, elim = elim)
}

#' Pyramid synthesis (reconstruction) from coarse/detail
#'
#' @param p list from pyramid_analysis
#' @return reconstructed signal
#' @export
pyramid_synthesis <- function(p) {
  n <- length(p$coarse) + length(p$detail)
  x <- numeric(n)
  x[p$keep] <- p$coarse
  x[p$elim] <- p$detail
  x
}
