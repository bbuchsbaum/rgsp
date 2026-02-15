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

# Tree multiresolution -------------------------------------------------------

.tree_depths_parents <- function(adjacency, root = 1L) {
  if (!inherits(adjacency, "dgCMatrix")) {
    adjacency <- methods::as(adjacency, "dgCMatrix")
  }
  n <- nrow(adjacency)
  if (ncol(adjacency) != n) stop("adjacency must be square")
  root <- as.integer(root)[1]
  if (root < 1L || root > n) stop("root must be in [1, n]")

  p <- adjacency@p
  i <- adjacency@i
  x <- adjacency@x

  depth <- rep.int(-1L, n)
  parent <- integer(n)
  w_parent <- numeric(n)

  queue <- integer(n)
  head <- 1L
  tail <- 1L
  queue[tail] <- root
  depth[root] <- 0L
  parent[root] <- 0L
  w_parent[root] <- NA_real_

  while (head <= tail) {
    v <- queue[head]
    head <- head + 1L

    if (p[v] < p[v + 1L]) {
      idx <- (p[v] + 1L):p[v + 1L]
      nbr <- i[idx] + 1L
      w <- x[idx]

      for (k in seq_along(nbr)) {
        u <- nbr[k]
        if (depth[u] < 0L) {
          depth[u] <- depth[v] + 1L
          parent[u] <- v
          w_parent[u] <- w[k]
          tail <- tail + 1L
          queue[tail] <- u
        }
      }
    }
  }

  if (any(depth < 0L)) stop("tree_multiresolution(): graph is not connected")
  list(depth = depth, parent = parent, w_parent = w_parent)
}

#' Tree multiresolution by even-depth downsampling
#'
#' Builds a multiresolution sequence from a tree by keeping even-depth vertices
#' (with respect to a chosen root) and coarsening by connecting each kept node
#' to its grandparent in the original tree. Edge weights are combined using a
#' simple rule.
#'
#' @param g A `gsp_graph` representing an undirected, connected tree.
#' @param levels Number of coarsening levels.
#' @param root Root vertex index for defining depths (1-based).
#' @param reduction_method How to combine weights along length-2 paths:
#'   `"resistance_distance"` uses series resistance (`1/(1/w1 + 1/w2)`),
#'   `"sum"` uses `w1 + w2`, and `"unweighted"` sets all weights to 1.
#'
#' @return A list with:
#'   \describe{
#'     \item{graphs}{list of graphs at each level (level 0..levels)}
#'     \item{keeps}{list of kept indices at each level (in previous graph)}
#'     \item{roots}{root index at each level (in that level's graph)}
#'   }
#' @export
tree_multiresolution <- function(g,
                                 levels,
                                 root = 1L,
                                 reduction_method = c("resistance_distance", "sum", "unweighted")) {
  if (!inherits(g, "gsp_graph")) stop("'g' must be a gsp_graph")
  if (isTRUE(g$directed)) stop("tree_multiresolution() requires an undirected graph")

  reduction_method <- match.arg(reduction_method)
  levels <- as.integer(levels)[1]
  if (levels < 0) stop("levels must be >= 0")

  A0 <- Matrix::drop0(g$adjacency)
  edges0 <- .edge_list_undirected(A0, min_weight = 0)
  if (length(edges0$w) != (g$n - 1L)) {
    stop("tree_multiresolution(): graph must be a tree (|E| = n-1); got |E|=", length(edges0$w), ", n=", g$n)
  }
  if (!.graph_is_connected(A0)) stop("tree_multiresolution(): graph is not connected")

  graphs <- list(g)
  keeps <- list()
  roots <- integer(levels + 1L)
  roots[1] <- as.integer(root)[1]

  gcur <- g
  root_cur <- roots[1]

  for (lev in seq_len(levels)) {
    A <- Matrix::drop0(gcur$adjacency)
    tp <- .tree_depths_parents(A, root = root_cur)

    keep <- which(tp$depth %% 2L == 0L)
    keep <- sort(keep)
    if (!(root_cur %in% keep)) stop("tree_multiresolution(): root not kept (unexpected)")
    keeps[[lev]] <- keep

    pos <- integer(gcur$n)
    pos[keep] <- seq_along(keep)
    new_root <- pos[root_cur]
    roots[lev + 1L] <- new_root

    # Build edges: connect each kept non-root node to its grandparent
    non_root_keep <- setdiff(keep, root_cur)
    if (length(non_root_keep) == 0L) {
      # Degenerate: only root remains
      adj_new <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(1, 1))
      coords_new <- if (!is.null(gcur$coords)) gcur$coords[root_cur, , drop = FALSE] else NULL
      graphs[[lev + 1L]] <- new_graph(adj_new, coords = coords_new, normalized = gcur$normalized, directed = FALSE)
      gcur <- graphs[[lev + 1L]]
      root_cur <- 1L
      next
    }

    parent1 <- tp$parent[non_root_keep]
    gp <- tp$parent[parent1]
    if (any(parent1 == 0L) || any(gp == 0L)) stop("tree_multiresolution(): encountered missing parent/grandparent")

    w1 <- tp$w_parent[non_root_keep]
    w2 <- tp$w_parent[parent1]

    new_w <- switch(
      reduction_method,
      unweighted = rep(1, length(non_root_keep)),
      sum = w1 + w2,
      resistance_distance = {
        denom <- (1 / w1) + (1 / w2)
        out <- 1 / denom
        out[!is.finite(out)] <- 0
        out
      }
    )

    iu <- pos[non_root_keep]
    iv <- pos[gp]

    adj_new <- Matrix::sparseMatrix(
      i = c(iu, iv),
      j = c(iv, iu),
      x = c(new_w, new_w),
      dims = c(length(keep), length(keep))
    )
    Matrix::diag(adj_new) <- 0
    adj_new <- Matrix::drop0(adj_new)
    adj_new <- Matrix::forceSymmetric(adj_new, uplo = "U")

    coords_new <- if (!is.null(gcur$coords)) gcur$coords[keep, , drop = FALSE] else NULL
    graphs[[lev + 1L]] <- new_graph(adj_new, coords = coords_new, normalized = gcur$normalized, directed = FALSE)

    gcur <- graphs[[lev + 1L]]
    root_cur <- new_root
    if (gcur$n <= 1L) break
  }

  list(graphs = graphs, keeps = keeps, roots = roots)
}
