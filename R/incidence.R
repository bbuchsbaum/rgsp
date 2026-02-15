#' Incidence matrix and edge list
#' @param g gsp_graph
#' @param oriented logical; if TRUE, use i<j orientation with +w at i, -w at j
#' @return list with sparse incidence matrix `B` (m x n) and vectors `rows`, `cols`, `w`
#' @export
graph_incidence <- function(g, oriented = TRUE) {
  adj <- g$adjacency
  rows <- adj@i
  cols <- rep(seq_along(adj@p[-length(adj@p)]) - 1L, diff(adj@p))
  w <- adj@x
  if (oriented) {
    keep <- rows < cols
    rows <- rows[keep]; cols <- cols[keep]; w <- w[keep]
  }
  m <- length(rows)
  n <- g$n
  i <- c(seq_len(m), seq_len(m))
  j <- c(rows + 1L, cols + 1L)
  x <- c(w, -w)
  B <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(m, n))
  list(B = B, rows = rows, cols = cols, w = w)
}

#' Graph gradient (edge differences)
#'
#' Computes weighted edge differences using the oriented incidence matrix.
#'
#' @param g gsp_graph
#' @param signal numeric vector/matrix of node values (n or n x m)
#' @param oriented logical; pass to graph_incidence (TRUE uses i<j orientation)
#' @return matrix of size m x m_signals (or vector length m) with edge gradients
#' @export
grad <- function(g, signal, oriented = TRUE) {
  if (is.null(dim(signal))) signal <- matrix(signal, ncol = 1)
  inc <- graph_incidence(g, oriented = oriented)
  out <- inc$B %*% signal
  if (ncol(out) == 1) drop(out) else out
}

#' Graph divergence (negative incidence transpose)
#'
#' Computes divergence of an edge signal back to nodes: -B^T * edge_signal.
#'
#' @param g gsp_graph
#' @param edge_signal numeric vector/matrix of edge values (m or m x k)
#' @param oriented logical; pass to graph_incidence (must match grad orientation)
#' @return matrix of size n x k (or vector length n)
#' @export
div <- function(g, edge_signal, oriented = TRUE) {
  if (is.null(dim(edge_signal))) edge_signal <- matrix(edge_signal, ncol = 1)
  inc <- graph_incidence(g, oriented = oriented)
  out <- -Matrix::t(inc$B) %*% edge_signal
  if (ncol(out) == 1) drop(out) else out
}

#' Precompute differential operators
#'
#' Returns a list containing the incidence matrix and closures for grad/div.
#'
#' @param g gsp_graph
#' @param oriented logical; orientation used for incidence
#' @return list with components \code{B} (incidence), \code{grad}, \code{div}
#' @export
compute_differential_operator <- function(g, oriented = TRUE) {
  inc <- graph_incidence(g, oriented = oriented)
  list(
    B = inc$B,
    grad = function(x) {
      if (is.null(dim(x))) x <- matrix(x, ncol = 1)
      out <- inc$B %*% x
      if (ncol(out) == 1) drop(out) else out
    },
    div = function(edge_x) {
      if (is.null(dim(edge_x))) edge_x <- matrix(edge_x, ncol = 1)
      out <- -Matrix::t(inc$B) %*% edge_x
      if (ncol(out) == 1) drop(out) else out
    }
  )
}
