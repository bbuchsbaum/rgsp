# Cache invalidation utilities

#' Invalidate cached spectra (eigenpairs, lambda_max, Chebyshev coeffs)
#' @param g a gsp_graph object
#' @return the graph with spectral cache cleared
#' @export
graph_cache_invalidate_spectra <- function(g) {
  if (is.environment(g$cache)) {
    to_drop <- ls(envir = g$cache, all.names = TRUE)
    to_drop <- to_drop[grepl("^(eig_|lmax_|cheby_)", to_drop)]
    if (length(to_drop) > 0) {
      rm(list = to_drop, envir = g$cache)
    }
  } else if (!is.null(g$cache)) {
    g$cache[grepl("^eig_|^lmax_|^cheby_", names(g$cache))] <- NULL
  }

  cheby_cache_clear()
  g
}

#' Invalidate cached Laplacians and degrees
#' @param g a gsp_graph object
#' @return the graph with structural cache cleared
#' @export
graph_cache_invalidate_struct <- function(g) {
  if (is.environment(g$cache)) {
    to_drop <- ls(envir = g$cache, all.names = TRUE)
    to_drop <- to_drop[grepl("^(lap_|degree$|laplacian_sp$)", to_drop)]
    if (length(to_drop) > 0) {
      rm(list = to_drop, envir = g$cache)
    }
  } else if (!is.null(g$cache)) {
    g$cache[grepl("^lap_|^degree", names(g$cache))] <- NULL
  }
  g$degree <- NULL
  g$laplacian <- NULL
  g
}

#' Invalidate all graph caches (structural + spectral)
#'
#' Convenience helper that clears Laplacian/degree, spectral radius, eigenpairs,
#' frame bounds, and any hash used for cache keys. Intended for use after
#' modifying adjacency/weights.
#' @param g a gsp_graph object
#' @return the graph with caches cleared
#' @export
graph_invalidate_cache <- function(g) {
  g <- graph_cache_invalidate_struct(g)
  g <- graph_cache_invalidate_spectra(g)
  g$frame_bounds <- NULL
  g$hash <- NULL
  g
}

#' Clear global Chebyshev coefficient cache
#' @return NULL invisibly
#' @export
cheby_cache_clear <- function() {
  rm(list = ls(envir = .cheby_cache), envir = .cheby_cache)
  invisible(NULL)
}

#' Validate signal dimensions against graph
#' @param g graph
#' @param signal vector or matrix
#' @keywords internal
validate_signal_dims <- function(g, signal) {
  n <- g$n
  if (is.null(dim(signal))) {
    if (length(signal) != n) stop("signal length must match number of nodes")
  } else {
    if (nrow(signal) != n) stop("signal rows must match number of nodes")
  }
}
