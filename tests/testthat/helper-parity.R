# Shared helpers for golden-value parity tests
# Loaded automatically by testthat for all test-parity-golden-*.R files

# ---------------------------------------------------------------------------
# Golden value loading
# ---------------------------------------------------------------------------

#' Load a golden reference JSON file from inst/extdata/parity/
#' @param filename basename of the JSON file (e.g., "golden_graphs.json")
#' @return parsed list
load_golden <- function(filename) {

  # During R CMD check, system.file works; during devtools::test, try inst/ directly

  path <- system.file("extdata", "parity", filename, package = "rgsp")
  if (path == "") {
    # Fallback for devtools::test
    path <- file.path("../../inst/extdata/parity", filename)
  }
  if (!file.exists(path)) {
    return(NULL)
  }
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

#' Skip test if golden file is missing
skip_if_no_golden <- function(filename) {
  ref <- load_golden(filename)
  if (is.null(ref)) {
    testthat::skip(paste0("Golden file not found: ", filename))
  }
  ref
}

# ---------------------------------------------------------------------------
# API mapping functions (PyGSP -> rgsp)
# ---------------------------------------------------------------------------

#' Convert PyGSP heat scale to rgsp t parameter
#' PyGSP: exp(-scale * lambda / lmax)
#' rgsp:  exp(-t * lambda)
#' Therefore: t = scale / lmax
pygsp_heat_to_rgsp <- function(scale, lmax) {
  scale / lmax
}

#' Convert PyGSP MexicanHat scale to rgsp sigma
#' PyGSP bandpass: scale * x * exp(-scale * x)
#' rgsp:           x * exp(-x / sigma)
#' Matching: sigma = 1 / scale
pygsp_mexhat_scale_to_rgsp_sigma <- function(scale) {
  1 / scale
}

# ---------------------------------------------------------------------------
# Tolerance profiles
# ---------------------------------------------------------------------------

tol_eigenvalues       <- 1e-10   # Exact computation, same graph
tol_gft               <- 1e-10   # Eigenvector sign ambiguity resolved via abs sort
tol_filter_exact      <- 1e-8    # Floating point accumulation
tol_filter_chebyshev  <- 1e-3    # Polynomial approximation error
tol_correlation       <- 0.95    # When tolerances differ too much
tol_cheby_coeffs      <- 1e-6    # Different quadrature resolution
tol_differential      <- 1e-10   # Exact sparse matrix
tol_frame_bounds      <- 0.05    # Estimation method differences
tol_tikhonov          <- 1e-3    # Iterative solver convergence
tol_kron              <- 1e-8    # Schur complement, exact
tol_lmax              <- 0.02    # PyGSP inflates by 1.01

# ---------------------------------------------------------------------------
# Custom comparison utilities
# ---------------------------------------------------------------------------

#' Compare two vectors up to element-wise sign flips
#' Useful for eigenvector comparisons where sign is ambiguous
expect_equal_up_to_sign <- function(actual, expected, tolerance = 1e-10) {
  # Check if either sign matches
  diff_pos <- max(abs(actual - expected))
  diff_neg <- max(abs(actual + expected))
  best_diff <- min(diff_pos, diff_neg)
  testthat::expect_true(
    best_diff < tolerance,
    label = paste0("Vectors differ by ", signif(best_diff, 4),
                   " (tolerance: ", tolerance, ")")
  )
}

#' Check that two vectors are highly correlated
expect_correlated <- function(actual, expected, min_cor = 0.95) {
  r <- cor(as.numeric(actual), as.numeric(expected))
  testthat::expect_true(
    r > min_cor,
    label = paste0("Correlation = ", round(r, 4), " (min: ", min_cor, ")")
  )
}

#' Compare sorted absolute values (for GFT with sign ambiguity)
expect_abs_sorted_equal <- function(actual, expected, tolerance = 1e-10) {
  testthat::expect_equal(
    sort(abs(as.numeric(actual))),
    sort(abs(as.numeric(expected))),
    tolerance = tolerance
  )
}

# ---------------------------------------------------------------------------
# Graph constructor mapping (name -> rgsp constructor)
# ---------------------------------------------------------------------------

golden_graph_constructors <- list(
  Ring_10         = function() graph_ring(10),
  Path_8          = function() graph_path(8),
  Grid2d_3_4      = function() graph_grid2d(3, 4),
  FullConnected_6 = function() graph_complete(6),
  Torus_4_5       = function() graph_grid2d(4, 5, periodic = TRUE)
)
