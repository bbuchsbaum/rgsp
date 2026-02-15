# Spectral Graph Wavelet Transform (SGWT)
#
# Provides clean forward/inverse transforms with proper frame handling,
# automatic scale selection, and batch processing of multiple signals.
#
# Reference: Hammond, Vandergheynst, Gribonval (2011).
# "Wavelets on Graphs via Spectral Graph Theory"
# Applied and Computational Harmonic Analysis, 30(2), 129-150.

#' Automatic scale selection for graph wavelets
#'
#' Generates logarithmically spaced scales covering the graph spectrum.
#'
#' @param g gsp_graph object
#' @param n_scales number of scales (default: 4)
#' @param lmax spectral radius (computed if NULL
#' @return numeric vector of scales
#' @keywords internal
auto_wavelet_scales <- function(g, n_scales = 4, lmax = NULL) {
  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)
  # Scales from fine (high freq) to coarse (low freq)
  # Mexican hat peak at lambda = sigma, so scale from lmax down to lmax/2^n_scales
  lmax / (2 ^ seq(0, n_scales - 1))
}

# Build wavelet filter bank (kernels + names) ---------------------------------

sgwt_filter_bank <- function(g, wavelet = "mexican_hat", scales = NULL,
                             include_lowpass = TRUE, lmax = NULL) {
  if (is.null(lmax)) {
    lmax <- lambda_max(g, normalized = g$normalized)
  }

  if (is.null(scales)) {
    scales <- auto_wavelet_scales(g, n_scales = 4, lmax = lmax)
  }
  scales <- sort(scales, decreasing = TRUE)  # coarse to fine

  kernel_generator <- if (is.function(wavelet)) {
    wavelet
  } else {
    switch(wavelet,
           mexican_hat = kernel_mexican_hat,
           meyer = function(s) kernel_meyer(lmax / (2 * s), lmax / s),
           heat = kernel_heat,
           stop("Unknown wavelet type: ", wavelet)
    )
  }

  wavelet_kernels <- lapply(scales, kernel_generator)

  if (include_lowpass) {
    lowpass_cutoff <- min(scales) / 2
    lowpass_kernel <- kernel_heat(1 / lowpass_cutoff)
    kernels <- c(list(lowpass = lowpass_kernel), wavelet_kernels)
    kernel_names <- c("lowpass", paste0("scale_", scales))
  } else {
    kernels <- wavelet_kernels
    kernel_names <- paste0("scale_", scales)
  }

  list(
    kernels = kernels,
    kernel_names = kernel_names,
    scales = scales,
    include_lowpass = include_lowpass,
    wavelet_name = if (is.function(wavelet)) "custom" else wavelet,
    lmax = lmax
  )
}

sgwt_cheby_coeffs <- function(kernels, K, lmax, jackson = FALSE) {
  lapply(kernels, function(kern) cheby_coeffs(kern, K = K, lmax = lmax, jackson = jackson))
}

.sgwt_apply_cheby <- function(L, x, cheby_list, lmax, band_names = NULL) {
  n_nodes <- nrow(L)
  n_signals <- ncol(x)
  n_bands <- length(cheby_list)

  coeffs <- array(NA_real_, dim = c(n_nodes, n_bands, n_signals))
  dimnames(coeffs) <- list(
    node = NULL,
    band = band_names,
    signal = if (n_signals > 1) paste0("sig_", seq_len(n_signals)) else NULL
  )

  for (b in seq_len(n_bands)) {
    coeffs[, b, ] <- chebyshev_filter_cpp(L, x, cheby_list[[b]], lmax)
  }
  coeffs
}

#' Spectral Graph Wavelet Transform (SGWT)
#'
#' Decompose graph signals into multi-scale wavelet coefficients using
#' spectral graph wavelets. Supports batch processing of multiple signals.
#'
#' @param g gsp_graph object (Laplacian will be computed if needed)
#' @param x numeric vector (single signal) or matrix (n_nodes x n_signals)
#' @param scales numeric vector of wavelet scales. If NULL, automatically
#'   selected to span the graph spectrum (default: NULL)
#' @param wavelet wavelet type: "mexican_hat" (default), "meyer", "heat",
#'   or a custom kernel-generating function that takes a scale and returns
#'   a kernel function
#' @param include_lowpass logical; include scaling function to capture
#'   low-frequency residual? Required for perfect reconstruction (default: TRUE)
#' @param K integer; Chebyshev polynomial order for approximation (default: 30)
#' @param lmax spectral radius; computed automatically if NULL
#'
#' @return Object of class "sgwt_coeffs" with components:
#'   \describe{
#'     \item{coeffs}{array with dimensions n_nodes x n_scales (+1 if lowpass) x n_signals}
#'     \item{scales}{scales used for decomposition}
#'     \item{wavelet}{wavelet type used}
#'     \item{include_lowpass}{whether lowpass band is included}
#'     \item{frame_bounds}{estimated frame bounds (A, B)}
#'     \item{lmax}{spectral radius used}
#'     \item{K}{Chebyshev order used}
#'     \item{g}{reference to the graph (for inverse transform)}
#'   }
#'
#' @details
#' The Spectral Graph Wavelet Transform (Hammond et al. 2011) decomposes signals
#' by "spatial frequency" on the graph topology. Fine scales capture local,
#' node-specific variation while coarse scales capture smooth, global patterns.
#'
#' The transform satisfies a frame condition: there exist constants A, B > 0 such that
#' \deqn{A \|x\|^2 \le \sum_s \|W_s x\|^2 \le B \|x\|^2}
#' When A = B (tight frame), perfect reconstruction is possible via simple summation.
#' When A != B, iterative reconstruction or pseudo-inverse is needed.
#'
#' @references
#' Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011). Wavelets on graphs
#' via spectral graph theory. Applied and Computational Harmonic Analysis, 30(2), 129-150.
#'
#' @export
#' @seealso [isgwt()], [filter_signal()], [kernel_mexican_hat()]
#'
#' @examples
#' # Basic usage
#' g <- graph_ring(30)
#' signal <- rnorm(30)
#' W <- sgwt(g, signal, scales = c(1, 2), K = 15)
#' print(W)
#'
sgwt <- function(g, x,
                               scales = NULL,
                               wavelet = "mexican_hat",
                               include_lowpass = TRUE,
                               K = 30,
                               lmax = NULL) {
  # Input validation
  if (!inherits(g, "gsp_graph")) {
    stop("'g' must be a gsp_graph object")
  }

  # Handle vector vs matrix input
  single_signal <- is.null(dim(x))
  if (single_signal) {
    x <- matrix(x, ncol = 1)
  }

  n_nodes <- nrow(x)
  n_signals <- ncol(x)

  validate_signal_dims(g, x)

  # Compute spectral radius if needed
  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)

  # Build wavelet bank + Chebyshev coefficients (graph-agnostic cache key)
  fb <- sgwt_filter_bank(g, wavelet = wavelet, scales = scales,
                         include_lowpass = include_lowpass, lmax = lmax)
  cheb_list <- sgwt_cheby_coeffs(fb$kernels, K = K, lmax = lmax)

  # Get Laplacian once
  L <- graph_laplacian(g, g$normalized)

  # Apply each kernel to all signals at once (vectorized over signals)
  coeffs <- .sgwt_apply_cheby(L, x, cheb_list, lmax, band_names = fb$kernel_names)

  # Estimate frame bounds (optional but useful for diagnostics)
  frame_bounds <- estimate_frame_bounds(fb$kernels, lmax, n_samples = 100)

  # Return structured object
  result <- structure(
    list(
      coeffs = if (single_signal) coeffs[, , 1, drop = FALSE] else coeffs,
      scales = fb$scales,
      wavelet = fb$wavelet_name,
      include_lowpass = fb$include_lowpass,
      frame_bounds = frame_bounds,
      lmax = lmax,
      K = K,
      g = g,
      n_signals = n_signals,
      single_signal = single_signal
    ),
    class = "sgwt_coeffs"
  )

  result
}

# Directional variants -------------------------------------------------------

#' Directional SGWT via graph gradient (Riesz-style)
#'
#' Computes directional wavelet coefficients by applying the graph gradient
#' to standard SGWT coefficients. This produces edge-based coefficients that
#' capture directional structure in the signal.
#'
#' @param g gsp_graph object
#' @param x numeric vector (single signal) or matrix (n_nodes x n_signals)
#' @param scales numeric vector of wavelet scales. If NULL, automatically selected
#' @param wavelet wavelet type: "heat" (default), "mexican_hat", "meyer"
#' @param include_lowpass logical; include lowpass band (default: TRUE)
#' @param K integer; Chebyshev polynomial order (default: 30)
#' @param lmax spectral radius; computed if NULL
#'
#' @return Object of class "dsgwt_riesz" with components:
#'   \describe{
#'     \item{node_coeffs}{array of node-space wavelet coefficients}
#'     \item{edge_coeffs}{array of edge-space (gradient) coefficients}
#'     \item{incidence}{graph incidence matrix}
#'     \item{gradient}{graph gradient operator}
#'     \item{g}{reference to the graph}
#'     \item{scales}{scales used}
#'     \item{wavelet}{wavelet type}
#'   }
#'
#' @export
#' @seealso [dsgwt_riesz_magnitude()], [dsgwt_steer()], [sgwt()]
dsgwt_riesz <- function(g, x,
                        scales = NULL,
                        wavelet = "heat",
                        include_lowpass = TRUE,
                        K = 30,
                        lmax = NULL) {
  W <- sgwt(g, x,
            scales = scales,
            wavelet = wavelet,
            include_lowpass = include_lowpass,
            K = K,
            lmax = lmax)

  D <- graph_gradient(g)
  B <- graph_incidence(g)$B

  node_coeffs <- W$coeffs
  n_nodes <- dim(node_coeffs)[1]
  n_bands <- dim(node_coeffs)[2]
  n_signals <- dim(node_coeffs)[3]

  coeffs_mat <- matrix(node_coeffs, nrow = n_nodes, ncol = n_bands * n_signals)
  edge_mat <- D %*% coeffs_mat
  edge_coeffs <- array(edge_mat, dim = c(nrow(D), n_bands, n_signals))

  structure(
    list(
      node_coeffs = node_coeffs,
      edge_coeffs = edge_coeffs,
      incidence = B,
      gradient = D,
      g = g,
      scales = W$scales,
      wavelet = W$wavelet,
      include_lowpass = W$include_lowpass,
      lmax = W$lmax,
      K = W$K,
      dim_order = c("node", "band", "signal"),
      n_signals = n_signals,
      single_signal = W$single_signal
    ),
    class = "dsgwt_riesz"
  )
}

#' Node-wise magnitude from directional (Riesz) coefficients
#'
#' Computes the magnitude of directional coefficients at each node by
#' aggregating edge coefficients via the incidence matrix.
#'
#' @param obj dsgwt_riesz object from [dsgwt_riesz()]
#'
#' @return array of node-wise magnitudes with dimensions (n_nodes x n_bands x n_signals)
#'
#' @export
#' @seealso [dsgwt_riesz()]
dsgwt_riesz_magnitude <- function(obj) {
  if (!inherits(obj, "dsgwt_riesz")) stop("obj must be a dsgwt_riesz")
  B <- obj$incidence
  edge_coeffs <- obj$edge_coeffs
  dimc <- dim(edge_coeffs)
  edge_mat <- matrix(edge_coeffs, nrow = dimc[1], ncol = dimc[2] * dimc[3])
  energy_edge <- edge_mat^2
  agg <- Matrix::t(abs(B))
  node_energy <- agg %*% energy_edge
  node_energy <- sqrt(node_energy)
  array(node_energy, dim = c(nrow(agg), dimc[2], dimc[3]))
}

#' Directional SGWT via steered Laplacians
#'
#' Computes directional wavelet coefficients using orientation-steered
#' graph Laplacians. Each direction produces a separate set of coefficients
#' capturing structure aligned with that orientation.
#'
#' @param g gsp_graph object with coordinates
#' @param x numeric vector (single signal) or matrix (n_nodes x n_signals)
#' @param n_directions integer; number of steering directions (default: 3)
#' @param orientations optional matrix of orientations (n_dir x d); if NULL,
#'   evenly spaced on the circle
#' @param scales numeric vector of wavelet scales. If NULL, automatically selected
#' @param wavelet wavelet type: "meyer" (default), "mexican_hat", "heat"
#' @param include_lowpass logical; include lowpass band (default: TRUE)
#' @param K integer; Chebyshev polynomial order (default: 30)
#' @param alpha_isotropic numeric in (0,1]; mixing with isotropic graph (default: 0.1)
#' @param p integer; power for directional weighting (default: 4)
#' @param absorb_sign logical; absorb sign of alignment (default: TRUE)
#' @param min_degree integer; minimum degree for steering (default: 3)
#' @param min_weight numeric; minimum edge weight (default: 1e-6)
#' @param lmax spectral radius; computed if NULL
#' @param check_lmax logical; warn if directional lmax deviates >10% (default: FALSE)
#' @param memory_saving logical; stream signals to reduce peak memory (default: FALSE)
#'
#' @return Object of class "dsgwt_steer" with components:
#'   \describe{
#'     \item{coeffs}{array of coefficients (n_nodes x n_bands x n_dir x n_signals)}
#'     \item{orientations}{matrix of orientation vectors}
#'     \item{n_dir}{number of directions}
#'     \item{scales}{scales used}
#'     \item{wavelet}{wavelet type}
#'     \item{g}{reference to the graph}
#'   }
#'
#' @export
#' @seealso [dsgwt_riesz()], [sgwt()], [graph_steerable_family()]
dsgwt_steer <- function(g, x,
                        n_directions = 3L,
                        orientations = NULL,
                        scales = NULL,
                        wavelet = "meyer",
                        include_lowpass = TRUE,
                        K = 30,
                        alpha_isotropic = 0.1,
                        p = 4L,
                        absorb_sign = TRUE,
                        min_degree = 3L,
                        min_weight = 1e-6,
                        lmax = NULL,
                        check_lmax = FALSE,
                        memory_saving = FALSE) {
  if (!inherits(g, "gsp_graph")) stop("'g' must be a gsp_graph object")

  single_signal <- is.null(dim(x))
  if (single_signal) x <- matrix(x, ncol = 1)
  validate_signal_dims(g, x)
  n_nodes <- nrow(x)
  n_signals <- ncol(x)

  lmax_base <- lmax %||% lambda_max(g, normalized = g$normalized)

  fb <- sgwt_filter_bank(g, wavelet = wavelet, scales = scales,
                         include_lowpass = include_lowpass, lmax = lmax_base)
  cheb_list <- sgwt_cheby_coeffs(fb$kernels, K = K, lmax = lmax_base)

  graphs <- graph_steerable_family(
    g,
    n_directions = n_directions,
    orientations = orientations,
    alpha_isotropic = alpha_isotropic,
    p = p,
    absorb_sign = absorb_sign,
    min_degree = min_degree,
    min_weight = min_weight,
    check_degrees = TRUE
  )

  if (check_lmax) {
    lmax_dirs <- vapply(graphs, function(gth) lambda_max(gth, normalized = gth$normalized), numeric(1))
    if (max(abs(lmax_dirs - lmax_base) / lmax_base) > 0.1) {
      warning("dsgwt_steer(): some orientations have lmax deviating >10% from base; consider adjusting alpha_isotropic or p")
    }
  }

  res_list <- lapply(graphs, function(gth) {
    Lth <- graph_laplacian(gth, gth$normalized)
    if (memory_saving) {
      # stream over signals to reduce peak memory
      out <- array(NA_real_, dim = c(n_nodes, length(fb$kernels), ncol(x)))
      for (s in seq_len(ncol(x))) {
        out[, , s] <- .sgwt_apply_cheby(Lth, x[, s, drop = FALSE], cheb_list, lmax_base,
                                        band_names = fb$kernel_names)[, , 1]
      }
      out
    } else {
      .sgwt_apply_cheby(Lth, x, cheb_list, lmax_base, band_names = fb$kernel_names)
    }
  })

  n_bands <- dim(res_list[[1]])[2]
  coeffs <- array(NA_real_, dim = c(n_nodes, n_bands, length(res_list), n_signals))
  dimnames(coeffs) <- list(
    node = NULL,
    band = fb$kernel_names,
    direction = paste0("dir_", seq_along(res_list)),
    signal = if (n_signals > 1) paste0("sig_", seq_len(n_signals)) else NULL
  )
  for (d in seq_along(res_list)) {
    coeffs[, , d, ] <- res_list[[d]]
  }

  orientation_mat <- do.call(rbind, lapply(graphs, function(gth) gth$orientation))

  structure(
    list(
      coeffs = if (single_signal) coeffs[, , , 1, drop = FALSE] else coeffs,
      orientations = orientation_mat,
      n_dir = length(res_list),
      scales = fb$scales,
      wavelet = fb$wavelet_name,
      include_lowpass = fb$include_lowpass,
      lmax = lmax_base,
      K = K,
      dim_order = c("node", "band", "direction", "signal"),
      g = g,
      n_signals = n_signals,
      single_signal = single_signal,
      alpha_isotropic = alpha_isotropic,
      p = p,
      absorb_sign = absorb_sign
    ),
    class = "dsgwt_steer"
  )
}

#' Directional wavelets on a k-NN graph (convenience)
#'
#' Build a k-NN graph from coordinates and apply directional graph wavelets
#' (steered Laplacians or Riesz/gradient) in one call.
#'
#' @param coords numeric matrix of node coordinates (n x d)
#' @param signal vector or matrix of signals (n x m)
#' @param method `"steer"` (orientation-steered Laplacians) or `"riesz"` (gradient-based)
#' @param k neighbors per node for k-NN construction
#' @param weight k-NN edge weighting scheme (`"heat"`, `"binary"`, `"distance"`)
#' @param sigma optional heat bandwidth for k-NN weights
#' @param sym k-NN symmetrization (`"union"` or `"mutual"`)
#' @param normalized logical; compute normalized Laplacian
#' @param seed optional RNG seed
#' @param ... passed to [dsgwt_steer()] or [dsgwt_riesz()] (e.g., `n_directions`, `orientations`, `scales`, `wavelet`, `K`)
#'
#' @return object of class `dsgwt_steer` or `dsgwt_riesz`
#' @export
directional_wavelet_knn <- function(coords, signal,
                                    method = c("steer", "riesz"),
                                    k = 6,
                                    weight = c("heat", "binary", "distance"),
                                    sigma = NULL,
                                    sym = c("union", "mutual"),
                                    normalized = FALSE,
                                    seed = NULL,
                                    ...) {
  method <- match.arg(method)
  g <- graph_knn(coords,
                 k = k,
                 weight = weight,
                 sigma = sigma,
                 sym = sym,
                 normalized = normalized,
                 seed = seed)
  if (identical(method, "steer")) {
    dsgwt_steer(g, signal, ...)
  } else {
    dsgwt_riesz(g, signal, ...)
  }
}

#' Estimate Frame Bounds for a Filter Bank
#'
#' Computes the frame bounds A and B by evaluating sum of squared filter
#' responses across the spectrum.
#'
#' @param kernels list of kernel functions
#' @param lmax spectral radius
#' @param n_samples number of spectral samples
#' @return named numeric vector c(A = lower, B = upper)
#' @keywords internal
estimate_frame_bounds <- function(kernels, lmax, n_samples = 100) {
  lambdas <- seq(0, lmax, length.out = n_samples)

  # Sum of squared responses at each lambda
  frame_response <- numeric(n_samples)
  for (kern in kernels) {
    frame_response <- frame_response + kern(lambdas)^2
  }

  c(A = min(frame_response), B = max(frame_response))
}

#' Compute frame diagonal via Hutchinson or exact eigendecomposition
#'
#' @param g graph
#' @param kernels list of kernel functions
#' @param method "diag_hutch" (default) or "exact"
#' @param R Hutchinson probe count (diag_hutch only)
#' @param K Chebyshev order (used by filtering)
#' @param lmax optional spectral radius
#' @param jackson use Jackson damping for Chebyshev filters
#' @return list with elements diag (length n), A, B, method
#' @export
compute_frame <- function(g, kernels, method = c("diag_hutch", "exact"),
                          R = 32, K = 30, lmax = NULL, jackson = FALSE) {
  method <- match.arg(method)
  n <- g$n
  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)

  if (method == "exact") {
    eig <- graph_eigenpairs(g)
    lam <- eig$values
    weights <- Reduce("+", lapply(kernels, function(k) (k(lam))^2))
    diag_frame <- as.numeric((eig$vectors^2) %*% weights)
  } else {
    L <- graph_laplacian(g, normalized = g$normalized)
    cheb_list <- lapply(kernels, function(kern) {
      cheby_coeffs(function(lambda) kern(lambda)^2, K = K, lmax = lmax, jackson = jackson)
    })
    diag_accum <- numeric(n)
    for (r in seq_len(R)) {
      z <- sample(c(-1, 1), n, replace = TRUE)
      Sz <- Reduce(`+`, Map(function(coeffs) {
        chebyshev_filter_cpp(L, matrix(z, ncol = 1), coeffs, lmax)
      }, cheb_list))
      Sz_vec <- drop(Sz)
      diag_accum <- diag_accum + z * Sz_vec
    }
    diag_frame <- diag_accum / R
  }
  list(
    diag = diag_frame,
    A = min(diag_frame),
    B = max(diag_frame),
    method = method
  )
}

#' Compute graph spectrogram (vertex-frequency energy)
#'
#' Applies a bank of window kernels and returns per-node energy at each scale.
#'
#' @param g graph
#' @param signal numeric vector/matrix (n or n x m)
#' @param scales numeric vector of window scales
#' @param window function taking a scale and returning a kernel; default heat
#' @param K Chebyshev order
#' @param lmax optional spectral radius
#' @param jackson logical; use Jackson damping
#' @return matrix of size n x length(scales) (or array n x scales x m)
#' @export
compute_spectrogram <- function(g, signal, scales,
                                window = function(s) kernel_heat(s),
                                K = 30, lmax = NULL, jackson = FALSE) {
  if (is.null(dim(signal))) signal <- matrix(signal, ncol = 1)
  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)
  out_list <- lapply(scales, function(s) {
    kern <- window(s)
    filter_signal(g, signal, kern, K = K, lmax = lmax, jackson = jackson)
  })
  arr <- array(unlist(out_list), dim = c(nrow(signal), length(scales), ncol(signal)))
  energy <- arr^2
  if (ncol(signal) == 1) {
    energy <- matrix(energy, nrow = nrow(signal), ncol = length(scales))
    colnames(energy) <- paste0("s_", scales)
  }
  energy
}


#' Inverse Spectral Graph Wavelet Transform
#'
#' Reconstruct graph signals from SGWT coefficients. Supports scale-selective
#' reconstruction by specifying which bands to include or exclude.
#'
#' @param W object of class "sgwt_coeffs" from [sgwt()]
#' @param bands which bands to use for reconstruction. Options:
#'   - NULL (default): use all bands
#'   - numeric vector: scale values to include (e.g., c(5, 10))
#'   - character vector: band names (e.g., c("lowpass", "scale_5"))
#'   - "lowpass": only the lowpass band
#'   - "wavelets": all wavelet bands (exclude lowpass)
#' @param exclude_bands bands to exclude (same format as bands). Ignored if
#'   bands is specified.
#' @param method reconstruction method:
#'   - "auto" (default): choose based on frame bounds
#'   - "sum": direct summation (exact for tight frames)
#'   - "iterative": conjugate gradient (for non-tight frames)
#' @param tol tolerance for iterative methods (default: 1e-6)
#' @param maxiter maximum iterations for iterative methods (default: 100)
#'
#' @return numeric vector or matrix of reconstructed signals (same shape as input to sgwt)
#'
#' @details
#' For tight frames (A = B), reconstruction is simply:
#' \deqn{x = (1/A) \sum_s W_s^* W_s x}
#'
#' For non-tight frames, iterative methods solve the normal equations.
#'
#' When using scale-selective reconstruction (bands or exclude_bands), the
#' result is a partial reconstruction containing only the selected frequency
#' components. This is useful for analyzing scale-specific contributions.
#'
#' @references
#' Hammond, D. K., Vandergheynst, P., & Gribonval, R. (2011). Wavelets on graphs
#' via spectral graph theory. Applied and Computational Harmonic Analysis, 30(2), 129-150.
#'
#' @export
#' @seealso [sgwt()], [sgwt_node_contribution()]
#'
#' @examples
#' g <- graph_sensor(100, seed = 42)
#' signal <- rnorm(100)
#' W <- sgwt(g, signal, scales = c(2, 5, 10, 20))
#'
#' # Full reconstruction
#' rec_full <- isgwt(W)
#'
#' # Reconstruct using only coarse scales (smooth component)
#' rec_coarse <- isgwt(W, bands = c(10, 20))
#'
#' # Reconstruct using only fine scales (detail component)
#' rec_fine <- isgwt(W, bands = c(2, 5))
#'
#' # Exclude lowpass (wavelets only)
#' rec_wavelets <- isgwt(W, bands = "wavelets")
#'
#' # Only lowpass
#' rec_lowpass <- isgwt(W, bands = "lowpass")
#'
isgwt <- function(W,
                  bands = NULL,
                  exclude_bands = NULL,
                  method = "auto",
                  tol = 1e-6,
                  maxiter = 100) {
  if (!inherits(W, "sgwt_coeffs")) {
    stop("'W' must be an sgwt_coeffs object from sgwt()")
  }

 # Extract components
  coeffs <- W$coeffs
  g <- W$g
  scales <- W$scales
  lmax <- W$lmax
  K <- W$K
  frame_bounds <- W$frame_bounds
  include_lowpass <- W$include_lowpass
  wavelet <- W$wavelet

  # Dimensions
  dims <- dim(coeffs)
  n_nodes <- dims[1]
  n_bands <- dims[2]
  n_signals <- if (length(dims) == 3) dims[3] else 1

  # Build band names
  if (include_lowpass) {
    band_names <- c("lowpass", paste0("scale_", scales))
  } else {
    band_names <- paste0("scale_", scales)
  }

  # Determine which bands to use
  use_bands <- seq_len(n_bands)

  if (!is.null(bands)) {
    if (is.character(bands) && length(bands) == 1) {
      if (bands == "lowpass") {
        if (!include_lowpass) stop("No lowpass band in this transform")
        use_bands <- 1
      } else if (bands == "wavelets") {
        use_bands <- if (include_lowpass) 2:n_bands else seq_len(n_bands)
      } else {
        # Single band name
        use_bands <- match(bands, band_names)
        if (any(is.na(use_bands))) {
          stop("Band '", bands[is.na(use_bands)][1], "' not found. Available: ",
               paste(band_names, collapse = ", "))
        }
      }
    } else if (is.character(bands)) {
      # Multiple band names
      use_bands <- match(bands, band_names)
      if (any(is.na(use_bands))) {
        stop("Band '", bands[is.na(use_bands)][1], "' not found. Available: ",
             paste(band_names, collapse = ", "))
      }
    } else if (is.numeric(bands)) {
      # Scale values
      use_bands <- match(paste0("scale_", bands), band_names)
      if (any(is.na(use_bands))) {
        stop("Scale ", bands[is.na(use_bands)][1], " not found. Available: ",
             paste(scales, collapse = ", "))
      }
    }
  } else if (!is.null(exclude_bands)) {
    if (is.character(exclude_bands) && length(exclude_bands) == 1 && exclude_bands == "lowpass") {
      exclude_idx <- 1
    } else if (is.character(exclude_bands)) {
      exclude_idx <- match(exclude_bands, band_names)
    } else if (is.numeric(exclude_bands)) {
      exclude_idx <- match(paste0("scale_", exclude_bands), band_names)
    }
    exclude_idx <- exclude_idx[!is.na(exclude_idx)]
    use_bands <- setdiff(seq_len(n_bands), exclude_idx)
  }

  if (length(use_bands) == 0) {
    stop("No bands selected for reconstruction")
  }

  # Rebuild kernels (same as forward transform)
  kernel_generator <- switch(wavelet,
    mexican_hat = kernel_mexican_hat,
    meyer = function(s) kernel_meyer(lmax / (2 * s), lmax / s),
    heat = kernel_heat,
    custom = stop("Cannot invert custom wavelet without stored kernels")
  )

  wavelet_kernels <- lapply(scales, kernel_generator)

  if (include_lowpass) {
    lowpass_cutoff <- min(scales) / 2
    lowpass_kernel <- kernel_heat(1 / lowpass_cutoff)
    all_kernels <- c(list(lowpass_kernel), wavelet_kernels)
  } else {
    all_kernels <- wavelet_kernels
  }

  # For scale-selective reconstruction, use "sum" method with no normalization
  # (partial reconstruction doesn't satisfy frame equations)
  partial_reconstruction <- length(use_bands) < n_bands

  # Choose method based on frame bounds (only for full reconstruction)
  if (method == "auto") {
    if (partial_reconstruction) {
      method <- "sum"
    } else {
      ratio <- frame_bounds["B"] / frame_bounds["A"]
      method <- if (ratio < 1.1 && frame_bounds["A"] > 0.5) "sum" else "iterative"
    }
  }

  # Get Laplacian
  L <- graph_laplacian(g, g$normalized)
  # Precompute Chebyshev coefficients for kernels and squared kernels
  cheb_all <- lapply(all_kernels, function(kern) cheby_coeffs(kern, K = K, lmax = lmax))
  cheb_sq_all <- lapply(all_kernels, function(kern) {
    cheby_coeffs(function(lambda) kern(lambda)^2, K = K, lmax = lmax)
  })

  # Reconstruction (only over selected bands)
  if (method == "sum") {
    # Direct summation with frame normalization
    # x_rec = (1/A) * sum_b W_b^T W_b x = (1/A) * sum_b W_b coeffs_b
    # Since W_b is self-adjoint (symmetric kernel), W_b^T = W_b

    A <- frame_bounds["A"]
    if (A < 1e-10 || partial_reconstruction) {
      # Don't normalize for partial reconstruction
      A <- 1
    }

    x_rec <- matrix(0, nrow = n_nodes, ncol = n_signals)

    for (b in use_bands) {
      cheby_c <- cheb_all[[b]]

      # Apply kernel to wavelet coefficients (adjoint = same kernel for real symmetric)
      if (n_signals == 1) {
        band_coeffs <- matrix(coeffs[, b, ], ncol = 1)
      } else {
        band_coeffs <- coeffs[, b, ]
      }
      x_rec <- x_rec + chebyshev_filter_cpp(L, band_coeffs, cheby_c, lmax)
    }

    if (!partial_reconstruction) {
      x_rec <- x_rec / A
    }

  } else if (method == "iterative") {
    # Conjugate gradient on normal equations
    # Solve: (sum_b W_b^T W_b) x = sum_b W_b^T coeffs_b

    x_rec <- matrix(0, nrow = n_nodes, ncol = n_signals)

    # Right-hand side: sum_b W_b^T coeffs_b (over selected bands)
    rhs <- matrix(0, nrow = n_nodes, ncol = n_signals)
    for (b in use_bands) {
      cheby_c <- cheb_all[[b]]
      if (n_signals == 1) {
        band_coeffs <- matrix(coeffs[, b, ], ncol = 1)
      } else {
        band_coeffs <- coeffs[, b, ]
      }
      rhs <- rhs + chebyshev_filter_cpp(L, band_coeffs, cheby_c, lmax)
    }

    # CG iteration for each signal
    for (s in seq_len(n_signals)) {
      x <- numeric(n_nodes)
      r <- rhs[, s]
      p <- r
      rsold <- sum(r^2)

      for (iter in seq_len(maxiter)) {
        # Apply frame operator: Ap = sum_b W_b^T W_b p = sum_b W_b^2 p (over selected bands)
        Ap <- numeric(n_nodes)
        for (b in use_bands) {
          cheby_c <- cheb_sq_all[[b]]
          Ap <- Ap + chebyshev_filter_cpp(L, matrix(p, ncol = 1), cheby_c, lmax)
        }

        alpha <- rsold / sum(p * Ap)
        x <- x + alpha * p
        r <- r - alpha * Ap
        rsnew <- sum(r^2)

        if (sqrt(rsnew) < tol) break

        p <- r + (rsnew / rsold) * p
        rsold <- rsnew
      }

      x_rec[, s] <- x
    }

  } else {
    stop("Unknown method: ", method)
  }

  # Return in same shape as input
  if (W$single_signal) {
    drop(x_rec)
  } else {
    x_rec
  }
}


# S3 Methods ----------------------------------------------------------------

#' Print SGWT Coefficients
#' @param x sgwt_coeffs object
#' @param ... ignored
#' @method print sgwt_coeffs
#' @export
print.sgwt_coeffs <- function(x, ...) {
  dims <- dim(x$coeffs)
  n_nodes <- dims[1]
  n_bands <- dims[2]
  n_signals <- if (length(dims) == 3) dims[3] else 1

  cat("Spectral Graph Wavelet Transform (SGWT) Coefficients\n")
  cat("====================================================\n")
  cat("  Nodes:", n_nodes, "\n")
  cat("  Scales:", paste(round(x$scales, 2), collapse = ", "), "\n")
  cat("  Bands:", n_bands,
      if (x$include_lowpass) "(including lowpass)" else "(no lowpass)", "\n")
  cat("  Signals:", n_signals, "\n")
  cat("  Wavelet:", x$wavelet, "\n")
  cat("  Frame bounds: A =", round(x$frame_bounds["A"], 3),
      ", B =", round(x$frame_bounds["B"], 3), "\n")

  ratio <- x$frame_bounds["B"] / x$frame_bounds["A"]
  if (ratio < 1.1) {
    cat("  Frame: tight (ratio =", round(ratio, 3), ")\n")
  } else {
    cat("  Frame: non-tight (ratio =", round(ratio, 3),
        ") - use iterative reconstruction\n")
  }

  invisible(x)
}

#' Print method for dsgwt_riesz
#' @param x dsgwt_riesz object
#' @param ... additional arguments (ignored)
#' @return invisibly returns the object
#' @method print dsgwt_riesz
#' @export
print.dsgwt_riesz <- function(x, ...) {
  dims_node <- dim(x$node_coeffs)
  dims_edge <- dim(x$edge_coeffs)
  cat("Directional SGWT (Riesz / gradient)\n")
  cat("===================================\n")
  cat("  Nodes:", dims_node[1], "Edges:", dims_edge[1], "\n")
  cat("  Bands:", dims_node[2],
      if (x$include_lowpass) "(including lowpass)" else "(no lowpass)", "\n")
  cat("  Signals:", dims_node[3], "\n")
  cat("  Wavelet:", x$wavelet, "\n")
  invisible(x)
}

#' Print method for dsgwt_steer
#' @param x dsgwt_steer object
#' @param ... additional arguments (ignored)
#' @return invisibly returns the object
#' @method print dsgwt_steer
#' @export
print.dsgwt_steer <- function(x, ...) {
  dims <- dim(x$coeffs)
  cat("Directional SGWT (steered Laplacians)\n")
  cat("=====================================\n")
  cat("  Nodes:", dims[1], "Bands:", dims[2],
      "Directions:", dims[3], "Signals:", dims[4], "\n")
  cat("  Wavelet:", x$wavelet,
      "  alpha_isotropic:", x$alpha_isotropic,
      "  p:", x$p, "\n")
  invisible(x)
}


#' Summarize Wavelet Coefficients
#' @param object sgwt_coeffs object
#' @param ... ignored
#' @method summary sgwt_coeffs
#' @export
summary.sgwt_coeffs <- function(object, ...) {
  dims <- dim(object$coeffs)
  n_bands <- dims[2]

 # Energy per band - handle 2D vs 3D arrays
  if (length(dims) == 3) {
    # [nodes, bands, signals] -> collapse over nodes and signals
    band_energy <- apply(object$coeffs, 2, function(b) mean(b^2))
  } else {
    # [nodes, bands] -> collapse over nodes
    band_energy <- colMeans(object$coeffs^2)
  }

  band_names <- if (object$include_lowpass) {
    c("lowpass", paste0("scale_", round(object$scales, 2)))
  } else {
    paste0("scale_", round(object$scales, 2))
  }

  cat("Wavelet Energy Distribution\n")
  cat("===========================\n")
  total_energy <- sum(band_energy)
  for (i in seq_len(n_bands)) {
    pct <- 100 * band_energy[i] / total_energy
    bar_len <- max(1, round(pct / 2))
    bar <- paste(rep("#", bar_len), collapse = "")
    cat(sprintf("  %-15s: %6.1f%% %s\n", band_names[i], pct, bar))
  }

  invisible(list(band_energy = setNames(band_energy, band_names)))
}


#' Plot Wavelet Coefficients
#'
#' Visualize the wavelet decomposition across scales.
#'
#' @param x sgwt_coeffs object
#' @param signal_idx which signal to plot (default: 1)
#' @param type "heatmap" (default), "energy", or "coeffs"
#' @param ... additional arguments passed to plotting functions
#' @method plot sgwt_coeffs
#' @export
plot.sgwt_coeffs <- function(x, signal_idx = 1, type = "energy", ...) {
  dims <- dim(x$coeffs)
  n_bands <- dims[2]

  band_names <- if (x$include_lowpass) {
    c("lowpass", paste0("s=", round(x$scales, 1)))
  } else {
    paste0("s=", round(x$scales, 1))
  }

  if (type == "energy") {
    # Bar plot of energy per band
    if (length(dims) == 3) {
      coeffs_sig <- x$coeffs[, , signal_idx]
    } else {
      coeffs_sig <- x$coeffs
    }

    band_energy <- apply(coeffs_sig, 2, function(b) sum(b^2))

    barplot(band_energy,
            names.arg = band_names,
            main = "Wavelet Energy by Scale",
            xlab = "Scale",
            ylab = "Energy",
            col = "steelblue",
            ...)

  } else if (type == "coeffs") {
    # Coefficient magnitudes per band
    if (length(dims) == 3) {
      coeffs_sig <- x$coeffs[, , signal_idx]
    } else {
      coeffs_sig <- x$coeffs
    }

    old_par <- par(mfrow = c(1, n_bands), mar = c(2, 2, 2, 1))
    on.exit(par(old_par))

    for (b in seq_len(n_bands)) {
      hist(coeffs_sig[, b],
           main = band_names[b],
           xlab = "", ylab = "",
           col = "steelblue",
           border = "white")
    }

  } else if (type == "heatmap") {
    # Nodes x bands heatmap
    if (length(dims) == 3) {
      coeffs_sig <- x$coeffs[, , signal_idx]
    } else {
      coeffs_sig <- x$coeffs
    }

    image(t(coeffs_sig),
          main = "Wavelet Coefficients (nodes x scales)",
          xlab = "Scale", ylab = "Node",
          col = hcl.colors(50, "RdBu", rev = TRUE),
          axes = FALSE,
          ...)
    axis(1, at = seq(0, 1, length.out = n_bands), labels = band_names)
  }

  invisible(x)
}


#' Extract SGWT Coefficients at a Specific Scale
#'
#' @param W sgwt_coeffs object from [sgwt()]
#' @param scale numeric scale value or "lowpass"
#' @return matrix of coefficients (n_nodes x n_signals)
#' @export
sgwt_coeffs_at_scale <- function(W, scale) {
  if (scale == "lowpass") {
    if (!W$include_lowpass) {
      stop("No lowpass band in this transform")
    }
    band_idx <- 1
  } else {
    band_idx <- which(W$scales == scale)
    if (length(band_idx) == 0) {
      stop("Scale ", scale, " not found. Available: ",
           paste(W$scales, collapse = ", "))
    }
    if (W$include_lowpass) band_idx <- band_idx + 1
  }

  if (length(dim(W$coeffs)) == 3) {
    W$coeffs[, band_idx, ]
  } else {
    W$coeffs[, band_idx]
  }
}


#' Compute SGWT Wavelet Atom at a Node
#'
#' Computes the wavelet atom (localized basis function) centered at a specific
#' node. This shows what the wavelet "sees" from that node's perspective —
#' its spatial extent and how it weights neighboring nodes.
#'
#' @param g gsp_graph object
#' @param node integer; node index (1-based)
#' @param scale numeric; wavelet scale
#' @param wavelet wavelet type: "mexican_hat" (default), "meyer", "heat"
#' @param K integer; Chebyshev order (default: 30)
#' @param lmax spectral radius; computed if NULL
#'
#' @return numeric vector of length n (the atom values at each node)
#'
#' @details
#' The wavelet atom at node i and scale s is defined as:
#' \deqn{\psi_{s,i} = W_s \delta_i}
#' where \eqn{\delta_i} is the Kronecker delta at node i.
#'
#' Atoms are localized around the center node, with spatial extent determined
#' by the scale parameter. Larger scales produce more spread-out atoms.
#'
#' @export
#' @seealso [sgwt()], [sgwt_node_contribution()]
#'
#' @examples
#' g <- graph_sensor(100, seed = 42)
#'
#' # Wavelet atom at node 42, scale 5
#' atom <- sgwt_atom(g, node = 42, scale = 5)
#'
#' # Compare atoms at different scales
#' atom_fine <- sgwt_atom(g, node = 42, scale = 2)
#' atom_coarse <- sgwt_atom(g, node = 42, scale = 10)
#'
#' # Atoms are localized around the center node
#' which.max(abs(atom))  # should be 42 or nearby
#'
sgwt_atom <- function(g, node, scale,
                      wavelet = "mexican_hat",
                      K = 30,
                      lmax = NULL) {
  if (!inherits(g, "gsp_graph")) {
    stop("'g' must be a gsp_graph object")
  }

  n <- g$n
  if (node < 1 || node > n) {
    stop("node must be between 1 and ", n)
  }

  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)

  # Build kernel for this scale
  kernel_generator <- switch(wavelet,
    mexican_hat = kernel_mexican_hat,
    meyer = function(s) kernel_meyer(lmax / (2 * s), lmax / s),
    heat = kernel_heat,
    stop("Unknown wavelet type: ", wavelet)
  )

  kern <- kernel_generator(scale)

  # Delta signal at node
  delta <- rep(0, n)
  delta[node] <- 1

  # Apply wavelet filter to delta
  atom <- filter_signal(g, delta, kern, K = K, lmax = lmax)

  attr(atom, "node") <- node
  attr(atom, "scale") <- scale
  attr(atom, "wavelet") <- wavelet

 atom
}


#' Compute Per-Node Scale Contributions
#'
#' Analyzes how much each wavelet scale contributes to the signal at each node.
#' Returns a matrix showing the energy or coefficient magnitude from each
#' scale at each node.
#'
#' @param W sgwt_coeffs object from [sgwt()]
#' @param type what to compute:
#'   - "energy" (default): squared coefficients (energy contribution)
#'   - "abs": absolute coefficient values
#'   - "raw": raw coefficients (can be negative)
#' @param normalize logical; if TRUE, normalize so each node's contributions
#'   sum to 1 (shows relative scale importance per node)
#'
#' @return matrix of size (n_nodes x n_bands) with scale contributions per node
#'
#' @details
#' This function answers: "For node i, how much of its signal comes from
#' each wavelet scale?"
#'
#' With `normalize = TRUE`, each row sums to 1, showing the relative importance
#' of each scale at that node. This is useful for identifying which nodes are
#' dominated by fine-scale (local) vs coarse-scale (global) structure.
#'
#' @export
#' @seealso [sgwt()], [sgwt_atom()], [isgwt()]
#'
#' @examples
#' g <- graph_sensor(100, seed = 42)
#' signal <- rnorm(100)
#' W <- sgwt(g, signal, scales = c(2, 5, 10, 20))
#'
#' # Energy contribution per node per scale
#' contrib <- sgwt_node_contribution(W)
#' dim(contrib)  # [100, 5] - 100 nodes, 5 bands
#'
#' # Which scale dominates each node?
#' contrib_norm <- sgwt_node_contribution(W, normalize = TRUE)
#' dominant_scale <- apply(contrib_norm, 1, which.max)
#'
#' # Find nodes dominated by fine scales (local structure)
#' fine_scale_nodes <- which(dominant_scale >= 4)
#'
#' # Find nodes dominated by coarse scales (global structure)
#' coarse_scale_nodes <- which(dominant_scale <= 2)
#'
sgwt_node_contribution <- function(W, type = "energy", normalize = FALSE) {
  if (!inherits(W, "sgwt_coeffs")) {
    stop("'W' must be an sgwt_coeffs object from sgwt()")
  }

  coeffs <- W$coeffs
  dims <- dim(coeffs)
  n_nodes <- dims[1]
  n_bands <- dims[2]
  n_signals <- if (length(dims) == 3) dims[3] else 1

  # Build band names
  if (W$include_lowpass) {
    band_names <- c("lowpass", paste0("scale_", W$scales))
  } else {
    band_names <- paste0("scale_", W$scales)
  }

  # Compute contribution based on type
  if (n_signals == 1) {
    # Single signal: [n_nodes, n_bands]
    contrib <- switch(type,
      energy = coeffs[, , 1]^2,
      abs = abs(coeffs[, , 1]),
      raw = coeffs[, , 1],
      stop("Unknown type: ", type)
    )
  } else {
    # Multiple signals: average over signals
    contrib <- switch(type,
      energy = apply(coeffs^2, c(1, 2), mean),
      abs = apply(abs(coeffs), c(1, 2), mean),
      raw = apply(coeffs, c(1, 2), mean),
      stop("Unknown type: ", type)
    )
  }

  # Ensure 2D matrix
  contrib <- matrix(contrib, nrow = n_nodes, ncol = n_bands)
  colnames(contrib) <- band_names

  # Normalize if requested
 if (normalize) {
    row_sums <- rowSums(contrib)
    row_sums[row_sums == 0] <- 1  # avoid division by zero
    contrib <- contrib / row_sums
  }

  contrib
}


#' Reconstruct Signal at Specific Nodes
#'
#' Convenience function to reconstruct and extract signal values at specific
#' nodes, optionally using only certain scales.
#'
#' @param W sgwt_coeffs object from [sgwt()]
#' @param nodes integer vector; node indices to extract (default: all)
#' @param bands bands to use for reconstruction (passed to [isgwt()])
#'
#' @return numeric vector (if single signal and single node) or matrix
#'
#' @export
#' @seealso [isgwt()], [sgwt_node_contribution()]
#'
#' @examples
#' g <- graph_sensor(100, seed = 42)
#' signal <- rnorm(100)
#' W <- sgwt(g, signal, scales = c(2, 5, 10, 20))
#'
#' # Get reconstructed value at node 42
#' val <- sgwt_reconstruct_at(W, nodes = 42)
#'
#' # Get coarse-scale component at nodes 1-10
#' coarse_vals <- sgwt_reconstruct_at(W, nodes = 1:10, bands = c(10, 20))
#'
sgwt_reconstruct_at <- function(W, nodes = NULL, bands = NULL) {
  if (!inherits(W, "sgwt_coeffs")) {
    stop("'W' must be an sgwt_coeffs object from sgwt()")
  }

  # Full reconstruction (or partial if bands specified)
  rec <- isgwt(W, bands = bands)

  # Extract nodes
  if (is.null(nodes)) {
    return(rec)
  }

  if (any(nodes < 1) || any(nodes > W$g$n)) {
    stop("nodes must be between 1 and ", W$g$n)
  }

  if (is.matrix(rec)) {
    rec[nodes, , drop = FALSE]
  } else {
    rec[nodes]
  }
}


#' Wavelet-space statistic mapped back to nodes
#'
#' Runs the SGWT, aggregates coefficients across signals with a statistic
#' (e.g., mean/median), and reconstructs the aggregated result to the node
#' domain. Useful for quickly visualizing group-level patterns.
#'
#' @param g gsp_graph object
#' @param X numeric vector or matrix (n_nodes x n_signals)
#' @param stat function applied across signals for each coefficient
#'   (default: [mean])
#' @param bands optional bands to use in the reconstruction (passed to [isgwt()])
#' @param wavelet wavelet type (passed to [sgwt()])
#' @param include_lowpass logical; include scaling function (default: TRUE)
#' @param K Chebyshev order (default: 30)
#' @param lmax spectral radius (computed if NULL)
#' @param normalize logical; if TRUE, divide aggregated coeffs by their max
#'   absolute value before reconstruction
#' @param ... extra arguments passed to `stat`
#'
#' @return list with elements:
#'   \describe{
#'     \item{map}{node-level reconstruction of the aggregated statistic}
#'     \item{coeffs}{aggregated coefficients (n_nodes x n_bands)}
#'     \item{sgwt}{`sgwt_coeffs` object with aggregated coeffs}
#'   }
#' @export
#'
#' @examples
#' g <- graph_ring(20)
#' X <- matrix(rnorm(20 * 5), nrow = 20)  # 5 signals on 20 nodes
#' res <- wavelet_stat_map(g, X, stat = mean)
#' plot(res$map, type = "h")
#'
wavelet_stat_map <- function(g, X, stat = mean, bands = NULL,
                             wavelet = "mexican_hat",
                             include_lowpass = TRUE,
                             K = 30, lmax = NULL,
                             normalize = FALSE, ...) {
  W <- sgwt(g, X,
            wavelet = wavelet,
            include_lowpass = include_lowpass,
            K = K, lmax = lmax)

  coeffs <- W$coeffs
  dims <- dim(coeffs)

  # Aggregate across signals (if multiple)
  if (length(dims) == 3) {
    agg_coeffs <- apply(coeffs, c(1, 2), stat, ...)
  } else {
    agg_coeffs <- coeffs
  }

  if (normalize) {
    max_abs <- max(abs(agg_coeffs))
    if (max_abs > 0) {
      agg_coeffs <- agg_coeffs / max_abs
    }
  }

  # Ensure coeffs retain 3D shape expected by isgwt
  if (length(dim(agg_coeffs)) == 2) {
    agg_coeffs <- array(agg_coeffs,
                        dim = c(nrow(agg_coeffs), ncol(agg_coeffs), 1))
  }

  # Reuse metadata but mark as single-signal
  W$coeffs <- agg_coeffs
  W$n_signals <- 1
  W$single_signal <- TRUE

  map <- isgwt(W, bands = bands)

  list(
    map = map,
    coeffs = agg_coeffs,
    sgwt = W
  )
}


#' Tidy wavelet coefficients
#'
#' Returns coefficients as a long or wide data frame for downstream
#' data manipulation and plotting.
#'
#' @param W sgwt_coeffs object from [sgwt()]
#' @param long logical; return long format (default TRUE). Wide format
#'   is available only for single-signal coefficients.
#'
#' @return data.frame
#' @export
#'
#' @examples
#' g <- graph_ring(10)
#' X <- matrix(rnorm(10 * 3), nrow = 10)
#' W <- sgwt(g, X, scales = c(2, 4))
#' df_long <- sgwt_coeffs_tidy(W)
#'
sgwt_coeffs_tidy <- function(W, long = TRUE) {
  if (!inherits(W, "sgwt_coeffs")) {
    stop("'W' must be an sgwt_coeffs object from sgwt()")
  }

  coeffs <- W$coeffs
  dims <- dim(coeffs)

  if (length(dims) == 2) {
    n_nodes <- dims[1]
    n_bands <- dims[2]
    n_signals <- 1
    arr <- array(coeffs, dim = c(n_nodes, n_bands, 1))
  } else {
    n_nodes <- dims[1]
    n_bands <- dims[2]
    n_signals <- dims[3]
    arr <- coeffs
  }

  band_names <- if (W$include_lowpass) {
    c("lowpass", paste0("scale_", W$scales))
  } else {
    paste0("scale_", W$scales)
  }

  scale_vals <- if (W$include_lowpass) {
    c(NA, W$scales)
  } else {
    W$scales
  }

  if (long) {
    data.frame(
      node = rep(seq_len(n_nodes), times = n_bands * n_signals),
      band = rep(band_names, each = n_nodes, times = n_signals),
      scale = rep(scale_vals, each = n_nodes, times = n_signals),
      signal = rep(seq_len(n_signals), each = n_nodes * n_bands),
      coeff = as.numeric(arr),
      stringsAsFactors = FALSE
    )
  } else {
    if (n_signals != 1) {
      stop("Wide format only supported when there is a single signal")
    }
    out <- data.frame(node = seq_len(n_nodes), coeffs, check.names = FALSE)
    colnames(out) <- c("node", band_names)
    out
  }
}


#' Band energy summary (data-frame friendly)
#'
#' Computes average energy per band and returns a tidy data frame
#' with percentages for quick plotting.
#'
#' @param W sgwt_coeffs object from [sgwt()]
#'
#' @return data.frame with columns band, scale, energy, pct
#' @export
#'
#' @examples
#' g <- graph_ring(12)
#' x <- rnorm(12)
#' W <- sgwt(g, x, scales = c(1, 2, 4))
#' band_energy_df(W)
#'
band_energy_df <- function(W) {
  if (!inherits(W, "sgwt_coeffs")) {
    stop("'W' must be an sgwt_coeffs object from sgwt()")
  }

  coeffs <- W$coeffs
  dims <- dim(coeffs)

  if (length(dims) == 3) {
    band_energy <- apply(coeffs^2, 2, mean)
  } else {
    band_energy <- colMeans(coeffs^2)
  }

  band_names <- if (W$include_lowpass) {
    c("lowpass", paste0("scale_", W$scales))
  } else {
    paste0("scale_", W$scales)
  }

  scale_vals <- if (W$include_lowpass) {
    c(NA, W$scales)
  } else {
    W$scales
  }

  total_energy <- sum(band_energy)
  pct <- if (total_energy == 0) rep(0, length(band_energy)) else
    100 * band_energy / total_energy

  data.frame(
    band = band_names,
    scale = scale_vals,
    energy = band_energy,
    pct = pct,
    stringsAsFactors = FALSE
  )
}
