# Filter kernels and Chebyshev-filtered application

# internal cache for Chebyshev coefficients
.cheby_cache <- new.env(parent = emptyenv())

cheby_cache_key <- function(kernel, K, lmax, M, jackson) {
  # Sample kernel at a few points to create unique fingerprint

  # This captures closure environment variables that affect kernel behavior
  test_lambdas <- c(0.1, 0.5, 1.0) * lmax
  kernel_fingerprint <- paste(signif(kernel(test_lambdas), 8), collapse = ",")
  paste(
    K,
    signif(lmax, 8),
    M,
    jackson,
    kernel_fingerprint,
    sep = "|"
  )
}

trim_cheby_coeffs <- function(coeffs, rel_tol = 1e-12, abs_tol = 1e-14) {
  if (!length(coeffs)) {
    return(coeffs)
  }
  scale <- max(abs(coeffs))
  keep <- which(abs(coeffs) > max(abs_tol, scale * rel_tol))
  if (!length(keep)) {
    return(coeffs[1])
  }
  coeffs[seq_len(max(keep))]
}

#' Chebyshev coefficients for a spectral kernel
#' @param kernel function of lambda
#' @param K polynomial order (>=1)
#' @param lmax spectral radius
#' @param M oversampling points (default 2*K)
#' @keywords internal
cheby_coeffs <- function(kernel, K, lmax, M = NULL, jackson = FALSE, cache = TRUE) {
  stopifnot(K >= 1, lmax > 0)
  M <- M %||% (2 * K)
  if (cache) {
    key <- cheby_cache_key(kernel, K, lmax, M, jackson)
    if (exists(key, envir = .cheby_cache, inherits = FALSE)) {
      return(get(key, envir = .cheby_cache, inherits = FALSE))
    }
  }
  n <- 0:(M - 1)
  theta <- pi * (n + 0.5) / M
  x <- cos(theta)              # in [-1,1]
  lambdas <- (x + 1) * lmax / 2
  f <- kernel(lambdas)
  coeffs <- numeric(K)
  for (k in 0:(K - 1)) {
    c_k <- sum(f * cos(k * theta)) * 2 / M
    if (k == 0) c_k <- c_k / 2
    coeffs[k + 1] <- c_k
  }
  if (jackson) {
    coeffs <- coeffs * jackson_weights(K)
  }
  coeffs <- trim_cheby_coeffs(coeffs)
  if (cache) {
    assign(key, coeffs, envir = .cheby_cache)
  }
  coeffs
}

# Jackson damping weights for Chebyshev series (first kind)
# Non-negative, monotonically decreasing, enforce uniform convergence.
# @keywords internal
jackson_weights <- function(K) {
  k <- 0:(K - 1)
  N <- K + 1
  w <- ((N - k) * cos(pi * k / N) + sin(pi * k / N) / tan(pi / N)) / N
  as.numeric(w)
}

# Built-in kernels -----------------------------------------------------------

#' Heat kernel exp(-t * lambda)
#' @param t diffusion time
#' @return kernel function
#' @export
kernel_heat <- function(t) {
  force(t)
  function(lambda) exp(-t * lambda)
}

#' Exponential kernel exp(- (lambda / alpha)^beta)
#' @param alpha scale parameter
#' @param beta shape parameter
#' @return kernel function
#' @export
kernel_exponential <- function(alpha = 1, beta = 1) {
  force(alpha); force(beta)
  function(lambda) exp(- (pmax(lambda, 0) / alpha) ^ beta)
}

#' Mexican hat-like kernel lambda * exp(-lambda / sigma)
#' @param sigma scale parameter
#' @return kernel function
#' @export
kernel_mexican_hat <- function(sigma = 1) {
  force(sigma)
  function(lambda) lambda * exp(-lambda / sigma)
}

#' Gabor kernel exp(-(lambda-mu)^2 / (2 sigma^2))
#' @param mu center frequency
#' @param sigma bandwidth
#' @return kernel function
#' @export
kernel_gabor <- function(mu, sigma) {
  force(mu); force(sigma)
  function(lambda) exp(- (lambda - mu)^2 / (2 * sigma^2))
}

#' Rectangular pass-band
#' @param low lower cutoff frequency
#' @param high upper cutoff frequency
#' @return kernel function
#' @export
kernel_rectangle <- function(low, high) {
  force(low); force(high)
  function(lambda) as.numeric(lambda >= low & lambda <= high)
}

#' Half-cosine transition band kernel
#' Smoothly transitions from 1 to 0 between low and high.
#' @param low lower cutoff frequency
#' @param high upper cutoff frequency
#' @return kernel function
#' @export
kernel_half_cosine <- function(low, high) {
  force(low); force(high)
  function(lambda) {
    y <- numeric(length(lambda))
    y[lambda <= low] <- 1
    inband <- lambda > low & lambda < high
    y[inband] <- 0.5 * (1 + cos(pi * (lambda[inband] - low) / (high - low)))
    y
  }
}

#' Tight frame-style band kernel (PyGSP-style)
#' Uses half-cosine for transition; lowpass to band to highpass coverage.
#' @param low lower cutoff
#' @param high upper cutoff
#' @export
kernel_tight_band <- function(low, high) {
  force(low); force(high)
  function(lambda) {
    y <- numeric(length(lambda))
    in1 <- lambda <= low
    in2 <- lambda >= high
    mid <- lambda > low & lambda < high
    y[in1] <- 1
    y[in2] <- 0
    y[mid] <- cos(pi/2 * (lambda[mid] - low)/(high - low))
    y
  }
}

#' Wave (cosine) kernel for wave equation propagation
#' h(lambda) = cos(t * sqrt(lambda))
#' @param t time parameter
#' @export
kernel_wave <- function(t) {
  force(t)
  function(lambda) cos(t * sqrt(pmax(lambda, 0)))
}

# Filter banks / wavelets ----------------------------------------------------
#' Smooth Meyer-like kernel (transition band)
#' @param low low cutoff
#' @param high high cutoff
#' @export
kernel_meyer <- function(low, high) {
  force(low); force(high)
  mid <- (low + high) / 2
  width <- (high - low) / 2
  sigma <- function(x) {
    u <- (x - low) / (high - low)
    u <- pmin(1, pmax(0, u))
    u^4 * (35 - 84 * u + 70 * u^2 - 20 * u^3)
  }
  function(lambda) {
    w <- abs(lambda - mid)
    mask <- lambda <= low | lambda >= high
    res <- numeric(length(lambda))
    res[mask] <- 0
    res[!mask] <- cos(pi / 2 * sigma(lambda[!mask]))
    res
  }
}

# Filter banks / wavelets ----------------------------------------------------

#' Meyer wavelet bank (lowpass + bandpasses)
#' @param scales multiplicative scales (e.g., 2^(0:J-1)), or a gsp_graph
#' @param lmax spectral radius
#' @param n_scales number of scales (if scales is NULL or a graph)
#' @return list of kernel functions
#' @export
meyer_wavelet_bank <- function(scales = NULL, lmax = NULL, n_scales = NULL) {
  # Accept a graph for convenience
  if (inherits(scales, "gsp_graph")) {
    g <- scales
    scales <- NULL
    if (is.null(lmax)) lmax <- lambda_max(g, normalized = g$normalized)
  }
  if (is.null(scales)) {
    stopifnot(!is.null(lmax))
    if (is.null(n_scales)) n_scales <- 6
    # PyGSP-like default: dyadic down to roughly 4/(3*lmax)
    scales <- (4 / (3 * lmax)) * 2 ^ seq(n_scales - 1, 0, by = -1)
  }
  scales <- sort(scales)
  nb <- length(scales)
  kernels <- vector("list", nb + 1)
  # lowpass
  low_cut <- lmax / max(scales) / 2
  kernels[[1]] <- kernel_rectangle(0, low_cut)
  for (i in seq_len(nb)) {
    s <- scales[i]
    low <- lmax / (2 * s)
    high <- lmax / s
    kernels[[i + 1]] <- kernel_meyer(low, high)
  }
  names(kernels) <- c("lowpass", paste0("band_", seq_len(nb)))
  kernels
}

#' Mexican-hat wavelet bank
#' @param scales scales vector, or a gsp_graph
#' @param n_scales number of scales (if scales is a graph)
#' @param lmax spectral radius
#' @param lpfactor factor for low-frequency bound
#' @return list of kernel functions
#' @export
mexican_hat_wavelet_bank <- function(scales, n_scales = NULL, lmax = NULL, lpfactor = 20) {
  # Allow calling with a graph for PyGSP parity: mexican_hat_wavelet_bank(g, n_scales)
  if (inherits(scales, "gsp_graph")) {
    g <- scales
    scales <- NULL
    if (is.null(lmax)) {
      lmax <- lambda_max(g, normalized = g$normalized)
    }
  }
  if (is.null(scales)) {
    stopifnot(!is.null(lmax), !is.null(n_scales))
    lmin <- lmax / lpfactor
    scales <- exp(seq(log(lmin), log(lmax), length.out = n_scales))
  }
  lapply(scales, function(s) kernel_mexican_hat(s))
}

#' Logarithmically spaced scales (PyGSP-compatible helper)
#'
#' @param lmin smallest non-zero eigenvalue proxy
#' @param lmax spectral radius
#' @param n_scales number of scales
#' @param t1 lower kernel support parameter
#' @param t2 upper kernel support parameter
#' @keywords internal
compute_log_scales <- function(lmin, lmax, n_scales, t1 = 1, t2 = 2) {
  stopifnot(lmin > 0, lmax > 0, n_scales >= 1)
  scale_min <- t1 / lmax
  scale_max <- t2 / lmin
  exp(seq(log(scale_max), log(scale_min), length.out = n_scales))
}

#' Abspline (A–B cubic spline) wavelet filter bank
#'
#' Port of `pygsp.filters.Abspline` as a list of kernel functions.
#'
#' @param lmax spectral radius, or a `gsp_graph`
#' @param Nf number of filters (default 6)
#' @param lpfactor low-pass factor (default 20)
#' @param scales optional vector of scales (length `Nf-1`)
#' @return list of `Nf` kernel functions (lowpass + wavelets)
#' @export
abspline_filter_bank <- function(lmax, Nf = 6, lpfactor = 20, scales = NULL) {
  if (inherits(lmax, "gsp_graph")) {
    g <- lmax
    lmax <- lambda_max(g, normalized = g$normalized)
  }
  lmax <- as.numeric(lmax)[1]
  Nf <- as.integer(Nf)[1]
  stopifnot(is.finite(lmax), lmax > 0, Nf >= 2, lpfactor > 0)

  # A–B cubic spline kernel (PyGSP defaults: alpha=beta=2, t1=1, t2=2)
  kernel_abspline3 <- local({
    alpha <- 2
    beta <- 2
    t1 <- 1
    t2 <- 2
    M <- rbind(
      c(1, t1, t1^2, t1^3),
      c(1, t2, t2^2, t2^3),
      c(0, 1, 2 * t1, 3 * t1^2),
      c(0, 1, 2 * t2, 3 * t2^2)
    )
    v <- c(1, 1, alpha, -beta / 2)
    a <- solve(M, v)

    function(x) {
      x <- pmax(x, 0)
      r <- numeric(length(x))
      r1 <- x <= t1
      r2 <- x > t1 & x < t2
      r3 <- x >= t2
      if (any(r1)) r[r1] <- x[r1]^alpha * t1^(-alpha)
      if (any(r2)) {
        x2 <- x[r2]
        r[r2] <- a[1] + a[2] * x2 + a[3] * x2^2 + a[4] * x2^3
      }
      if (any(r3)) r[r3] <- x[r3]^(-beta) * t2^beta
      r
    }
  })

  gl <- function(x) exp(-x^4)
  lmin <- lmax / lpfactor

  if (is.null(scales)) {
    scales <- compute_log_scales(lmin, lmax, Nf - 1)
  }
  if (length(scales) != (Nf - 1)) stop("len(scales) must be Nf-1")
  scales <- as.numeric(scales)

  gb <- function(x) kernel_abspline3(x)

  # gamma_l = max_{x in [1,2]} gb(x)
  grid <- seq(1, 2, length.out = 2048)
  gamma_l <- max(gb(grid))
  lminfac <- 0.6 * lmin

  kernels <- vector("list", Nf)
  kernels[[1]] <- local({
    gl0 <- gl
    gam <- gamma_l
    lmf <- lminfac
    function(lambda) gam * gl0(lambda / lmf)
  })

  for (k in seq_len(Nf - 1)) {
    kernels[[k + 1]] <- local({
      s <- scales[k]
      function(lambda) gb(s * lambda)
    })
  }

  names(kernels) <- c("lowpass", paste0("wavelet_", seq_len(Nf - 1)))
  attr(kernels, "scales") <- scales
  attr(kernels, "lpfactor") <- lpfactor
  kernels
}

#' Itersine tight-frame filter bank
#'
#' Port of `pygsp.filters.Itersine` as a list of kernel functions.
#'
#' @param lmax spectral radius, or a `gsp_graph`
#' @param Nf number of filters from 0..lmax (default 6)
#' @param overlap overlap parameter (default 2)
#' @return list of `Nf` kernel functions
#' @export
itersine_filter_bank <- function(lmax, Nf = 6, overlap = 2) {
  if (inherits(lmax, "gsp_graph")) {
    g <- lmax
    lmax <- lambda_max(g, normalized = g$normalized)
  }
  lmax <- as.numeric(lmax)[1]
  Nf <- as.integer(Nf)[1]
  overlap <- as.numeric(overlap)[1]
  stopifnot(is.finite(lmax), lmax > 0, Nf >= 2, overlap > 0)

  scale <- lmax / (Nf - overlap + 1) * overlap

  base_kernel <- function(x) {
    y <- cos(x * pi)^2
    y <- sin(0.5 * pi * y)
    y * as.numeric(x >= -0.5 & x <= 0.5)
  }

  kernels <- vector("list", Nf)
  for (k in seq_len(Nf)) {
    kernels[[k]] <- local({
      kk <- k
      sc <- scale
      ov <- overlap
      function(lambda) {
        x <- lambda / sc - (kk - ov / 2) / ov
        base_kernel(x) * sqrt(2 / ov)
      }
    })
  }
  names(kernels) <- paste0("band_", seq_len(Nf))
  attr(kernels, "overlap") <- overlap
  kernels
}

#' Simple tight-frame wavelet filter bank
#'
#' Port of `pygsp.filters.SimpleTight` as a list of kernel functions.
#'
#' @param lmax spectral radius, or a `gsp_graph`
#' @param Nf number of filters covering `[0, lmax]` (default 6)
#' @param scales optional scales (length `Nf-1`)
#' @return list of `Nf` kernel functions (scaling + wavelets)
#' @export
simpletight_filter_bank <- function(lmax, Nf = 6, scales = NULL) {
  if (inherits(lmax, "gsp_graph")) {
    g <- lmax
    lmax <- lambda_max(g, normalized = g$normalized)
  }
  lmax <- as.numeric(lmax)[1]
  Nf <- as.integer(Nf)[1]
  stopifnot(is.finite(lmax), lmax > 0, Nf >= 2)

  if (is.null(scales)) {
    scales <- (1 / (2 * lmax)) * 2^(seq(Nf - 2, 0, by = -1))
  }
  if (length(scales) != (Nf - 1)) stop("len(scales) must be Nf-1")
  scales <- as.numeric(scales)

  kernel_simpletight <- function(x, type = c("sf", "wavelet")) {
    type <- match.arg(type)
    l1 <- 0.25
    l2 <- 0.5
    l3 <- 1.0
    h <- function(z) sin(pi * z / 2)^2
    r <- numeric(length(x))
    r1 <- x < l1
    r2 <- x >= l1 & x < l2
    r3 <- x >= l2 & x < l3

    if (type == "sf") {
      r[r1] <- 1
      if (any(r2)) r[r2] <- sqrt(pmax(0, 1 - h(4 * x[r2] - 1)^2))
    } else {
      if (any(r2)) r[r2] <- h(4 * (x[r2] - 1 / 4))
      if (any(r3)) r[r3] <- sqrt(pmax(0, 1 - h(2 * x[r3] - 1)^2))
    }
    r
  }

  kernels <- vector("list", Nf)
  kernels[[1]] <- local({
    s0 <- scales[1]
    function(lambda) kernel_simpletight(s0 * lambda, "sf")
  })
  for (k in seq_len(Nf - 1)) {
    kernels[[k + 1]] <- local({
      s <- scales[k]
      function(lambda) kernel_simpletight(s * lambda, "wavelet")
    })
  }

  names(kernels) <- c("scaling", paste0("wavelet_", seq_len(Nf - 1)))
  attr(kernels, "scales") <- scales
  kernels
}

#' Tight frame filter bank (lowpass + bandpasses)
#'
#' Creates a partition-of-unity filter bank using half-cosine / half-sine
#' transitions at each edge.  For K+1 edges the bank contains K filters
#' whose squared responses sum to 1 everywhere.
#'
#' @param edges monotone increasing vector of cutoff edges (len >=2)
#' @export
tight_frame_bank <- function(edges) {
  stopifnot(length(edges) >= 2)
  edges <- sort(edges)
  K <- length(edges) - 1
  kernels <- vector("list", K)
  for (k in seq_len(K)) {
    kernels[[k]] <- local({
      kk <- k; KK <- K; ee <- edges
      function(lambda) {
        y <- numeric(length(lambda))
        if (kk == 1L) {
          # Lowpass: flat 1 below e[2], cos fall in [e[2], e[3]]
          y[lambda <= ee[2]] <- 1
          if (KK >= 2L) {
            tr <- lambda > ee[2] & lambda < ee[3]
            y[tr] <- cos(pi / 2 * (lambda[tr] - ee[2]) / (ee[3] - ee[2]))
          }
        } else if (kk == KK) {
          # Highpass: sin rise in [e[K], e[K+1]], flat 1 above e[K+1]
          tr <- lambda > ee[KK] & lambda < ee[KK + 1L]
          y[tr] <- sin(pi / 2 * (lambda[tr] - ee[KK]) / (ee[KK + 1L] - ee[KK]))
          y[lambda >= ee[KK + 1L]] <- 1
        } else {
          # Middle band: sin rise in [e[k], e[k+1]], cos fall in [e[k+1], e[k+2]]
          lo <- ee[kk]; mid_e <- ee[kk + 1L]; hi <- ee[kk + 2L]
          rise <- lambda > lo & lambda < mid_e
          y[rise] <- sin(pi / 2 * (lambda[rise] - lo) / (mid_e - lo))
          y[lambda == mid_e] <- 1
          fall <- lambda > mid_e & lambda < hi
          y[fall] <- cos(pi / 2 * (lambda[fall] - mid_e) / (hi - mid_e))
        }
        y
      }
    })
  }
  names(kernels) <- paste0("band_", seq_len(K))
  kernels
}

#' Tight frame bank presets (Regular/Held/Simoncelli/Papadakis style)
#' @param lmax spectral radius
#' @param J number of dyadic bands
#' @export
tight_frame_regular <- function(lmax, J = 3) {
  if (inherits(lmax, "gsp_graph")) {
    lmax <- lambda_max(lmax, normalized = lmax$normalized)
  }
  edges <- lmax / 2 ^ seq(J, 0, by = -1)
  tight_frame_bank(c(0, edges))
}

#' @rdname tight_frame_regular
#' @export
tight_frame_held <- function(lmax, J = 3) {
  tight_frame_regular(lmax, J)
}

#' @rdname tight_frame_regular
#' @export
tight_frame_simoncelli <- function(lmax, J = 3) {
  tight_frame_regular(lmax, J)
}

#' @rdname tight_frame_regular
#' @export
tight_frame_papadakis <- function(lmax, J = 3) {
  tight_frame_regular(lmax, J)
}

# Filtering helpers ----------------------------------------------------------

#' Apply a spectral kernel via Chebyshev approximation
#' @param g graph
#' @param signal vector or matrix (n x m)
#' @param kernel function of lambda
#' @param K Chebyshev order
#' @param lmax optional precomputed lambda_max
#' @param jackson logical; apply Jackson damping
#' @param threads optional number of OpenMP threads
#' @param precision numeric precision (currently only "double")
#' @param strategy filtering strategy: "auto", "cheby", "lanczos", or "exact"
#' @return filtered signal (vector or matrix)
#' @export
filter_signal <- function(g, signal, kernel, K = 30, lmax = NULL, jackson = FALSE,
                          threads = NULL,
                          precision = c("double"),
                          strategy = c("auto", "cheby", "lanczos", "exact")) {
  validate_signal_dims(g, signal)
  if (is.null(dim(signal))) signal <- matrix(signal, ncol = 1)
  precision <- match.arg(precision)
  strategy <- match.arg(strategy)

  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)

  # threads control
  threads <- threads %||% getOption("rgsp.threads", NULL)
  if (!is.null(threads)) {
    set_omp_threads_cpp(as.integer(threads))
  }

  # Choose strategy
  if (strategy == "auto") {
    if (filter_signal_exact_auto_ok(g)) {
      strategy <- "exact"
    } else if (!jackson && K >= 60) {
      strategy <- "lanczos"
    } else {
      strategy <- "cheby"
    }
  }

  if (strategy == "lanczos") {
    return(filter_signal_lanczos(g, signal, kernel, K = K))
  }

  if (strategy == "exact") {
    return(filter_signal_exact(g, signal, kernel, lmax = lmax))
  }

  coeffs <- cheby_coeffs(kernel, K = K, lmax = lmax, jackson = jackson)
  if (!isTRUE(g$normalized) && !isTRUE(g$directed) &&
      !is.null(g$graph_type) && g$graph_type %in% c("ring", "path")) {
    out <- chebyshev_filter_ring_cpp(
      signal,
      coeffs,
      lmax,
      periodic = identical(g$graph_type, "ring")
    )
    return(if (ncol(out) == 1) drop(out) else out)
  }
  L_sp <- g$cache$laplacian_sp
  if (is.null(L_sp)) {
    L <- graph_laplacian(g, g$normalized)
    L_sp <- laplacian_sp_xptr_cpp(L)
    g$cache$laplacian_sp <- L_sp
  }
  prec_flag <- if (precision == "float") 1L else 0L
  out <- chebyshev_filter_omp_cpp(L_sp, signal, coeffs, lmax, threads = threads %||% -1L, precision = prec_flag)
  if (ncol(out) == 1) drop(out) else out
}

filter_signal_exact_auto_ok <- function(g) {
  if (isTRUE(g$normalized) || isTRUE(g$directed)) {
    return(FALSE)
  }

  if (identical(g$graph_type, "ring")) {
    return(TRUE)
  }

  if (identical(g$graph_type, "path")) {
    return(TRUE)
  }

  if (identical(g$graph_type, "grid2d")) {
    return(!isTRUE(g$graph_params$periodic))
  }

  FALSE
}

path_dct_forward <- function(signal) {
  if (is.null(dim(signal))) signal <- matrix(signal, ncol = 1)
  n <- nrow(signal)
  ext <- rbind(signal, signal[n:1, , drop = FALSE])
  coeffs <- stats::mvfft(ext)
  k <- 0:(n - 1)
  phase <- exp((-1i * pi * k) / (2 * n))
  out <- (sqrt(2 / n) / 2) * sweep(coeffs[1:n, , drop = FALSE], 1, phase, `*`)
  out[1, ] <- out[1, ] / sqrt(2)
  Re(out)
}

path_dct_inverse <- function(coeffs) {
  if (is.null(dim(coeffs))) coeffs <- matrix(coeffs, ncol = 1)
  n <- nrow(coeffs)
  k <- 0:(n - 1)
  phase <- exp((1i * pi * k) / (2 * n))
  scaled <- coeffs
  scaled[1, ] <- scaled[1, ] * sqrt(2)
  spectrum <- (2 / sqrt(2 / n)) * sweep(scaled, 1, phase, `*`)
  ext <- matrix(0 + 0i, nrow = 2 * n, ncol = ncol(coeffs))
  ext[1:n, ] <- spectrum
  if (n > 1) {
    ext[(n + 2):(2 * n), ] <- Conj(spectrum[n:2, , drop = FALSE])
  }
  signal <- stats::mvfft(ext, inverse = TRUE) / (2 * n)
  Re(signal[1:n, , drop = FALSE])
}

path_spectral_data <- function(n) {
  j <- 0:(n - 1)
  k <- 0:(n - 1)
  U <- sqrt(2 / n) * cos(outer(j + 0.5, k) * pi / n)
  U[, 1] <- U[, 1] / sqrt(2)
  lambda <- 2 - 2 * cos(pi * k / n)
  list(values = lambda, vectors = U)
}

grid2d_spectral_data <- function(nrow, ncol) {
  row_spec <- path_spectral_data(nrow)
  col_spec <- path_spectral_data(ncol)
  list(
    row_vectors = row_spec$vectors,
    row_values = row_spec$values,
    col_vectors = col_spec$vectors,
    col_values = col_spec$values
  )
}

#' Apply a spectral kernel exactly
#'
#' Uses graph-specific fast transforms where available and falls back to the
#' eigendecomposition otherwise.
#'
#' @param g graph
#' @param signal vector or matrix (n x m)
#' @param kernel function of eigenvalues
#' @param lmax optional spectral radius
#' @return filtered signal (vector or matrix)
#' @export
filter_signal_exact <- function(g, signal, kernel, lmax = NULL) {
  if (is.null(dim(signal))) signal <- matrix(signal, ncol = 1)

  if (!isTRUE(g$normalized) && !isTRUE(g$directed) &&
      identical(g$graph_type, "ring")) {
    n <- nrow(signal)
    idx <- 0:(n - 1)
    lambda <- 2 - 2 * cos(2 * pi * idx / n)
    H <- kernel(lambda)
    Xhat <- stats::mvfft(signal)
    Yhat <- sweep(Xhat, 1, H, `*`)
    out <- stats::mvfft(Yhat, inverse = TRUE) / n
    out <- Re(out)
    return(if (ncol(out) == 1) drop(out) else out)
  }

  if (!isTRUE(g$normalized) && !isTRUE(g$directed) &&
      identical(g$graph_type, "path")) {
    n <- nrow(signal)
    lambda <- 2 - 2 * cos(pi * (0:(n - 1)) / n)
    coeffs <- path_dct_forward(signal)
    out <- path_dct_inverse(kernel(lambda) * coeffs)
    return(if (ncol(out) == 1) drop(out) else out)
  }

  if (!isTRUE(g$normalized) && !isTRUE(g$directed) &&
      identical(g$graph_type, "grid2d") &&
      !isTRUE(g$graph_params$periodic)) {
    spec <- g$cache$grid2d_spec
    if (is.null(spec)) {
      spec <- grid2d_spectral_data(g$graph_params$nrow, g$graph_params$ncol)
      g$cache$grid2d_spec <- spec
    }
    Ur <- spec$row_vectors
    Uc <- spec$col_vectors
    Hr <- spec$row_values
    Hc <- spec$col_values
    H <- outer(Hr, Hc, `+`)

    out <- matrix(0, nrow = nrow(signal), ncol = ncol(signal))
    for (i in seq_len(ncol(signal))) {
      Xi <- matrix(signal[, i], nrow = g$graph_params$nrow, ncol = g$graph_params$ncol, byrow = TRUE)
      Ci <- crossprod(Ur, Xi) %*% Uc
      Yi <- Ur %*% (kernel(H) * Ci) %*% t(Uc)
      out[, i] <- as.vector(t(Yi))
    }
    return(if (ncol(out) == 1) drop(out) else out)
  }

  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  H <- kernel(eig$values)
  coeffs <- t(U) %*% signal
  out <- U %*% (H * coeffs)
  if (ncol(out) == 1) drop(out) else out
}

#' Apply a spectral kernel via Lanczos approximation
#'
#' Approximates f(L) x using a K-step Lanczos tridiagonal projection.
#' Typically more accurate than Chebyshev for ill-conditioned spectra.
#'
#' @param g graph
#' @param signal vector or matrix (n x m)
#' @param kernel function of eigenvalues
#' @param K Lanczos steps (<= n)
#' @export
filter_signal_lanczos <- function(g, signal, kernel, K = 30) {
  if (is.null(dim(signal))) signal <- matrix(signal, ncol = 1)
  L <- graph_laplacian(g, normalized = g$normalized)
  n <- nrow(L)
  K <- min(as.integer(K), n)
  res <- lanczos_filter_cpp(L, signal, K, kernel)
  if (ncol(res) == 1) drop(res) else res
}

#' Apply a bank of kernels
#' @param g graph
#' @param signal vector or matrix
#' @param kernels named list of kernel functions
#' @param K Chebyshev order
#' @param lmax optional spectral radius
#' @param jackson logical; apply Jackson damping
#' @return list of filtered signals (same order as kernels)
#' @export
filter_bank_apply <- function(g, signal, kernels, K = 30, lmax = NULL, jackson = FALSE) {
  validate_signal_dims(g, signal)
  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)
  lapply(kernels, function(kern) filter_signal(g, signal, kern, K = K, lmax = lmax, jackson = jackson))
}

# Wavelets -------------------------------------------------------------------

#' Heat wavelet bank at given scales
#' @param scales numeric vector of scales
#' @export
wavelet_heat_bank <- function(scales) {
  lapply(scales, function(s) kernel_heat(s))
}

#' Compute heat wavelet transform
#' @param g graph
#' @param signal vector/matrix
#' @param scales scales for heat kernels
#' @param K Chebyshev order for filtering
#' @param lmax optional spectral radius
#' @param jackson logical; apply Jackson damping
#' @return list of filtered signals, one per scale
#' @export
wavelet_heat_transform <- function(g, signal, scales, K = 30, lmax = NULL, jackson = FALSE) {
  validate_signal_dims(g, signal)
  bank <- wavelet_heat_bank(scales)
  res <- filter_bank_apply(g, signal, bank, K = K, lmax = lmax, jackson = jackson)
  names(res) <- paste0("scale_", scales)
  res
}

#' Heat wavelet transform on a k-NN graph (convenience)
#'
#' Build a k-nearest-neighbor graph from coordinates and apply the heat wavelet
#' transform in one step.
#'
#' @param coords numeric matrix of node coordinates (n x d)
#' @param signal vector or matrix of signals (n x m)
#' @param scales numeric vector of heat scales
#' @param k neighbors per node (passed to [graph_knn()])
#' @param weight weighting scheme for k-NN edges (`"heat"`, `"binary"`, `"distance"`)
#' @param sigma optional bandwidth for `weight = "heat"`; defaults to median k-NN distance
#' @param sym symmetrization mode for k-NN graph (`"union"` or `"mutual"`)
#' @param normalized logical; compute normalized Laplacian
#' @param seed optional RNG seed
#' @param K Chebyshev order
#' @param lmax optional spectral radius
#' @param jackson logical; apply Jackson damping
#'
#' @return list of filtered signals, one per scale (as in [wavelet_heat_transform()])
#' @export
wavelet_heat_knn <- function(coords, signal, scales,
                             k = 6,
                             weight = c("heat", "binary", "distance"),
                             sigma = NULL,
                             sym = c("union", "mutual"),
                             normalized = FALSE,
                             seed = NULL,
                             K = 30,
                             lmax = NULL,
                             jackson = FALSE) {
  g <- graph_knn(coords,
                 k = k,
                 weight = weight,
                 sigma = sigma,
                 sym = sym,
                 normalized = normalized,
                 seed = seed)
  wavelet_heat_transform(g, signal, scales = scales, K = K, lmax = lmax, jackson = jackson)
}

# Gabor / Modulation ---------------------------------------------------------

#' Gabor filter bank
#' @param mus center frequencies
#' @param sigma bandwidth
#' @export
gabor_filter_bank <- function(mus, sigma) {
  lapply(mus, function(mu) kernel_gabor(mu, sigma))
}

#' Modulation / Gabor vertex-frequency transform
#' @param g graph
#' @param signal vector or matrix (n x m)
#' @param mus vector of center frequencies
#' @param sigma bandwidth of gabor kernel
#' @param K Chebyshev order
#' @param lmax optional spectral radius
#' @param jackson logical; use Jackson damping
#' @return list of filtered signals (one per mu)
#' @export
modulation_transform <- function(g, signal, mus, sigma, K = 30, lmax = NULL, jackson = FALSE) {
  bank <- gabor_filter_bank(mus, sigma)
  filter_bank_apply(g, signal, bank, K = K, lmax = lmax, jackson = jackson)
}

# Diffusion / heat propagation ----------------------------------------------

#' Heat diffusion of a signal for one or multiple times
#' @param g graph
#' @param signal vector/matrix
#' @param t scalar or vector of times
#' @param K Chebyshev order
#' @param lmax optional spectral radius
#' @return filtered signal(s), or list if multiple times
#' @export
heat_propagate <- function(g, signal, t, K = 30, lmax = NULL) {
  times <- as.numeric(t)
  lmax <- lmax %||% lambda_max(g, normalized = g$normalized)
  out_list <- lapply(times, function(tt) {
    filter_signal(g, signal, kernel_heat(tt), K = K, lmax = lmax)
  })
  if (length(times) == 1) {
    out_list[[1]]
  } else {
    names(out_list) <- paste0("t_", times)
    out_list
  }
}
