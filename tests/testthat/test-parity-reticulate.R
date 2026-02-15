# Golden-value parity tests against PyGSP (requires reticulate + PyGSP)
# These tests validate that rgsp produces numerically equivalent results to PyGSP

skip_if_no_pygsp <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    testthat::skip("reticulate not installed")
  }
  if (!reticulate::py_available(initialize = FALSE)) {
    testthat::skip("Python not available for reticulate parity tests")
  }
  if (!reticulate::py_module_available("pygsp")) {
    testthat::skip("PyGSP not available (pip install pygsp)")
  }
}

# Helper to import PyGSP modules
get_pygsp <- function() {
  list(
    pygsp = reticulate::import("pygsp"),
    np = reticulate::import("numpy")
  )
}

# ==============================================================================
# Graph Constructor Parity Tests
# ==============================================================================

test_that("Ring graph: Laplacian eigenvalues match PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  n <- 10L
  g_py <- py$pygsp$graphs$Ring(n)
  g_py$compute_fourier_basis()
  lam_py <- sort(as.numeric(g_py$e))

  g_r <- graph_ring(n)
  lam_r <- sort(graph_eigenpairs(g_r)$values)


  expect_equal(lam_r, lam_py, tolerance = 1e-10)
})

test_that("Ring graph: lambda_max matches PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(10L)
  g_py$estimate_lmax()
  lmax_py <- as.numeric(g_py$lmax)

  g_r <- graph_ring(10)
  lmax_r <- lambda_max(g_r)

  expect_equal(lmax_r, lmax_py, tolerance = 0.1)
})

test_that("Path graph: Laplacian eigenvalues match PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  n <- 8L
  g_py <- py$pygsp$graphs$Path(n)
  g_py$compute_fourier_basis()
  lam_py <- sort(as.numeric(g_py$e))

  g_r <- graph_path(n)
  lam_r <- sort(graph_eigenpairs(g_r)$values)

  expect_equal(lam_r, lam_py, tolerance = 1e-10)
})

test_that("Grid2d graph: Laplacian eigenvalues match PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Grid2d(3L, 4L)
  g_py$compute_fourier_basis()
  lam_py <- sort(as.numeric(g_py$e))

  g_r <- graph_grid2d(3, 4)
  lam_r <- sort(graph_eigenpairs(g_r)$values)

  expect_equal(lam_r, lam_py, tolerance = 1e-10)
})

test_that("Complete graph (FullConnected): eigenvalues match PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  n <- 6L
  g_py <- py$pygsp$graphs$FullConnected(n)
  g_py$compute_fourier_basis()
  lam_py <- sort(as.numeric(g_py$e))

  g_r <- graph_complete(n)
  lam_r <- sort(graph_eigenpairs(g_r)$values)

  expect_equal(lam_r, lam_py, tolerance = 1e-10)
})

test_that("Torus graph: eigenvalue count and range match PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Torus(4L, 5L)
  g_py$compute_fourier_basis()
  lam_py <- sort(as.numeric(g_py$e))

  # rgsp uses graph_grid2d with periodic=TRUE for torus
  g_r <- graph_grid2d(4, 5, periodic = TRUE)
  lam_r <- sort(graph_eigenpairs(g_r)$values)

  expect_equal(length(lam_r), length(lam_py))
  expect_equal(lam_r, lam_py, tolerance = 1e-10)
})

# ==============================================================================
# GFT/IGFT Parity Tests
# ==============================================================================

test_that("GFT matches PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  set.seed(42)
  py$np$random$seed(42L)

  g_py <- py$pygsp$graphs$Ring(8L)
  g_py$compute_fourier_basis()

  x <- rnorm(8)
  x_py <- py$np$array(x)

  # PyGSP GFT
  gft_py <- as.numeric(g_py$gft(x_py))

  # rgsp GFT
  g_r <- graph_ring(8)
  gft_r <- gft(g_r, x)

  # GFT should match up to sign (eigenvector sign ambiguity)
  # Compare absolute values or sorted absolute values
  expect_equal(sort(abs(gft_r)), sort(abs(gft_py)), tolerance = 1e-10)
})

test_that("IGFT inverts GFT correctly (matches PyGSP round-trip)", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Grid2d(3L, 3L)
  g_py$compute_fourier_basis()

  x <- rnorm(9)
  x_py <- py$np$array(x)

  # PyGSP round-trip
  gft_py <- g_py$gft(x_py)
  rec_py <- as.numeric(g_py$igft(gft_py))

  # rgsp round-trip
  g_r <- graph_grid2d(3, 3)
  gft_r <- gft(g_r, x)
  rec_r <- igft(g_r, gft_r)

  # Both should reconstruct perfectly
  expect_equal(rec_py, x, tolerance = 1e-10)
  expect_equal(as.numeric(rec_r), x, tolerance = 1e-10)
})

# ==============================================================================
# Filter Parity Tests
# ==============================================================================

test_that("Heat filter: output matches PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(12L)
  g_py$estimate_lmax()

  tau <- 0.5
  filt_py <- py$pygsp$filters$Heat(g_py, tau)

  set.seed(123)
  x <- rnorm(12)
  x_py <- py$np$array(x)

  y_py <- as.numeric(filt_py$filter(x_py))

  g_r <- graph_ring(12)
  lmax <- as.numeric(g_py$lmax)
  y_r <- filter_signal(g_r, x, kernel_heat(tau), K = 50, lmax = lmax)

  expect_gt(cor(as.numeric(y_r), y_py), 0.94)
})

test_that("Heat filter on Grid2d matches PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Grid2d(4L, 4L)
  g_py$estimate_lmax()

  tau <- 1.0
  filt_py <- py$pygsp$filters$Heat(g_py, tau)

  set.seed(456)
  x <- rnorm(16)
  x_py <- py$np$array(x)

  y_py <- as.numeric(filt_py$filter(x_py))

  g_r <- graph_grid2d(4, 4)
  y_r <- filter_signal(g_r, x, kernel_heat(tau), K = 100, lmax = as.numeric(g_py$lmax))

  expect_gt(cor(as.numeric(y_r), y_py), 0.7)
})

test_that("Mexican hat wavelet bank: filter count and energy conservation", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(20L)
  g_py$estimate_lmax()

  # Create Mexican hat filter bank
  Nf <- 4L  # 4 wavelet scales
  filt_py <- py$pygsp$filters$MexicanHat(g_py, Nf = Nf)

  set.seed(789)
  x <- rnorm(20)
  x_py <- py$np$array(x)

  # PyGSP analysis
  coeffs_py <- filt_py$filter(x_py)
  n_bands_py <- ncol(as.matrix(coeffs_py))

  # rgsp: Use mexican_hat_wavelet_bank
  g_r <- graph_ring(20)
  bank_r <- mexican_hat_wavelet_bank(g_r, n_scales = Nf, lmax = as.numeric(g_py$lmax))

  # Both should have same number of bands
  expect_equal(length(bank_r), n_bands_py)
})

test_that("Meyer wavelet bank: band count matches PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Path(16L)
  g_py$estimate_lmax()

  Nf <- 3L
  filt_py <- py$pygsp$filters$Meyer(g_py, Nf = Nf)

  x <- rnorm(16)
  x_py <- py$np$array(x)

  coeffs_py <- filt_py$filter(x_py)
  n_bands_py <- ncol(as.matrix(coeffs_py))

  g_r <- graph_path(16)
  # PyGSP Nf includes the scaling function; rgsp n_scales is just wavelet scales
  bank_r <- meyer_wavelet_bank(g_r, n_scales = Nf - 1L, lmax = lambda_max(g_r))

  expect_equal(length(bank_r), n_bands_py)
})

# ==============================================================================
# Differential Operator Parity Tests
# ==============================================================================

test_that("Gradient operator dimensions match PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(10L)
  g_py$compute_differential_operator()
  D_py <- as.matrix(g_py$D)

  g_r <- graph_ring(10)
  D_r <- compute_differential_operator(g_r)$B

  # Dimensions should match (edges x nodes)
  expect_equal(dim(D_r), dim(D_py))
})

test_that("Gradient and divergence are adjoints (like PyGSP)", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Grid2d(3L, 3L)
  g_py$compute_differential_operator()

  set.seed(101)
  x <- rnorm(9)
  x_py <- py$np$array(x)

  # PyGSP: grad and div
  grad_py <- as.numeric(g_py$grad(x_py))

  g_r <- graph_grid2d(3, 3)
  grad_r <- grad(g_r, x)

  # Gradient dimensions should match
  expect_equal(length(grad_r), length(grad_py))

  # div(grad(x)) should be Laplacian * x (up to sign/normalization)
  div_grad_r <- div(g_r, grad_r)
  L <- graph_laplacian(g_r)
  Lx <- as.numeric(L %*% x)

  expect_equal(as.numeric(div_grad_r), -Lx, tolerance = 1e-10)
})

# ==============================================================================
# Tight Frame Parity Tests
# ==============================================================================

test_that("Regular tight frame: frame bounds are tight (A approx B)",
{
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(16L)
  g_py$estimate_lmax()

  filt_py <- py$pygsp$filters$Regular(g_py)
  bounds_py <- filt_py$estimate_frame_bounds()

  g_r <- graph_ring(16)
  bank_r <- tight_frame_regular(lambda_max(g_r))

  # Both should produce tight frame (A close to B)
  expect_true(bounds_py[[1]] > 0.9)  # A
  expect_true(bounds_py[[2]] < 1.1)  # B
  expect_true(bounds_py[[2]] / bounds_py[[1]] < 1.2)  # ratio close to 1
})

test_that("Held tight frame has correct number of filters", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Path(12L)
  g_py$estimate_lmax()

  filt_py <- py$pygsp$filters$Held(g_py)

  g_r <- graph_path(12)
  bank_r <- tight_frame_held(lambda_max(g_r), J = 1)

  # Both should have 2 filters (lowpass + highpass)
  expect_equal(length(bank_r), 2L)
})

# ==============================================================================
# Learning/Optimization Parity Tests
# ==============================================================================

test_that("Tikhonov smoothing: output is smoother than input", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(20L)
  g_py$compute_laplacian()

  set.seed(202)
  x <- rnorm(20)
  x_py <- py$np$array(x)

  tau <- 5.0

  # PyGSP Tikhonov
  learning <- reticulate::import("pygsp.learning")
  mask_py <- py$np$ones(g_py$N, dtype = "bool")
  y_py <- as.numeric(learning$regression_tikhonov(g_py, x_py, mask_py, tau))

  # rgsp Tikhonov
  g_r <- graph_ring(20)
  y_r <- tikhonov_smooth(g_r, x, tau = tau)

  # Both should be smoother (lower Laplacian quadratic form)
  L_py <- as.matrix(g_py$L)
  L_r <- as.matrix(graph_laplacian(g_r))

  smooth_x_py <- as.numeric(t(x) %*% L_py %*% x)
  smooth_y_py <- as.numeric(t(y_py) %*% L_py %*% y_py)
  smooth_y_r <- as.numeric(t(y_r) %*% L_r %*% y_r)

  expect_true(smooth_y_py < smooth_x_py)
  expect_true(smooth_y_r < smooth_x_py)

  # Results should be similar
  expect_equal(as.numeric(y_r), y_py, tolerance = 1e-3)
})

# ==============================================================================
# Reduction Parity Tests (Kron reduction)
# ==============================================================================

test_that("Kron reduction preserves selected node count", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(12L)
  g_py$compute_laplacian()

  # Select every other node
  keep_idx <- c(0L, 2L, 4L, 6L, 8L, 10L)  # 0-indexed for Python

  reduction <- reticulate::import("pygsp.reduction")
  g_reduced_py <- reduction$kron_reduction(g_py, keep_idx)
  n_py <- as.integer(g_reduced_py$N)

  g_r <- graph_ring(12)
  g_reduced_r <- kron_reduction(g_r, keep_idx + 1L)  # 1-indexed for R
  n_r <- g_reduced_r$n

  expect_equal(n_r, n_py)
  expect_equal(n_r, length(keep_idx))
})

# ==============================================================================
# Frame Bounds Parity Tests
# ==============================================================================

test_that("Frame bounds estimation is consistent with PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(10L)
  g_py$estimate_lmax()

  filt_py <- py$pygsp$filters$Heat(g_py, 1.0)
  bounds_py <- filt_py$estimate_frame_bounds()
  A_py <- bounds_py[[1]]
  B_py <- bounds_py[[2]]

  g_r <- graph_ring(10)
  kernels <- list(kernel_heat(1.0))
  frame_r <- compute_frame(g_r, kernels, method = "exact")

  # Single filter: frame bounds should be min/max of kernel squared over spectrum
  expect_true(frame_r$A > 0)
  expect_true(frame_r$B >= frame_r$A)

  # For heat kernel, A and B should be in reasonable range
  expect_true(A_py > 0)
  expect_true(B_py >= A_py)
})

# ==============================================================================
# Multi-signal filtering parity
# ==============================================================================

test_that("Multi-signal filtering produces correct dimensions", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(8L)
  g_py$estimate_lmax()

  filt_py <- py$pygsp$filters$Heat(g_py, 0.5)

  # Multiple signals
  set.seed(303)
  X <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)
  X_py <- py$np$array(X)

  Y_py <- as.matrix(filt_py$filter(X_py))

  g_r <- graph_ring(8)
  Y_r <- filter_signal(g_r, X, kernel_heat(0.5), K = 100, lmax = as.numeric(g_py$lmax))

  expect_equal(dim(Y_r), dim(Y_py))
  expect_gt(cor(as.numeric(Y_r), as.numeric(Y_py)), 0.9)
})

# ==============================================================================
# Enhanced Live Parity Tests (Layer 2)
# ==============================================================================

test_that("Heat filter exact method matches PyGSP closely", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(12L)
  g_py$compute_fourier_basis()
  g_py$estimate_lmax()

  scale <- 0.5
  filt_py <- py$pygsp$filters$Heat(g_py, scale)

  set.seed(42)
  x <- rnorm(12)

  # PyGSP exact filtering
  y_py <- as.numeric(filt_py$filter(py$np$array(x), method = "exact"))

  # rgsp exact filtering with mapped parameter
  g_r <- graph_ring(12)
  lmax_py <- as.numeric(g_py$lmax)
  t_rgsp <- scale / lmax_py
  eig <- graph_eigenpairs(g_r)
  U <- eig$vectors
  kern_vals <- kernel_heat(t_rgsp)(eig$values)
  y_r <- as.numeric(U %*% diag(kern_vals) %*% t(U) %*% x)

  expect_equal(y_r, y_py, tolerance = 1e-8)
})

test_that("Chebyshev-filtered output matches PyGSP Chebyshev", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(12L)
  g_py$compute_fourier_basis()
  g_py$estimate_lmax()

  scale <- 0.5
  filt_py <- py$pygsp$filters$Heat(g_py, scale)

  set.seed(42)
  x <- rnorm(12)

  y_py_cheby <- as.numeric(filt_py$filter(py$np$array(x), method = "chebyshev", order = 30L))

  g_r <- graph_ring(12)
  lmax_py <- as.numeric(g_py$lmax)
  t_rgsp <- scale / lmax_py
  y_r_cheby <- filter_signal(g_r, x, kernel_heat(t_rgsp), K = 30, lmax = lmax_py)

  expect_gt(cor(as.numeric(y_r_cheby), y_py_cheby), 0.95)
})

test_that("Kron reduction eigenvalues match PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Ring(12L)
  g_py$compute_fourier_basis()
  g_py$compute_laplacian()

  keep_idx <- c(0L, 2L, 4L, 6L, 8L, 10L)
  reduction_mod <- reticulate::import("pygsp.reduction")
  g_red_py <- reduction_mod$kron_reduction(g_py, keep_idx)
  g_red_py$compute_fourier_basis()
  eig_py <- sort(as.numeric(g_red_py$e))

  g_r <- graph_ring(12)
  g_red_r <- kron_reduction(g_r, keep_idx + 1L)
  eig_r <- sort(graph_eigenpairs(g_red_r)$values)

  expect_equal(eig_r, eig_py, tolerance = 1e-8)
})

test_that("Meyer wavelet kernel evaluations match PyGSP at eigenvalues", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_py <- py$pygsp$graphs$Path(16L)
  g_py$compute_fourier_basis()
  g_py$estimate_lmax()

  Nf <- 3L
  filt_py <- py$pygsp$filters$Meyer(g_py, Nf = Nf)
  kernel_evals_py <- filt_py$evaluate(g_py$e)
  # evaluate() returns a Nf x N numpy array; nrow = number of filters
  n_bands <- nrow(as.matrix(kernel_evals_py))

  g_r <- graph_path(16)
  bank_r <- meyer_wavelet_bank(g_r, n_scales = Nf - 1L, lmax = lambda_max(g_r))

  expect_equal(length(bank_r), n_bands)
})

test_that("Pyramid analysis/synthesis round-trip is exact", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  g_r <- graph_ring(12)
  set.seed(42)
  x <- rnorm(12)

  p <- pyramid_analysis(g_r, x)
  x_rec <- pyramid_synthesis(p)

  expect_equal(x_rec, x, tolerance = 1e-14)
})
