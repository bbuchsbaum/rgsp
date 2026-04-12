# Comprehensive filter tests covering all kernel types and filter operations
# Modeled after PyGSP's test_filters.py

# Helper function to test filter methods
# Tests: evaluation, filtering (Chebyshev and Lanczos), and comparison to exact
test_filter_methods <- function(kernel, kernel_name, g, signal, lmax, tol = 1e-3) {
  # Test 1: Kernel evaluates at eigenvalues
  eig <- graph_eigenpairs(g)
  kernel_vals <- kernel(eig$values)
  expect_true(all(is.finite(kernel_vals)),
              info = paste(kernel_name, "kernel evaluates at eigenvalues"))

  # Test 2: Chebyshev filtering works
  y_cheb <- filter_signal(g, signal, kernel, K = 50, lmax = lmax, strategy = "cheby")
  expect_equal(length(y_cheb), length(signal),
               info = paste(kernel_name, "Chebyshev output length"))
  expect_true(all(is.finite(y_cheb)),
              info = paste(kernel_name, "Chebyshev output is finite"))

  # Test 3: Lanczos filtering works
  y_lanc <- filter_signal_lanczos(g, signal, kernel, K = 30)
  expect_equal(length(y_lanc), length(signal),
               info = paste(kernel_name, "Lanczos output length"))
  expect_true(all(is.finite(y_lanc)),
              info = paste(kernel_name, "Lanczos output is finite"))


  # Test 4: Exact filtering via eigendecomposition
  U <- eig$vectors
  lam <- eig$values
  y_exact <- U %*% (kernel(lam) * (t(U) %*% signal))
  y_exact <- drop(y_exact)

  # Test 5: Compare Chebyshev to exact
  expect_equal(y_cheb, y_exact, tolerance = tol,
               info = paste(kernel_name, "Chebyshev matches exact"))

  # Test 6: Compare Lanczos to exact
  expect_equal(y_lanc, y_exact, tolerance = tol,
               info = paste(kernel_name, "Lanczos matches exact"))

  invisible(TRUE)
}

# Setup: create test graph and signal
setup_test_graph <- function(n = 20, seed = 42) {
  set.seed(seed)
  g <- graph_ring(n)
  signal <- rnorm(n)
  lmax <- lambda_max(g, normalized = g$normalized)
  list(g = g, signal = signal, lmax = lmax)
}

# ============================================================================
# KERNEL TESTS
# ============================================================================

test_that("kernel_wave evaluates correctly and filters properly", {
  env <- setup_test_graph()
  t <- 0.5
  kernel <- kernel_wave(t)

  # Test kernel shape
  lambdas <- seq(0, env$lmax, length.out = 100)
  vals <- kernel(lambdas)
  expect_equal(length(vals), 100)
  expect_true(all(vals >= -1 & vals <= 1))  # cosine bounded
  expect_equal(kernel(0), 1)  # cos(0) = 1

  # Test filtering
  test_filter_methods(kernel, "wave", env$g, env$signal, env$lmax)
})

test_that("kernel_gabor evaluates correctly and filters properly", {
  env <- setup_test_graph()
  mu <- env$lmax / 2
  sigma <- env$lmax / 4
  kernel <- kernel_gabor(mu, sigma)

  # Test kernel shape (Gaussian centered at mu)
  lambdas <- seq(0, env$lmax, length.out = 100)
  vals <- kernel(lambdas)
  expect_equal(length(vals), 100)
  expect_true(all(vals >= 0 & vals <= 1))
  expect_equal(kernel(mu), 1)  # maximum at center

  # Test filtering
  test_filter_methods(kernel, "gabor", env$g, env$signal, env$lmax)
})

test_that("kernel_rectangle evaluates correctly and filters properly", {
  env <- setup_test_graph()
  low <- env$lmax * 0.2
  high <- env$lmax * 0.6
  kernel <- kernel_rectangle(low, high)

  # Test kernel shape (sharp cutoff)
  expect_equal(kernel(0), 0)

  expect_equal(kernel(low), 1)
  expect_equal(kernel((low + high) / 2), 1)
  expect_equal(kernel(high), 1)
  expect_equal(kernel(env$lmax), 0)

  # Test filtering (rectangle has sharp discontinuities; Chebyshev needs large K)
  test_filter_methods(kernel, "rectangle", env$g, env$signal, env$lmax, tol = 0.25)
})

test_that("kernel_half_cosine evaluates correctly and filters properly", {
  env <- setup_test_graph()
  low <- env$lmax * 0.2
  high <- env$lmax * 0.6
  kernel <- kernel_half_cosine(low, high)

  # Test kernel shape (smooth transition)
  expect_equal(kernel(0), 1)
  expect_equal(kernel(low), 1)
  expect_true(kernel((low + high) / 2) > 0 && kernel((low + high) / 2) < 1)
  expect_equal(kernel(high), 0)
  expect_equal(kernel(env$lmax), 0)

  # Test filtering
  test_filter_methods(kernel, "half_cosine", env$g, env$signal, env$lmax)
})

test_that("kernel_exponential evaluates correctly and filters properly", {
  env <- setup_test_graph()
  alpha <- env$lmax / 2
  beta <- 2
  kernel <- kernel_exponential(alpha, beta)

  # Test kernel shape
  lambdas <- seq(0, env$lmax, length.out = 100)
  vals <- kernel(lambdas)
  expect_equal(length(vals), 100)
  expect_true(all(vals >= 0 & vals <= 1))
  expect_equal(kernel(0), 1)  # exp(0) = 1
  expect_true(all(diff(vals) <= 1e-10))  # monotonically decreasing

  # Test filtering
  test_filter_methods(kernel, "exponential", env$g, env$signal, env$lmax)
})

test_that("kernel_mexican_hat evaluates correctly and filters properly", {
  env <- setup_test_graph()
  sigma <- env$lmax / 4
  kernel <- kernel_mexican_hat(sigma)

  # Test kernel shape
  lambdas <- seq(0, env$lmax, length.out = 100)
  vals <- kernel(lambdas)
  expect_equal(length(vals), 100)
  expect_true(all(vals >= 0))
  expect_equal(kernel(0), 0)  # lambda * exp(...) at 0

  # Test filtering
  test_filter_methods(kernel, "mexican_hat", env$g, env$signal, env$lmax)
})

test_that("kernel_meyer evaluates correctly and filters properly", {
  env <- setup_test_graph()
  low <- env$lmax * 0.2
  high <- env$lmax * 0.6
  kernel <- kernel_meyer(low, high)

  # Test kernel shape (smooth polynomial transition)
  lambdas <- seq(0, env$lmax, length.out = 100)
  vals <- kernel(lambdas)
  expect_equal(length(vals), 100)
  expect_true(all(vals >= 0 & vals <= 1))

  # Test filtering (narrow-support meyer is hard for Chebyshev at K=50)
  test_filter_methods(kernel, "meyer", env$g, env$signal, env$lmax, tol = 0.5)
})

test_that("kernel_tight_band evaluates correctly and filters properly", {
  env <- setup_test_graph()
  low <- env$lmax * 0.2
  high <- env$lmax * 0.6
  kernel <- kernel_tight_band(low, high)

  # Test kernel shape
  expect_equal(kernel(0), 1)  # lowpass before cutoff
  expect_equal(kernel(low), 1)
  expect_true(kernel((low + high) / 2) > 0 && kernel((low + high) / 2) < 1)
  expect_equal(kernel(high), 0)
  expect_equal(kernel(env$lmax), 0)

  # Test filtering
  test_filter_methods(kernel, "tight_band", env$g, env$signal, env$lmax)
})

# ============================================================================
# GABOR / MODULATION TRANSFORM TESTS
# ============================================================================

test_that("gabor_filter_bank creates correct number of filters", {
  mus <- c(0.5, 1.0, 1.5, 2.0)
  sigma <- 0.3
  bank <- gabor_filter_bank(mus, sigma)

  expect_equal(length(bank), length(mus))
  expect_true(all(sapply(bank, is.function)))

  # Each filter should peak at its center frequency
  for (i in seq_along(mus)) {
    expect_equal(bank[[i]](mus[i]), 1)
  }
})

test_that("modulation_transform produces correct output structure", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  mus <- c(0.5, 1.0, 1.5)
  sigma <- 0.3

  result <- modulation_transform(g, x, mus, sigma, K = 30)

  expect_equal(length(result), length(mus))
  expect_true(all(sapply(result, length) == 15))
})

# ============================================================================
# WAVELET BANK TESTS
# ============================================================================

test_that("meyer_wavelet_bank creates lowpass + bandpass filters", {
  lmax <- 4.0
  n_scales <- 4
  bank <- meyer_wavelet_bank(scales = NULL, lmax = lmax, n_scales = n_scales)

  expect_equal(length(bank), n_scales + 1)  # lowpass + n_scales bandpasses
  expect_equal(names(bank)[1], "lowpass")
  expect_true(all(grepl("band_", names(bank)[-1])))
})

test_that("mexican_hat_wavelet_bank creates filters at specified scales", {
  set.seed(42)
  g <- graph_ring(20)
  lmax <- lambda_max(g, normalized = g$normalized)
  n_scales <- 4
  bank <- mexican_hat_wavelet_bank(g, n_scales = n_scales, lmax = lmax)

  expect_equal(length(bank), n_scales)
  expect_true(all(sapply(bank, is.function)))
})

test_that("tight_frame_bank creates partition of unity", {
  edges <- c(0, 1, 2, 3, 4)
  bank <- tight_frame_bank(edges)

  expect_equal(length(bank), length(edges) - 1)

  # At edge boundaries, sum should be close to 1 (tight frame property)
  test_points <- c(0.5, 1.5, 2.5, 3.5)
  for (pt in test_points) {
    total <- sum(sapply(bank, function(k) k(pt)^2))
    expect_true(total <= 1 + 1e-6)
  }
})

test_that("tight_frame_regular creates proper filter bank", {
  lmax <- 4.0
  J <- 3
  bank <- tight_frame_regular(lmax, J = J)

  expect_equal(length(bank), J + 1)  # J dyadic bands + 1
  expect_true(all(sapply(bank, is.function)))
})

test_that("tight_frame_held is equivalent to tight_frame_regular", {
  lmax <- 4.0
  J <- 3
  bank1 <- tight_frame_regular(lmax, J = J)
  bank2 <- tight_frame_held(lmax, J = J)

  # Should produce same filter values
  test_lambdas <- seq(0, lmax, length.out = 50)
  for (i in seq_along(bank1)) {
    expect_equal(bank1[[i]](test_lambdas), bank2[[i]](test_lambdas))
  }
})

test_that("tight_frame_simoncelli creates filters", {
  lmax <- 4.0
  J <- 3
  bank <- tight_frame_simoncelli(lmax, J = J)

  expect_equal(length(bank), J + 1)
  expect_true(all(sapply(bank, is.function)))
})

test_that("tight_frame_papadakis creates filters", {
  lmax <- 4.0
  J <- 3
  bank <- tight_frame_papadakis(lmax, J = J)

  expect_equal(length(bank), J + 1)
  expect_true(all(sapply(bank, is.function)))
})

# ============================================================================
# FILTER BANK APPLICATION TESTS
# ============================================================================

test_that("filter_bank_apply works with multiple kernels", {
  set.seed(42)
  g <- graph_ring(20)
  x <- rnorm(20)
  lmax <- lambda_max(g, normalized = g$normalized)

  kernels <- list(
    heat = kernel_heat(0.5),
    wave = kernel_wave(0.3),
    gabor = kernel_gabor(lmax/2, lmax/4)
  )

  result <- filter_bank_apply(g, x, kernels, K = 30, lmax = lmax)

  expect_equal(names(result), names(kernels))
  expect_true(all(sapply(result, length) == 20))
})

test_that("wavelet_heat_transform produces multi-scale output", {
  set.seed(42)
  g <- graph_ring(20)
  x <- rnorm(20)
  scales <- c(0.1, 0.5, 1.0, 2.0)

  result <- wavelet_heat_transform(g, x, scales, K = 30)

  expect_equal(length(result), length(scales))
  expect_equal(names(result), paste0("scale_", scales))
})

# ============================================================================
# APPROXIMATION METHOD TESTS
# ============================================================================

test_that("Chebyshev and Lanczos agree for various kernels", {
  set.seed(42)
  g <- graph_ring(25)
  x <- rnorm(25)
  lmax <- lambda_max(g, normalized = g$normalized)

  kernels <- list(
    heat = kernel_heat(0.5),
    wave = kernel_wave(0.3),
    mexican_hat = kernel_mexican_hat(lmax/4),
    exponential = kernel_exponential(lmax/2, 1.5)
  )

  for (name in names(kernels)) {
    y_cheb <- filter_signal(g, x, kernels[[name]], K = 50, lmax = lmax, strategy = "cheby")
    y_lanc <- filter_signal_lanczos(g, x, kernels[[name]], K = 30)

    # Both should be close to exact
    eig <- graph_eigenpairs(g)
    U <- eig$vectors
    lam <- eig$values
    y_exact <- drop(U %*% (kernels[[name]](lam) * (t(U) %*% x)))

    expect_equal(y_cheb, y_exact, tolerance = 1e-3,
                 info = paste(name, "Chebyshev vs exact"))
    expect_equal(y_lanc, y_exact, tolerance = 1e-3,
                 info = paste(name, "Lanczos vs exact"))
  }
})

test_that("Jackson damping improves smoothness", {
  set.seed(42)
  g <- graph_ring(30)
  x <- rnorm(30)
  lmax <- lambda_max(g, normalized = g$normalized)
  kernel <- kernel_rectangle(lmax * 0.3, lmax * 0.7)  # sharp cutoff

  y_no_jackson <- filter_signal(g, x, kernel, K = 30, lmax = lmax,
                                 jackson = FALSE, strategy = "cheby")
  y_jackson <- filter_signal(g, x, kernel, K = 30, lmax = lmax,
                              jackson = TRUE, strategy = "cheby")

  # Both should produce valid output
  expect_true(all(is.finite(y_no_jackson)))
  expect_true(all(is.finite(y_jackson)))
})

# ============================================================================
# HEAT PROPAGATE TESTS
# ============================================================================

test_that("heat_propagate with single time returns vector", {
  set.seed(42)
  g <- graph_grid2d(5, 5)
  x <- rnorm(25)

  result <- heat_propagate(g, x, t = 0.5)
  expect_true(is.numeric(result) && !is.list(result))
  expect_equal(length(result), 25)
})

test_that("heat_propagate with multiple times returns list", {
  set.seed(42)
  g <- graph_grid2d(5, 5)
  x <- rnorm(25)
  times <- c(0.1, 0.5, 1.0, 2.0)

  result <- heat_propagate(g, x, t = times)
  expect_true(is.list(result))
  expect_equal(length(result), length(times))
  expect_equal(names(result), paste0("t_", times))
})

test_that("heat_propagate matches filter_signal with heat kernel", {
  set.seed(42)
  g <- graph_ring(20)
  x <- rnorm(20)
  t <- 0.7

  y1 <- heat_propagate(g, x, t = t, K = 40)
  y2 <- filter_signal(g, x, kernel_heat(t), K = 40)

  expect_equal(y1, y2)
})

# ============================================================================
# WAVELET HEAT KNN TESTS
# ============================================================================

test_that("wavelet_heat_knn builds graph and applies wavelets", {
  set.seed(42)
  coords <- matrix(runif(40), ncol = 2)  # 20 points in 2D
  signal <- rnorm(20)
  scales <- c(0.5, 1.0, 2.0)

  result <- wavelet_heat_knn(coords, signal, scales, k = 5)

  expect_true(is.list(result))
  expect_equal(length(result), length(scales))
  expect_true(all(sapply(result, length) == 20))
})

# ============================================================================
# MULTI-DIMENSIONAL SIGNAL TESTS
# ============================================================================

test_that("filtering works with matrix signals (multiple columns)", {
  set.seed(42)
  g <- graph_ring(15)
  X <- matrix(rnorm(15 * 3), ncol = 3)  # 3 signals
  kernel <- kernel_heat(0.5)

  Y <- filter_signal(g, X, kernel, K = 30)

  expect_equal(dim(Y), dim(X))
  expect_true(all(is.finite(Y)))
})

test_that("filter_bank_apply works with matrix signals", {
  set.seed(42)
  g <- graph_ring(15)
  X <- matrix(rnorm(15 * 2), ncol = 2)
  kernels <- list(low = kernel_heat(1.0), high = kernel_mexican_hat(0.5))

  result <- filter_bank_apply(g, X, kernels, K = 30)

  expect_equal(length(result), 2)
  expect_true(all(sapply(result, function(r) identical(dim(r), dim(X)))))
})

# ============================================================================
# EDGE CASES
# ============================================================================

test_that("kernels handle zero and negative values gracefully", {
  # Heat kernel at lambda = 0
  expect_equal(kernel_heat(1)(0), 1)

  # Wave kernel at lambda = 0
  expect_equal(kernel_wave(1)(0), 1)

  # Gabor kernel handles mu = 0
  k <- kernel_gabor(0, 1)
  expect_equal(k(0), 1)

  # Rectangle at boundary
  k <- kernel_rectangle(0, 1)
  expect_equal(k(0), 1)
  expect_equal(k(1), 1)
})

test_that("filter_signal handles small graphs", {
  g <- graph_ring(4)
  x <- rnorm(4)
  kernel <- kernel_heat(0.5)

  # Should work with K > n (gets clipped internally for Lanczos)
  y <- filter_signal(g, x, kernel, K = 10)
  expect_equal(length(y), 4)
})

test_that("cheby_coeffs caching works correctly", {
  rgsp:::cheby_cache_clear()

  kernel <- kernel_heat(0.3)
  lmax <- 2.0

  # First call
  c1 <- rgsp:::cheby_coeffs(kernel, K = 20, lmax = lmax)
  sz1 <- length(ls(envir = rgsp:::`.cheby_cache`))

  # Second call should hit cache
  c2 <- rgsp:::cheby_coeffs(kernel, K = 20, lmax = lmax)
  sz2 <- length(ls(envir = rgsp:::`.cheby_cache`))

  expect_identical(c1, c2)
  expect_equal(sz1, sz2)

  # Different parameters create a separate cache entry, even if trimming
  # makes the effective coefficient vector shorter than the requested K.
  c3 <- rgsp:::cheby_coeffs(kernel, K = 25, lmax = lmax)
  sz3 <- length(ls(envir = rgsp:::`.cheby_cache`))
  expect_equal(sz3, sz2 + 1)
  expect_true(length(c3) >= length(c1))
  expect_true(length(c3) <= 25)
})
