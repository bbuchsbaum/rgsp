test_that("sgwt followed by isgwt reconstructs ring signal (tightish frame)", {
  set.seed(123)
  g <- graph_ring(10)
  x <- rnorm(10)
  scales <- c(1, 2, 4)
  W <- sgwt(g, x, scales = scales, include_lowpass = TRUE, K = 30, wavelet = "heat")
  xr <- isgwt(W, method = "auto", tol = 1e-6, maxiter = 200)
  expect_equal(xr, x, tolerance = 2e-1)
})

test_that("sgwt handles matrix input and preserves shape", {
  set.seed(2)
  g <- graph_sensor(12, k = 3, seed = 7)
  X <- matrix(rnorm(12 * 3), nrow = 12, ncol = 3)
  W <- sgwt(g, X, scales = c(0.5, 1, 2), include_lowpass = FALSE, K = 20)
  coeffs <- W$coeffs
  expect_equal(dim(coeffs), c(12, 3, 3))
  expect_false(W$include_lowpass)
})

test_that("sgwt_coeffs_at_scale extracts correct band", {
  set.seed(4)
  g <- graph_ring(8)
  x <- rnorm(8)
  scales <- c(1, 3)
  W <- sgwt(g, x, scales = scales, include_lowpass = TRUE, K = 25)
  low <- sgwt_coeffs_at_scale(W, "lowpass")
  band1 <- sgwt_coeffs_at_scale(W, scales[1])
  expect_equal(length(low), 8)
  expect_equal(length(band1), 8)
})

test_that("compute_spectrogram handles multiple signals", {
  set.seed(5)
  g <- graph_ring(6)
  X <- matrix(rnorm(6 * 4), nrow = 6, ncol = 4)
  S <- compute_spectrogram(g, X, scales = c(0.5, 1.5), K = 15, jackson = TRUE)
  expect_equal(dim(S), c(6, 2, 4))
  expect_true(all(S >= 0))
})

test_that("compute_frame exact A,B nonnegative and ordered", {
  g <- graph_grid2d(3, 3)
  lmax <- lambda_max(g)
  kernels <- meyer_wavelet_bank(scales = c(1, 2), lmax = lmax)
  cf <- compute_frame(g, kernels, method = "exact")
  expect_true(cf$A >= 0)
  expect_true(cf$B >= cf$A)
})

test_that("graph_incidence with oriented FALSE returns double edges", {
  g <- graph_ring(4)
  inc <- graph_incidence(g, oriented = FALSE)
  inc_or <- graph_incidence(g, oriented = TRUE)
  expect_true(nrow(inc$B) >= 2 * nrow(inc_or$B))
})

test_that("grad/div with multiple signals round trip to near-zero", {
  set.seed(6)
  g <- graph_ring(7)
  X <- matrix(rnorm(7 * 3), nrow = 7)
  G <- grad(g, X)
  D <- div(g, G)
  expect_true(max(abs(Matrix::colSums(as.matrix(D)))) < 1e-6)
})

test_that("graph_random_geometric uses coords of correct size", {
  g <- graph_random_geometric(15, radius = 0.5, seed = 1)
  expect_equal(dim(g$coords), c(15, 2))
})

test_that("graph_grid3d coordinates ordered correctly", {
  g <- graph_grid3d(2, 2, 2)
  expect_equal(dim(g$coords), c(8, 3))
  expect_equal(unname(g$coords[1, ]), c(1, 1, 1))
})

test_that("graph_complete directed has correct edge count", {
  g <- graph_complete(5, directed = TRUE)
  # for directed no self loops: n*(n-1)
  expect_equal(length(g$adjacency@x), 5 * 4)
})

test_that("filter_bank_apply handles jackson damping flag", {
  g <- graph_ring(6)
  x <- rnorm(6)
  kernels <- list(low = kernel_rectangle(0, 1), high = kernel_rectangle(1, 3))
  out1 <- filter_bank_apply(g, x, kernels, K = 15, jackson = FALSE)
  out2 <- filter_bank_apply(g, x, kernels, K = 15, jackson = TRUE)
  expect_equal(length(out1), length(out2))
  expect_equal(names(out1), names(out2))
})

test_that("cheby_cache_clear empties cache", {
  rgsp:::cheby_cache_clear()
  g <- graph_ring(5)
  x <- rnorm(5)
  filter_signal(g, x, kernel_heat(0.3), K = 10)
  expect_true(length(ls(envir = rgsp:::`.cheby_cache`)) >= 1)
  cheby_cache_clear()
  expect_true(length(ls(envir = rgsp:::`.cheby_cache`)) == 0)
})

test_that("graph_sbm probabilities bounded", {
  block_sizes <- c(4, 4, 4)
  P <- matrix(c(0.2, 0.1, 0.05,
                0.1, 0.3, 0.02,
                0.05, 0.02, 0.25), 3, 3, byrow = TRUE)
  g <- graph_sbm(block_sizes, P, seed = 11)
  expect_true(Matrix::isSymmetric(g$adjacency))
  expect_equal(Matrix::diag(g$adjacency), rep(0, g$n))
})

test_that("grad/div agree with laplacian for scalar signal", {
  g <- graph_ring(9)
  x <- rnorm(9)
  gx <- grad(g, x)
  dx <- div(g, gx)
  Lx <- graph_laplacian(g) %*% x
  expect_equal(as.numeric(dx), -as.numeric(Lx), tolerance = 1e-6)
})

test_that("lambda_max cached invalidates after rescale", {
  g <- graph_ring(8)
  l1 <- lambda_max(g)
  g2 <- graph_rescale(g, 2)
  l2 <- lambda_max(g2)
  expect_true(l2 > l1 * 1.9)
})

test_that("heat_propagate with vector and multi-time returns expected shapes", {
  g <- graph_ring(6)
  x <- rnorm(6)
  out1 <- heat_propagate(g, x, t = 0.2)
  out2 <- heat_propagate(g, x, t = c(0.1, 0.3))
  expect_equal(length(out1), 6)
  expect_equal(length(out2), 2)
  expect_true(all(names(out2) == c("t_0.1", "t_0.3")))
})

test_that("filter_signal_lanczos handles matrix signals", {
  g <- graph_ring(5)
  X <- matrix(rnorm(10), nrow = 5)
  out <- filter_signal_lanczos(g, X, kernel_heat(0.4), K = 10)
  expect_equal(dim(out), dim(X))
})

test_that("wavelet_heat_transform preserves column count", {
  g <- graph_ring(7)
  X <- matrix(rnorm(14), nrow = 7, ncol = 2)
  res <- wavelet_heat_transform(g, X, scales = c(0.5, 1), K = 15)
  expect_equal(length(res), 2)
  expect_equal(dim(res[[1]]), dim(X))
})

test_that("compute_frame Hutchinson gets tighter with more probes", {
  g <- graph_ring(10)
  kernels <- mexican_hat_wavelet_bank(scales = c(1, 2, 4))
  cf10 <- compute_frame(g, kernels, method = "diag_hutch", R = 8, K = 20)
  cf50 <- compute_frame(g, kernels, method = "diag_hutch", R = 50, K = 20)
  expect_true(cf50$B - cf50$A <= cf10$B - cf10$A + 1e-6)
})
