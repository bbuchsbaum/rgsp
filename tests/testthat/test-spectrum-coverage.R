# Tests for spectrum.R coverage
# Covers: graph_eigenpairs, gft, igft, lambda_max, chebyshev_filter

test_that("graph_eigenpairs computes full eigendecomposition", {
  g <- graph_ring(10)
  eig <- graph_eigenpairs(g)

  expect_equal(length(eig$values), 10)
  expect_equal(dim(eig$vectors), c(10, 10))
  expect_true(all(is.finite(eig$values)))
  # Eigenvalues should be non-negative for Laplacian
  expect_true(all(eig$values >= -1e-10))
})

test_that("graph_eigenpairs caches results", {
  g <- graph_ring(10)
  eig1 <- graph_eigenpairs(g)
  eig2 <- graph_eigenpairs(g)
  expect_identical(eig1, eig2)
})

test_that("graph_eigenpairs with partial k", {
  g <- graph_ring(20)
  eig <- graph_eigenpairs(g, k = 5)
  expect_equal(length(eig$values), 5)
  expect_equal(ncol(eig$vectors), 5)
})

test_that("gft and igft are inverses", {
  set.seed(42)
  g <- graph_ring(15)
  signal <- rnorm(15)

  coeffs <- gft(g, signal)
  reconstructed <- igft(g, coeffs)

  expect_equal(drop(reconstructed), signal, tolerance = 1e-10)
})

test_that("gft works with matrix input", {
  set.seed(42)
  g <- graph_ring(10)
  X <- matrix(rnorm(10 * 3), nrow = 10)

  coeffs <- gft(g, X)
  expect_equal(dim(coeffs), c(10, 3))
})

test_that("gft errors on wrong signal length", {
  g <- graph_ring(10)
  expect_error(gft(g, rnorm(5)), "must match")
})

test_that("igft works with vector input", {
  g <- graph_ring(10)
  coeffs <- rnorm(10)
  result <- igft(g, coeffs)
  expect_equal(nrow(result), 10)
})

test_that("lambda_max returns positive value", {
  g <- graph_ring(15)
  lmax <- lambda_max(g)
  expect_true(lmax > 0)
  expect_true(is.finite(lmax))
})

test_that("lambda_max caches result", {
  g <- graph_ring(15)
  lmax1 <- lambda_max(g)
  lmax2 <- lambda_max(g)
  expect_identical(lmax1, lmax2)
})

test_that("chebyshev_filter works", {
  set.seed(42)
  g <- graph_ring(15)
  signal <- rnorm(15)
  kern <- kernel_heat(1)
  lmax <- lambda_max(g)
  coeffs <- cheby_coeffs(kern, K = 30, lmax = lmax)

  result <- chebyshev_filter(g, signal, coeffs, lambda_max_opt = lmax)
  expect_equal(length(result), 15)
  expect_true(all(is.finite(result)))
})

test_that("chebyshev_filter without precomputed lmax", {
  set.seed(42)
  g <- graph_ring(10)
  signal <- rnorm(10)
  kern <- kernel_heat(1)
  lmax <- lambda_max(g)
  coeffs <- cheby_coeffs(kern, K = 20, lmax = lmax)

  result <- chebyshev_filter(g, signal, coeffs)
  expect_equal(length(result), 10)
})
