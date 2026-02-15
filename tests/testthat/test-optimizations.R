test_that("filter_signal_lanczos matches exact spectral filter (vector and matrix)", {
  set.seed(11)
  g <- graph_ring(20)
  L <- as.matrix(graph_laplacian(g, g$normalized))
  eig <- eigen(L, symmetric = TRUE)
  kernel <- function(lambda) exp(-0.4 * lambda)

  signal_vec <- rnorm(20)
  signal_mat <- cbind(rnorm(20), rnorm(20))

  apply_exact <- function(sig) {
    fLam <- kernel(eig$values)
    eig$vectors %*% (diag(fLam) %*% crossprod(eig$vectors, sig))
  }

  exact_vec <- drop(apply_exact(signal_vec))
  exact_mat <- apply_exact(signal_mat)

  res_vec <- filter_signal_lanczos(g, signal_vec, kernel = kernel, K = 15)
  res_mat <- filter_signal_lanczos(g, signal_mat, kernel = kernel, K = 15)

  expect_equal(res_vec, exact_vec, tolerance = 1e-5)
  expect_equal(res_mat, exact_mat, tolerance = 1e-5)
})

test_that("compute_frame Hutchinson estimator tracks exact frame bounds", {
  set.seed(21)
  g <- graph_ring(30)
  kernels <- list(kernel_heat(0.3), kernel_heat(0.8))

  exact <- compute_frame(g, kernels, method = "exact")
  approx <- compute_frame(g, kernels, method = "diag_hutch", R = 80, K = 40)

  max_abs_err <- max(abs(approx$diag - exact$diag))
  expect_lt(max_abs_err, 0.2)
  expect_equal(mean(approx$diag), mean(exact$diag), tolerance = 0.05)
})

test_that("isgwt iterative reconstruction matches original signal", {
  set.seed(33)
  g <- graph_ring(25)
  x <- rnorm(25)

  W <- sgwt(g, x, scales = c(2, 4, 8), include_lowpass = TRUE, K = 25)
  rec <- isgwt(W, method = "iterative", tol = 1e-8, maxiter = 200)

  expect_equal(rec, x, tolerance = 1e-4)
})
