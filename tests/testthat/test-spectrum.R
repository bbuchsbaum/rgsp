test_that("GFT followed by IGFT recovers the signal", {
  g <- graph_ring(6)
  x <- rnorm(6)
  cfs <- gft(g, x)
  xr <- igft(g, cfs)
  expect_equal(drop(xr), x, tolerance = 1e-8)
})

test_that("Chebyshev filter matches manual recurrence", {
  g <- graph_ring(5)
  L <- graph_laplacian(g)
  lmax <- lambda_max(g)
  x <- matrix(rnorm(5), ncol = 1)
  coeffs <- c(1, 0.1, -0.05)

  # manual recurrence in R
  Ltilde <- (2 / lmax) * L - Matrix::Diagonal(5)
  T0 <- x
  T1 <- Ltilde %*% x
  T2 <- 2 * Ltilde %*% T1 - T0
  y_manual <- coeffs[1] * T0 + coeffs[2] * T1 + coeffs[3] * T2
  y_manual <- as.numeric(y_manual)

  y_cpp <- chebyshev_filter(g, x, coeffs, lambda_max_opt = lmax)
  expect_equal(drop(y_cpp), y_manual, tolerance = 1e-8)
})
