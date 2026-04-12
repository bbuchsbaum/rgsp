bench_median_ms <- function(fun, reps = 3L, warmup = 1L) {
  for (i in seq_len(warmup)) fun()
  median(replicate(reps, {
    gc()
    system.time(fun())[["elapsed"]]
  })) * 1000
}

cheby_reference_r <- function(L, X, coeffs, lmax) {
  K <- length(coeffs)
  s <- 2 / lmax
  Y <- coeffs[[1]] * X
  if (K == 1) return(Y)

  Z <- L %*% X
  T_prev <- X
  T_curr <- s * Z - X
  Y <- Y + coeffs[[2]] * T_curr
  if (K == 2) return(Y)

  for (k in 3:K) {
    Z <- L %*% T_curr
    Z <- (2 * s) * Z - 2 * T_curr - T_prev
    Y <- Y + coeffs[[k]] * Z
    T_prev <- T_curr
    T_curr <- Z
  }

  Y
}

test_that("gft cached path stays faster than rebuilding transpose", {
  skip_on_cran()

  set.seed(1)
  n <- 1024
  n_signals <- 8
  g <- graph_ring(n)
  X <- matrix(rnorm(n * n_signals), nrow = n, ncol = n_signals)

  eig <- graph_eigenpairs(g)
  invisible(gft(g, X))
  U <- eig$vectors

  t_cached <- bench_median_ms(function() gft(g, X))
  t_raw <- bench_median_ms(function() t(U) %*% X)

  expect_lte(t_cached / t_raw, 0.9)
})

test_that("Chebyshev filter stays faster than the pure-R recurrence baseline", {
  skip_on_cran()

  set.seed(1)
  n <- 1024
  n_signals <- 8
  lmax <- 4
  g <- graph_ring(n)
  X <- matrix(rnorm(n * n_signals), nrow = n, ncol = n_signals)
  L <- graph_laplacian(g, g$normalized)
  coeffs <- cheby_coeffs(kernel_heat(0.25), K = 30, lmax = lmax)

  t_cpp <- bench_median_ms(function() {
    filter_signal(g, X, kernel_heat(0.25), K = 30, lmax = lmax, strategy = "cheby")
  })
  t_ref <- bench_median_ms(function() cheby_reference_r(L, X, coeffs, lmax))

  expect_lte(t_cpp / t_ref, 0.8)
})
