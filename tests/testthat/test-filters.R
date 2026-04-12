test_that("heat filter via Chebyshev matches eigen filtering", {
  set.seed(42)
  g <- graph_ring(8)
  x <- rnorm(8)
  t <- 0.3
  K <- 40
  # exact via eigen
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  lam <- eig$values
  y_exact <- U %*% (exp(-t * lam) * (t(U) %*% x))
  y_exact <- drop(y_exact)

  y_cheb <- filter_signal(g, x, kernel_heat(t), K = K)
  expect_equal(y_cheb, y_exact, tolerance = 1e-4)
})

test_that("heat filter via Chebyshev matches eigen filtering for matrix signals", {
  set.seed(42)
  g <- graph_ring(128)
  X <- matrix(rnorm(128 * 8), nrow = 128, ncol = 8)
  t <- 0.25
  K <- 40
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  lam <- eig$values
  H <- diag(exp(-t * lam))
  y_exact <- U %*% H %*% t(U) %*% X

  y_cheb <- filter_signal(g, X, kernel_heat(t), K = K, strategy = "cheby")
  expect_equal(y_cheb, y_exact, tolerance = 1e-4)
})

test_that("exact ring filtering matches eigen filtering for matrix signals", {
  set.seed(42)
  g <- graph_ring(128)
  X <- matrix(rnorm(128 * 8), nrow = 128, ncol = 8)
  t <- 0.25
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  lam <- eig$values
  H <- diag(exp(-t * lam))
  y_exact <- U %*% H %*% t(U) %*% X

  y_ring <- filter_signal_exact(g, X, kernel_heat(t))
  expect_equal(y_ring, y_exact, tolerance = 1e-8)
})

test_that("auto strategy on ring uses the exact path", {
  set.seed(42)
  g <- graph_ring(64)
  X <- matrix(rnorm(64 * 4), nrow = 64, ncol = 4)
  t <- 0.25

  y_auto <- filter_signal(g, X, kernel_heat(t), strategy = "auto")
  y_exact <- filter_signal(g, X, kernel_heat(t), strategy = "exact")

  expect_equal(y_auto, y_exact, tolerance = 1e-10)
})

test_that("exact path filtering matches eigen filtering for matrix signals", {
  set.seed(42)
  g <- graph_path(96)
  X <- matrix(rnorm(96 * 6), nrow = 96, ncol = 6)
  t <- 0.2
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  lam <- eig$values
  H <- diag(exp(-t * lam))
  y_exact <- U %*% H %*% t(U) %*% X

  y_path <- filter_signal_exact(g, X, kernel_heat(t))
  expect_equal(y_path, y_exact, tolerance = 1e-8)
})

test_that("exact grid2d filtering matches eigen filtering for matrix signals", {
  set.seed(42)
  g <- graph_grid2d(6, 5)
  X <- matrix(rnorm(30 * 4), nrow = 30, ncol = 4)
  t <- 0.2
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  lam <- eig$values
  H <- diag(exp(-t * lam))
  y_exact <- U %*% H %*% t(U) %*% X

  y_grid <- filter_signal_exact(g, X, kernel_heat(t))
  expect_equal(y_grid, y_exact, tolerance = 1e-8)
})

test_that("auto strategy on path and grid2d uses the exact path", {
  set.seed(42)

  g_path <- graph_path(512)
  X_path <- matrix(rnorm(512 * 3), nrow = 512, ncol = 3)
  expect_equal(
    filter_signal(g_path, X_path, kernel_heat(0.2), strategy = "auto"),
    filter_signal(g_path, X_path, kernel_heat(0.2), strategy = "exact"),
    tolerance = 1e-10
  )

  g_grid <- graph_grid2d(4, 5)
  X_grid <- matrix(rnorm(20 * 3), nrow = 20, ncol = 3)
  expect_equal(
    filter_signal(g_grid, X_grid, kernel_heat(0.2), strategy = "auto"),
    filter_signal(g_grid, X_grid, kernel_heat(0.2), strategy = "exact"),
    tolerance = 1e-10
  )
})

test_that("heat filter via Lanczos matches eigen filtering", {
  set.seed(42)
  g <- graph_ring(12)
  x <- rnorm(12)
  t <- 0.3
  K <- 20
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  lam <- eig$values
  y_exact <- U %*% (exp(-t * lam) * (t(U) %*% x))
  y_exact <- drop(y_exact)

  y_lanc <- filter_signal_lanczos(g, x, kernel_heat(t), K = K)
  expect_equal(y_lanc, y_exact, tolerance = 1e-4)
})

test_that("Chebyshev coefficients are cached and reused", {
  rgsp:::cheby_cache_clear()
  g <- graph_ring(10)
  x <- rnorm(10)
  kern <- kernel_heat(0.2)
  # first call populates cache
  filter_signal(g, x, kern, K = 15)
  sz1 <- length(ls(envir = rgsp:::`.cheby_cache`))
  # second call should not increase cache size
  filter_signal(g, x, kern, K = 15)
  sz2 <- length(ls(envir = rgsp:::`.cheby_cache`))
  expect_equal(sz1, sz2)
})

test_that("Jackson damping applies monotone weights", {
  K <- 20
  w <- rgsp:::jackson_weights(K)
  expect_equal(length(w), K)
  expect_true(all(w >= 0))
  expect_true(max(w) <= 1)
  expect_true(all(diff(w[-1]) <= 1e-12)) # non-increasing after k=1
})

test_that("Chebyshev coefficients trim negligible heat terms", {
  coeffs <- rgsp:::cheby_coeffs(kernel_heat(0.25), K = 30, lmax = 4, cache = FALSE)
  expect_lt(length(coeffs), 30)
  expect_equal(length(coeffs), 10)
})

test_that("filter bank returns list with correct lengths", {
  g <- graph_ring(5)
  x <- rnorm(5)
  kernels <- list(low = kernel_rectangle(0, 1), high = kernel_rectangle(1, 4))
  out <- filter_bank_apply(g, x, kernels, K = 20)
  expect_equal(names(out), c("low", "high"))
  expect_equal(length(out$low), 5)
})

test_that("tight frame bank covers expected bands", {
  edges <- c(0, 1, 2)
  bank <- tight_frame_bank(edges)
  expect_equal(length(bank), 2)
  l <- c(0.5, 1.5, 2.5)
  vals <- sapply(bank, function(k) k(l))
  expect_true(all(vals[1, ] >= 0))
})
