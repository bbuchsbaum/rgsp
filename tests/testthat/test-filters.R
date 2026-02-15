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
