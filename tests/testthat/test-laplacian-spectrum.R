test_that("laplacian eigenvalues are nonnegative", {
  g <- graph_grid2d(3, 3)
  eig <- graph_eigenpairs(g)
  expect_true(all(eig$values >= -1e-8))
})

test_that("lambda_max caches and matches eigs", {
  g <- graph_ring(12)
  l1 <- lambda_max(g)
  l2 <- lambda_max(g) # cached
  expect_equal(l1, l2)
  eig <- graph_eigenpairs(g)
  expect_equal(max(eig$values), l1, tolerance = 1e-6)
})

test_that("gft followed by igft returns original signal", {
  g <- graph_ring(10)
  x <- rnorm(10)
  X <- gft(g, x)
  xr <- igft(g, X)
  expect_equal(drop(xr), x, tolerance = 1e-8)
})

test_that("lambda_max respects normalized flag", {
  g <- graph_grid2d(2, 2, normalized = TRUE)
  ln <- lambda_max(g, normalized = TRUE)
  lu <- lambda_max(g, normalized = FALSE)
  expect_true(ln <= 2 + 1e-6)
  expect_true(lu >= ln)
})
