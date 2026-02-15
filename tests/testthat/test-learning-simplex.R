test_that("classification_tikhonov_simplex returns simplex logits", {
  g <- graph_ring(30)

  # 3-class labels with missing values
  y <- rep(NA_integer_, g$n)
  y[c(1, 6, 11)] <- 1L
  y[c(3, 8, 13)] <- 2L
  y[c(5, 10, 15)] <- 3L
  mask <- !is.na(y)

  X <- classification_tikhonov_simplex(g, y, mask = mask, tau = 0.2, maxit = 200, tol = 1e-7)

  expect_true(is.matrix(X))
  expect_equal(nrow(X), g$n)
  expect_equal(ncol(X), 3L)

  # Simplex constraints
  expect_equal(rowSums(X), rep(1, g$n), tolerance = 1e-6)
  expect_true(all(X >= -1e-12))
})

