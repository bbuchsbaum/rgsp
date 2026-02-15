test_that("random_lowpass_basis returns orthonormal loadings", {
  set.seed(123)
  g <- graph_ring(20)
  res <- random_lowpass_basis(g, k = 5, K = 20, oversample = 5, seed = 123)

  U <- res$loadings
  expect_equal(dim(U), c(20, 5))

  gram <- crossprod(U)
  expect_lt(max(abs(gram - diag(5))), 1e-6)

  expect_length(res$lambdas, 5)
})

test_that("graph_compress projects data to correct shape", {
  set.seed(42)
  g <- graph_ring(30)
  X <- matrix(rnorm(20 * 30), nrow = 20, ncol = 30) # rows = samples, cols = nodes

  res <- graph_compress(g, X, k = 6, K = 15, oversample = 4, seed = 99)

  expect_equal(dim(res$basis), c(20, 6))
  expect_equal(dim(res$loadings), c(30, 6))

  X_hat <- res$basis %*% t(res$loadings)
  rel_err <- norm(X - X_hat, "F") / norm(X, "F")
  expect_lt(rel_err, 1)  # projection should not increase energy
})

test_that("graph_compress block_rows matches full computation", {
  set.seed(7)
  g <- graph_ring(15)
  X <- matrix(rnorm(40 * 15), nrow = 40, ncol = 15)

  res_full <- graph_compress(g, X, k = 4, K = 12, oversample = 3, seed = 7)
  res_block <- graph_compress(g, X, k = 4, K = 12, oversample = 3, seed = 7, block_rows = 10)

  diff <- max(abs(res_full$basis - res_block$basis))
  expect_lt(diff, 1e-10)
})

test_that("default_probe_block_cols matches heuristic", {
  expect_equal(default_probe_block_cols(1e5, 1000), 375L)
  expect_equal(default_probe_block_cols(10, 18, target_mb = 300), 18L)
})

test_that("random_lowpass_basis logs chosen block when verbose", {
  g <- graph_ring(10)
  expect_message(
    random_lowpass_basis(g, k = 2, K = 5, jackson = TRUE, verbose = TRUE, seed = 1),
    "block="
  )
})
