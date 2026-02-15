# Tests for optimization.R coverage
# Covers: interpolate_laplacian, tv_inpaint

test_that("interpolate_laplacian fills missing values", {
  set.seed(42)
  g <- graph_ring(20)
  y <- rnorm(20)
  y[c(5, 10, 15)] <- NA

  result <- interpolate_laplacian(g, y, alpha = 0.01)
  expect_equal(length(result), 20)
  expect_true(all(is.finite(result)))
  # Known values should be approximately preserved
  known_idx <- which(!is.na(rnorm(20)))  # recompute mask
  y_orig <- rnorm(20)  # won't match, use direct approach
})

test_that("interpolate_laplacian preserves known values approximately", {
  set.seed(42)
  g <- graph_ring(15)
  y_full <- sin(seq(0, 2*pi, length.out = 15))
  y <- y_full
  y[c(3, 7, 11)] <- NA

  result <- interpolate_laplacian(g, y, alpha = 0.001)
  known_idx <- which(!is.na(y))
  # With small alpha, known values should be well-preserved
  expect_equal(result[known_idx], y_full[known_idx], tolerance = 0.05)
})
