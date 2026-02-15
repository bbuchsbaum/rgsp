test_that("kron_reduction shrinks node count and preserves symmetry", {
  g <- graph_grid2d(3, 3)
  keep <- c(1, 3, 5, 7, 9)
  g2 <- kron_reduction(g, keep)
  expect_equal(g2$n, length(keep))
  expect_true(Matrix::isSymmetric(g2$adjacency))
})

test_that("graph_multiresolution returns hierarchy", {
  g <- graph_ring(12)
  h <- graph_multiresolution(g, levels = 2, keep_fraction = 0.5, seed = 1)
  expect_equal(length(h$graphs), 3)
  expect_true(h$graphs[[2]]$n <= g$n)
})

test_that("pyramid analysis/synthesis round-trips signal", {
  g <- graph_ring(6)
  x <- rnorm(6)
  p <- pyramid_analysis(g, x, keep = c(1, 3, 5))
  xr <- pyramid_synthesis(p)
  expect_equal(xr, x)
})
