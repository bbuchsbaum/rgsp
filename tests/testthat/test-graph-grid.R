test_that("grid2d torus degrees are correct", {
  g <- graph_grid2d(3, 4, periodic = TRUE)
  expect_equal(g$n, 12L)
  expect_true(Matrix::isSymmetric(g$adjacency))
  # torus: each node degree 4
  expect_equal(as.numeric(g$degree), rep(4, 12))
})

test_that("grid2d with no wrap has boundary degrees", {
  g <- graph_grid2d(2, 2, periodic = FALSE)
  expect_equal(sort(as.numeric(g$degree)), c(2,2,2,2)) # small grid all degree 2
})
