test_that("grid3d torus degrees are 6 for dims >=3", {
  g <- graph_grid3d(3, 3, 3, periodic = TRUE)
  expect_equal(g$n, 27L)
  expect_equal(as.numeric(g$degree), rep(6, 27))
})
