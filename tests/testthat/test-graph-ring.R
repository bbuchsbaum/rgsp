test_that("ring graph builds and degrees are correct", {
  g <- graph_ring(5)
  expect_s3_class(g, "gsp_graph")
  expect_equal(g$n, 5L)
  expect_equal(length(g$degree), 5L)
  expect_equal(as.numeric(g$degree), rep(2, 5))
  expect_true(Matrix::isSymmetric(g$adjacency))
})

test_that("path graph builds and end nodes have degree 1", {
  g <- graph_path(4)
  expect_equal(as.numeric(g$degree), c(1, 2, 2, 1))
})

test_that("random geometric graph builds", {
  set.seed(1)
  g <- graph_random_geometric(20, radius = 0.3)
  expect_s3_class(g, "gsp_graph")
  expect_equal(g$n, 20L)
})

test_that("rescale invalidates caches", {
  g <- graph_ring(4)
  g$lmax_pre <- lambda_max(g)
  g2 <- graph_rescale(g, 2)
  expect_false(identical(g2$cache$lmax_TRUE, g$lmax_pre))
})
