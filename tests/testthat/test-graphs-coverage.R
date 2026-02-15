# Tests for graphs.R coverage gaps
# Covers: graph_degree, graph_laplacian, print.gsp_graph, graph_cache_clear,
# graph_gradient, knn_graph, graph_sbm, graph_random_regular, graph_barabasi_albert

test_that("graph_degree computes degrees", {
  g <- graph_ring(10)
  d <- graph_degree(g)
  expect_equal(length(d), 10)
  expect_true(all(d == 2))  # ring graph has degree 2
})

test_that("graph_laplacian returns sparse matrix", {
  g <- graph_ring(10)
  L <- graph_laplacian(g)
  expect_true(inherits(L, "Matrix") || is.matrix(L))
  expect_equal(nrow(L), 10)
  expect_equal(ncol(L), 10)
})

test_that("print.gsp_graph works", {
  g <- graph_ring(10)
  out <- capture.output(print(g))
  expect_true(length(out) > 0)
})

test_that("graph_cache_clear removes cached values", {
  g <- graph_ring(10)
  lambda_max(g)  # populate cache
  g2 <- graph_cache_clear(g)
  expect_equal(length(g2$cache), 0)
})

test_that("graph_gradient computes edge differences", {
  g <- graph_ring(10)
  D <- graph_gradient(g)
  expect_true(nrow(D) > 0)
  expect_equal(ncol(D), 10)
})

test_that("knn_graph builds graph from coordinates", {
  set.seed(42)
  coords <- cbind(rnorm(20), rnorm(20))
  g <- knn_graph(coords, k = 4)
  expect_s3_class(g, "gsp_graph")
  expect_equal(g$n, 20)
})

test_that("graph_sbm builds stochastic block model graph", {
  set.seed(42)
  P <- matrix(c(0.5, 0.1, 0.1, 0.5), nrow = 2)
  g <- graph_sbm(block_sizes = c(10, 10), P = P)
  expect_s3_class(g, "gsp_graph")
  expect_equal(g$n, 20)
})

test_that("graph_random_regular builds regular graph", {
  set.seed(42)
  g <- graph_random_regular(n = 20, degree = 4, seed = 42)
  expect_s3_class(g, "gsp_graph")
  expect_equal(g$n, 20)
})

test_that("graph_barabasi_albert builds scale-free graph", {
  set.seed(42)
  g <- graph_barabasi_albert(n = 20, m = 2, seed = 42)
  expect_s3_class(g, "gsp_graph")
  expect_equal(g$n, 20)
})
