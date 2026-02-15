test_that("singleton graph builds and has zero laplacian", {
  A <- Matrix::Matrix(0, 1, 1, sparse = TRUE)
  g <- new_graph(A)
  expect_equal(g$n, 1L)
  expect_equal(as.numeric(g$laplacian), 0)
})

test_that("empty 3-node graph has zero degrees and laplacian", {
  A <- Matrix::Matrix(0, 3, 3, sparse = TRUE)
  g <- new_graph(A)
  expect_equal(as.numeric(g$degree), c(0, 0, 0))
  expect_true(all(g$laplacian == 0))
})

test_that("new_graph rejects nonsquare adjacency", {
  A <- Matrix::Matrix(0, 2, 3, sparse = TRUE)
  expect_error(new_graph(A), "square")
})

test_that("directed adjacency passes symmetry check", {
  A <- Matrix::sparseMatrix(i = c(1, 2), j = c(2, 3), x = 1, dims = c(3, 3))
  expect_error(new_graph(A, directed = FALSE), "symmetric")
  g <- new_graph(A, directed = TRUE)
  expect_true(g$directed)
})

test_that("graph_sparsify drops tiny weights", {
  g <- graph_ring(5)
  adj <- g$adjacency
  adj[1, 2] <- 1e-15
  adj[2, 1] <- 1e-15
  g$adjacency <- adj
  g2 <- graph_sparsify(g, tol = 1e-12)
  expect_true(length(g2$adjacency@x) < length(g$adjacency@x))
})

test_that("graph_rescale scales degrees and laplacian", {
  g <- graph_ring(4)
  g2 <- graph_rescale(g, 2)
  expect_equal(as.numeric(graph_degree(g2)), as.numeric(graph_degree(g)) * 2)
  expect_equal(Matrix::diag(graph_laplacian(g2)), Matrix::diag(graph_laplacian(g)) * 2)
})

test_that("graph_incidence oriented/unoriented sizes match expectations", {
  g <- graph_ring(5)
  inc1 <- graph_incidence(g, oriented = TRUE)
  inc2 <- graph_incidence(g, oriented = FALSE)
  expect_equal(nrow(inc1$B), g$n)
  expect_equal(nrow(inc2$B), length(g$adjacency@x))
})
