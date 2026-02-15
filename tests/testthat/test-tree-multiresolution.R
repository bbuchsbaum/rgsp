is_tree_graph <- function(g) {
  e <- length(rgsp:::.edge_list_undirected(g$adjacency)$w)
  isTRUE(e == (g$n - 1L)) && isTRUE(rgsp:::.graph_is_connected(g$adjacency))
}

test_that("tree_multiresolution coarsens a path graph", {
  g <- graph_path(9)
  expect_true(is_tree_graph(g))

  res <- tree_multiresolution(g, levels = 2, root = 1, reduction_method = "unweighted")
  expect_equal(length(res$graphs), 3)
  expect_equal(length(res$keeps), 2)
  expect_equal(length(res$roots), 3)

  expect_equal(res$graphs[[2]]$n, 5L)
  expect_equal(res$graphs[[3]]$n, 3L)

  expect_true(is_tree_graph(res$graphs[[2]]))
  expect_true(is_tree_graph(res$graphs[[3]]))
})

test_that("tree_multiresolution weight combination methods behave", {
  g <- graph_path(7)

  r_res <- tree_multiresolution(g, levels = 1, root = 1, reduction_method = "resistance_distance")
  s_res <- tree_multiresolution(g, levels = 1, root = 1, reduction_method = "sum")
  u_res <- tree_multiresolution(g, levels = 1, root = 1, reduction_method = "unweighted")

  w_r <- Matrix::summary(Matrix::triu(r_res$graphs[[2]]$adjacency, 1))$x
  w_s <- Matrix::summary(Matrix::triu(s_res$graphs[[2]]$adjacency, 1))$x
  w_u <- Matrix::summary(Matrix::triu(u_res$graphs[[2]]$adjacency, 1))$x

  expect_true(all(abs(w_r - 0.5) < 1e-12))
  expect_true(all(abs(w_s - 2.0) < 1e-12))
  expect_true(all(abs(w_u - 1.0) < 1e-12))
})

