test_that("graph_knn produces symmetric adjacency (union)", {
  coords <- matrix(c(0, 0,
                     1, 0,
                     2, 0,
                     0, 2,
                     2, 2), ncol = 2, byrow = TRUE)
  g <- graph_knn(coords, k = 2, sym = "union", weight = "binary", seed = 123)
  expect_true(Matrix::isSymmetric(g$adjacency))
  expect_equal(g$n, nrow(coords))
})

test_that("graph_knn mutual is subset of union", {
  coords <- matrix(runif(40), ncol = 2)
  g_union <- graph_knn(coords, k = 4, sym = "union", weight = "distance", seed = 1)
  g_mut   <- graph_knn(coords, k = 4, sym = "mutual", weight = "distance", seed = 1)
  # mutual adjacency should never have edges not present in union
  diff <- Matrix::drop0(g_mut$adjacency - pmin(g_mut$adjacency, g_union$adjacency))
  expect_equal(Matrix::nnzero(diff), 0)
  expect_true(Matrix::isSymmetric(g_mut$adjacency))
})

test_that("graph_knn heat weight auto sigma is finite", {
  coords <- matrix(runif(30), ncol = 3)
  g <- graph_knn(coords, k = 3, weight = "heat")
  expect_true(all(is.finite(g$adjacency@x)))
  expect_true(Matrix::isSymmetric(g$adjacency))
})
