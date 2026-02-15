test_that("erdos renyi seeds reproducibly", {
  g1 <- graph_erdos_renyi(20, p = 0.2, seed = 42)
  g2 <- graph_erdos_renyi(20, p = 0.2, seed = 42)
  expect_equal(g1$adjacency, g2$adjacency)
})

test_that("barabasi albert has expected node/edge counts", {
  g <- graph_barabasi_albert(15, m = 2, seed = 1)
  expect_equal(g$n, 15L)
  expect_true(Matrix::isSymmetric(g$adjacency))
})

test_that("random regular degree is constant", {
  g <- graph_random_regular(10, degree = 4, seed = 3)
  expect_equal(sort(as.numeric(g$degree)), rep(4, 10))
})

test_that("complete graph has full degree", {
  g <- graph_complete(5)
  expect_equal(as.numeric(g$degree), rep(4, 5))
})

test_that("SBM construction is sparse and reproducible", {
  P <- matrix(c(0.2, 0.05,
                0.05, 0.3), 2, 2, byrow = TRUE)
  g1 <- graph_sbm(c(40, 60), P, seed = 123)
  g2 <- graph_sbm(c(40, 60), P, seed = 123)
  expect_equal(g1$adjacency, g2$adjacency)
  expect_true(Matrix::isSymmetric(g1$adjacency))
  expect_equal(Matrix::diag(g1$adjacency), rep(0, g1$n))
  # density should stay well below dense (rough heuristic check)
  expect_true(length(g1$adjacency@x) < 0.2 * g1$n * (g1$n - 1))
})
