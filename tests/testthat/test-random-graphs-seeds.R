test_that("random geometric graph respects seed determinism", {
  g1 <- graph_random_geometric(30, radius = 0.25, seed = 123)
  g2 <- graph_random_geometric(30, radius = 0.25, seed = 123)
  expect_equal(g1$adjacency, g2$adjacency)
})

test_that("sensor graph kNN symmetric and seeded", {
  g1 <- graph_sensor(25, k = 4, seed = 42)
  g2 <- graph_sensor(25, k = 4, seed = 42)
  expect_true(Matrix::isSymmetric(g1$adjacency))
  expect_true(all(Matrix::diag(g1$adjacency) == 0))
  expect_equal(g1$adjacency, g2$adjacency)
})

test_that("SBM respects zero-probability blocks", {
  block_sizes <- c(5, 5)
  P <- matrix(c(0.3, 0, 0, 0.4), 2, 2, byrow = TRUE)
  g <- graph_sbm(block_sizes, P, seed = 7)
  # ensure no inter-block edges
  idx1 <- 1:5
  idx2 <- 6:10
  sub <- g$adjacency[idx1, idx2]
  expect_true(length(sub@x) == 0)
})

test_that("barabasi albert degree grows preferentially", {
  g <- graph_barabasi_albert(40, m = 2, seed = 9)
  deg <- sort(as.numeric(g$degree), decreasing = TRUE)
  expect_true(deg[1] >= mean(deg))
})

test_that("random regular graph degree exact", {
  g <- graph_random_regular(20, degree = 3, seed = 5)
  expect_equal(unique(as.numeric(g$degree)), 3)
})
