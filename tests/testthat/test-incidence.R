test_that("incidence has correct dimensions and weights", {
  g <- graph_ring(5)
  inc <- graph_incidence(g)
  B <- inc$B
  expect_equal(dim(B), c(5, 5))
  # each row should sum to zero
  expect_equal(Matrix::rowSums(B), rep(0, 5))
  # degree from incidence matches adjacency degree
  deg_inc <- as.numeric(Matrix::colSums(abs(B)))
  expect_equal(deg_inc, as.numeric(g$degree))
})

test_that("incidence uses upper-tri orientation", {
  g <- graph_ring(4)
  inc <- graph_incidence(g)
  expect_true(all(inc$rows < inc$cols))
})
