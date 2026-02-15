test_that("heat_propagate returns list for multiple times", {
  g <- graph_ring(5)
  x <- rnorm(5)
  res <- heat_propagate(g, x, t = c(0.1, 0.2), K = 20)
  expect_type(res, "list")
  expect_equal(length(res), 2)
  expect_equal(length(res[[1]]), 5)
})

test_that("heat_propagate single time returns vector", {
  g <- graph_ring(4)
  x <- rnorm(4)
  res <- heat_propagate(g, x, t = 0.1, K = 20)
  expect_type(res, "double")
  expect_equal(length(res), 4)
})

test_that("heat_propagate matches filter_signal heat kernel", {
  g <- graph_ring(6)
  x <- rnorm(6)
  t <- 0.2
  y1 <- heat_propagate(g, x, t = t, K = 40)
  y2 <- filter_signal(g, x, kernel_heat(t), K = 40)
  expect_equal(y1, y2, tolerance = 1e-8)
})
