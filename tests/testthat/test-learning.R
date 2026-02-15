test_that("regression_tikhonov returns finite values", {
  g <- graph_ring(10)
  y <- rnorm(10)
  x <- regression_tikhonov(g, y, tau = 0.5)
  expect_true(all(is.finite(x)))
})

test_that("classification_tikhonov produces one column per class", {
  g <- graph_ring(8)
  labels <- factor(c(1, 1, 2, 2, 3, 3, 1, 2))
  scores <- classification_tikhonov(g, labels, tau = 0.3)
  expect_equal(ncol(scores), 3)
  expect_true(all(colnames(scores) == c("class_1", "class_2", "class_3")))
  # Higher tau should smooth more (reduce variance across neighbors)
  scores2 <- classification_tikhonov(g, labels, tau = 1.0)
  var1 <- sum((scores[, 1] - mean(scores[, 1]))^2)
  var2 <- sum((scores2[, 1] - mean(scores2[, 1]))^2)
  expect_true(var2 < var1 + 1e-8)
})
