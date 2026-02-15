test_that("tikhonov_smooth returns finite values and reduces energy", {
  set.seed(3)
  g <- graph_ring(6)
  b <- rnorm(6)
  x <- tikhonov_smooth(g, b, tau = 0.5)
  expect_true(all(is.finite(x)))
  expect_true(sum(x^2) <= sum(b^2) + 1e-8)
})

test_that("prox_l1 shrinks to zero when tau large", {
  x <- matrix(c(-1, 0.5, 2), ncol = 1)
  y <- prox_l1(x, tau = 10)
  expect_equal(y, matrix(0, nrow = 3, ncol = 1))
})

test_that("proj_simplex projects negative input", {
  y <- c(-1, -2, -3)
  x <- proj_simplex(y)
  expect_true(all(x >= -1e-12))
  expect_equal(sum(x), 1, tolerance = 1e-8)
})

test_that("proj_l2_ball returns same vector if within radius", {
  y <- c(0.1, 0.2)
  x <- proj_l2_ball(y, r = 1)
  expect_equal(x, y)
})

test_that("tv_denoise reduces variation more with higher lambda", {
  set.seed(1)
  g <- graph_ring(20)
  y <- rnorm(20)
  x1 <- tv_denoise(g, y, lambda = 0.5, maxit = 80)
  x2 <- tv_denoise(g, y, lambda = 2.0, maxit = 80)
  expect_true(sum(abs(diff(x2))) <= sum(abs(diff(x1))) + 1e-6)
})

test_that("tv_inpaint leaves observed values intact", {
  set.seed(2)
  g <- graph_ring(12)
  y <- rnorm(12)
  mask <- c(2, 5, 9)
  y_miss <- y
  y_miss[mask] <- NA
  x <- tv_inpaint(g, y_miss, lambda = 1, maxit = 60)
  expect_equal(x[-mask], y[-mask], tolerance = 1e-6)
})
