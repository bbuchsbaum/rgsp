test_that("prox_tv1d shrinks toward piecewise constant", {
  set.seed(1)
  y <- c(1, 2, 3, 10, 11)
  x <- prox_tv1d(y, lambda = 3)
  tv_y <- sum(abs(diff(y)))
  tv_x <- sum(abs(diff(x)))
  expect_true(tv_x < tv_y)
})

test_that("tv_denoise reduces variation on ring", {
  set.seed(2)
  g <- graph_ring(10)
  y <- rnorm(10)
  x <- tv_denoise(g, y, lambda = 2, maxit = 120)
  expect_true(var(diff(x)) <= var(diff(y)) + 1e-6)
})

test_that("prox_box clamps values", {
  y <- c(-2, 0, 5)
  x <- prox_box(y, lower = -1, upper = 1)
  expect_equal(x, c(-1, 0, 1))
})

test_that("simplex projection sums to 1 and nonnegative", {
  y <- c(-1, 0.2, 3)
  x <- proj_simplex(y)
  expect_true(all(x >= -1e-12))
  expect_equal(sum(x), 1, tolerance = 1e-8)
})

test_that("l2 ball projection enforces radius", {
  y <- c(3, 4) # norm 5
  x <- proj_l2_ball(y, r = 2)
  expect_equal(round(sqrt(sum(x^2)), 6), 2)
})

test_that("tv_inpaint respects observed entries", {
  set.seed(3)
  g <- graph_ring(8)
  y <- rnorm(8)
  y_missing <- y
  y_missing[c(3, 5)] <- NA
  x <- tv_inpaint(g, y_missing, lambda = 1, maxit = 50)
  expect_equal(x[c(1, 2, 4, 6, 7, 8)], y[c(1, 2, 4, 6, 7, 8)], tolerance = 1e-8)
})
