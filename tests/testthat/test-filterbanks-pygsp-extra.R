test_that("Abspline/Itersine/SimpleTight filter banks have expected shapes", {
  g <- graph_ring(32)
  lmax <- lambda_max(g)

  abs_k <- abspline_filter_bank(lmax, Nf = 4)
  it_k <- itersine_filter_bank(lmax, Nf = 4, overlap = 2)
  st_k <- simpletight_filter_bank(lmax, Nf = 4)

  expect_equal(length(abs_k), 4)
  expect_equal(length(it_k), 4)
  expect_equal(length(st_k), 4)

  x <- rnorm(g$n)
  expect_equal(length(filter_bank_apply(g, x, abs_k, K = 30, lmax = lmax)), 4)
  expect_equal(length(filter_bank_apply(g, x, it_k, K = 30, lmax = lmax)), 4)
  expect_equal(length(filter_bank_apply(g, x, st_k, K = 30, lmax = lmax)), 4)
})

test_that("Itersine and SimpleTight are (approximately) tight partitions of unity", {
  g <- graph_ring(64)
  lmax <- lambda_max(g)
  lambdas <- seq(0, lmax, length.out = 401)

  it_k <- itersine_filter_bank(lmax, Nf = 6, overlap = 2)
  st_k <- simpletight_filter_bank(lmax, Nf = 6)

  sumsq <- function(kernels) {
    resp <- vapply(kernels, function(k) k(lambdas), numeric(length(lambdas)))
    rowSums(resp^2)
  }

  it_s <- sumsq(it_k)
  st_s <- sumsq(st_k)

  expect_lt(max(abs(it_s - 1)), 1e-6)
  expect_lt(max(abs(st_s - 1)), 1e-6)
})

