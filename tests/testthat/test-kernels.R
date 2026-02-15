test_that("kernel_exponential decays with lambda", {
  k <- kernel_exponential(alpha = 2, beta = 1)
  lam <- c(0, 1, 2, 4)
  vals <- k(lam)
  expect_true(all(vals <= 1))
  expect_true(vals[1] == 1)
  expect_true(isTRUE(all(diff(vals) <= 0)))
})

test_that("kernel_mexican_hat peaks near sigma", {
  k <- kernel_mexican_hat(sigma = 1.5)
  lam <- seq(0, 5, length.out = 50)
  vals <- k(lam)
  expect_true(max(vals) > vals[1])
  expect_true(which.max(vals) >= 5)
})

test_that("kernel_gabor centered at mu", {
  k <- kernel_gabor(mu = 2, sigma = 0.5)
  lam <- c(0, 1, 2, 3, 4)
  vals <- k(lam)
  expect_equal(which.max(vals), 3) # index of lambda=2
})

test_that("kernel_rectangle is binary within band", {
  k <- kernel_rectangle(1, 3)
  lam <- seq(0, 4, by = 0.5)
  vals <- k(lam)
  expect_true(all(vals %in% c(0, 1)))
  expect_true(all(vals[lam >= 1 & lam <= 3] == 1))
})

test_that("kernel_half_cosine transitions smoothly", {
  k <- kernel_half_cosine(1, 2)
  vals <- k(c(0.5, 1, 1.5, 2, 2.5))
  expect_equal(vals[1], 1)
  expect_equal(vals[2], 1)
  expect_equal(vals[5], 0)
  expect_true(vals[3] < 1 && vals[3] > 0)
})

test_that("tight_frame_bank builds correct number of bands", {
  edges <- c(0, 1, 2, 4)
  bank <- tight_frame_bank(edges)
  expect_equal(length(bank), 3)
  lam <- c(0.5, 1.5, 3)
  vals <- sapply(bank, function(k) k(lam))
  expect_true(all(vals >= 0 & vals <= 1))
})
