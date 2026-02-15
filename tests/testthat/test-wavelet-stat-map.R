test_that("wavelet_stat_map recovers mean signal", {
  set.seed(123)
  g <- graph_ring(15)
  X <- matrix(rnorm(15 * 4), nrow = 15)

  res <- wavelet_stat_map(g, X, stat = mean)

  expect_equal(length(res$map), 15)
  expect_lt(max(abs(res$map - rowMeans(X))), 1e-3)
})


test_that("sgwt_coeffs_tidy returns expected long rows", {
  set.seed(1)
  g <- graph_ring(10)
  X <- matrix(rnorm(10 * 3), nrow = 10)
  W <- sgwt(g, X, scales = c(2, 4))

  df <- sgwt_coeffs_tidy(W, long = TRUE)

  dims <- dim(W$coeffs)
  n_nodes <- dims[1]
  n_bands <- dims[2]
  n_signals <- if (length(dims) == 3) dims[3] else 1

  expect_equal(nrow(df), n_nodes * n_bands * n_signals)

  expected_bands <- if (W$include_lowpass) {
    c("lowpass", paste0("scale_", W$scales))
  } else {
    paste0("scale_", W$scales)
  }
  expect_setequal(unique(df$band), expected_bands)
})


test_that("band_energy_df percentages sum to 100", {
  set.seed(2)
  g <- graph_ring(12)
  x <- rnorm(12)
  W <- sgwt(g, x, scales = c(1, 3))

  df <- band_energy_df(W)

  expect_equal(nrow(df), dim(W$coeffs)[2])
  expect_true(all(df$pct >= 0))
  expect_lt(abs(sum(df$pct) - 100), 1e-6)
})
