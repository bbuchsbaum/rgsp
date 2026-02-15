test_that("dsgwt_steer handles directed graphs with random-walk Laplacian", {
  set.seed(42)
  g <- graph_erdos_renyi(12, p = 0.3, directed = TRUE, normalized = FALSE)
  coords <- matrix(rnorm(g$n * 2), ncol = 2)
  g$coords <- coords
  x <- rnorm(g$n)

  # Should run without error and produce finite coefficients
  res <- dsgwt_steer(g, x,
                     n_directions = 2,
                     scales = c(1.5, 3),
                     wavelet = "heat",
                     include_lowpass = TRUE,
                     K = 10,
                     alpha_isotropic = 0.3,
                     p = 2)
  expect_true(all(is.finite(res$coeffs)))
  expect_equal(dim(res$coeffs)[1], g$n)
})
