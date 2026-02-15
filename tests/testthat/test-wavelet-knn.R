test_that("wavelet_heat_knn runs and returns per-scale list", {
  set.seed(1)
  coords <- matrix(runif(20), ncol = 2)
  signal <- rnorm(nrow(coords))
  scales <- c(1, 2)

  res <- wavelet_heat_knn(coords, signal, scales = scales, k = 3, K = 5)
  expect_length(res, length(scales))
  expect_named(res, paste0("scale_", scales))
})

test_that("directional_wavelet_knn (steer) produces steer object", {
  set.seed(2)
  coords <- matrix(runif(20), ncol = 2)
  signal <- rnorm(nrow(coords))
  res <- directional_wavelet_knn(coords, signal,
                                 method = "steer",
                                 k = 3,
                                 scales = c(1, 2),
                                 wavelet = "mexican_hat",
                                 n_directions = 2,
                                 K = 5)
  expect_s3_class(res, "dsgwt_steer")
  expect_equal(dim(res$coeffs)[3], 2) # directions
})
