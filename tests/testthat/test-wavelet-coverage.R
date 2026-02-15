# Tests for wavelet_transform.R coverage
# Covers: sgwt, isgwt, sgwt_filter_bank, print/summary/plot methods,
# sgwt_coeffs_at_scale, sgwt_atom, sgwt_node_contribution,
# sgwt_reconstruct_at, sgwt_coeffs_tidy, band_energy_df,
# dsgwt_riesz, dsgwt_riesz_magnitude, directional_wavelet_knn

test_that("sgwt forward transform works with default params", {
  set.seed(42)
  g <- graph_ring(20)
  x <- rnorm(20)
  W <- sgwt(g, x, scales = c(1, 2, 4), K = 15)

  expect_s3_class(W, "sgwt_coeffs")
  expect_equal(dim(W$coeffs)[1], 20)  # n_nodes
  expect_equal(dim(W$coeffs)[2], 4)   # 3 scales + lowpass
  expect_equal(W$wavelet, "mexican_hat")
  expect_true(W$include_lowpass)
  expect_true(all(is.finite(W$coeffs)))
})

test_that("sgwt works with matrix input (multi-signal)", {
  set.seed(42)
  g <- graph_ring(15)
  X <- matrix(rnorm(15 * 3), nrow = 15)
  W <- sgwt(g, X, scales = c(1, 2), K = 15)

  expect_equal(dim(W$coeffs)[1], 15)
  expect_equal(dim(W$coeffs)[3], 3)   # 3 signals
  expect_false(W$single_signal)
})

test_that("sgwt rejects non-graph input", {
  expect_error(sgwt(list(), rnorm(5)), "must be a gsp_graph")
})

test_that("sgwt without lowpass", {
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), include_lowpass = FALSE, K = 15)
  expect_equal(dim(W$coeffs)[2], 2)
  expect_false(W$include_lowpass)
})

test_that("sgwt with different wavelet types", {
  g <- graph_ring(10)
  x <- rnorm(10)

  W_heat <- sgwt(g, x, scales = c(1, 2), wavelet = "heat", K = 15)
  expect_equal(W_heat$wavelet, "heat")

  W_custom <- sgwt(g, x, scales = c(1, 2),
                   wavelet = function(s) kernel_heat(s), K = 15)
  expect_equal(W_custom$wavelet, "custom")
})

test_that("sgwt_filter_bank builds correct banks", {
  g <- graph_ring(15)
  lmax <- lambda_max(g, normalized = g$normalized)

  fb <- sgwt_filter_bank(g, wavelet = "mexican_hat", scales = c(1, 2),
                         include_lowpass = TRUE, lmax = lmax)
  expect_equal(length(fb$kernels), 3)  # lowpass + 2 wavelet
  expect_true(fb$include_lowpass)
  expect_equal(fb$wavelet_name, "mexican_hat")

  fb2 <- sgwt_filter_bank(g, wavelet = "meyer", scales = c(1, 2),
                          include_lowpass = FALSE, lmax = lmax)
  expect_equal(length(fb2$kernels), 2)
  expect_false(fb2$include_lowpass)
})

test_that("sgwt_filter_bank errors on unknown wavelet", {
  g <- graph_ring(10)
  lmax <- lambda_max(g, normalized = g$normalized)
  expect_error(sgwt_filter_bank(g, wavelet = "bogus", scales = c(1),
                                lmax = lmax), "Unknown wavelet")
})

test_that("isgwt reconstructs signal (full, sum method)", {
  set.seed(42)
  g <- graph_ring(20)
  x <- rnorm(20)
  W <- sgwt(g, x, scales = c(1, 2, 4), K = 30)
  x_rec <- isgwt(W, method = "sum")

  expect_equal(length(x_rec), 20)
  expect_true(all(is.finite(x_rec)))
})

test_that("isgwt reconstructs signal (iterative method)", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  W <- sgwt(g, x, scales = c(1, 2), K = 20)
  x_rec <- isgwt(W, method = "iterative")

  expect_equal(length(x_rec), 15)
  expect_true(all(is.finite(x_rec)))
})

test_that("isgwt with band selection", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  W <- sgwt(g, x, scales = c(1, 2, 4), K = 15)

  rec_lp <- isgwt(W, bands = "lowpass")
  expect_equal(length(rec_lp), 15)

  rec_wav <- isgwt(W, bands = "wavelets")
  expect_equal(length(rec_wav), 15)

  rec_scale <- isgwt(W, bands = c(1, 2))
  expect_equal(length(rec_scale), 15)

  rec_name <- isgwt(W, bands = c("lowpass", "scale_1"))
  expect_equal(length(rec_name), 15)
})

test_that("isgwt with exclude_bands", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)

  rec <- isgwt(W, exclude_bands = "lowpass")
  expect_equal(length(rec), 15)

  rec2 <- isgwt(W, exclude_bands = c("scale_1"))
  expect_equal(length(rec2), 15)

  rec3 <- isgwt(W, exclude_bands = c(1))
  expect_equal(length(rec3), 15)
})

test_that("isgwt errors on invalid band names", {
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)

  expect_error(isgwt(W, bands = "nonexistent"), "not found")
  expect_error(isgwt(W, bands = c(999)), "not found")
})

test_that("isgwt rejects non-sgwt_coeffs input", {
  expect_error(isgwt(list()), "must be an sgwt_coeffs")
})

test_that("print.sgwt_coeffs works", {
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)
  out <- capture.output(print(W))
  expect_true(any(grepl("SGWT", out)))
  expect_true(any(grepl("Frame", out)))
})

test_that("summary.sgwt_coeffs works", {
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)
  out <- capture.output(res <- summary(W))
  expect_true(any(grepl("Energy", out)))
  expect_true(is.list(res))
})

test_that("plot.sgwt_coeffs energy type works", {
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)
  expect_no_error(plot(W, type = "energy"))
})

test_that("plot.sgwt_coeffs coeffs type works", {
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)
  expect_no_error(plot(W, type = "coeffs"))
})

test_that("plot.sgwt_coeffs heatmap type works", {
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)
  expect_no_error(plot(W, type = "heatmap"))
})

test_that("sgwt_coeffs_at_scale extracts correct bands", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  W <- sgwt(g, x, scales = c(1, 2, 4), K = 15)

  lp <- sgwt_coeffs_at_scale(W, "lowpass")
  expect_equal(length(lp), 15)

  s1 <- sgwt_coeffs_at_scale(W, 1)
  expect_equal(length(s1), 15)
})

test_that("sgwt_coeffs_at_scale errors on invalid scale", {
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)
  expect_error(sgwt_coeffs_at_scale(W, 999), "not found")
})

test_that("sgwt_atom computes localized atom", {
  set.seed(42)
  g <- graph_ring(20)
  atom <- sgwt_atom(g, node = 5, scale = 2, K = 15)
  expect_equal(length(atom), 20)
  expect_true(all(is.finite(atom)))
  expect_equal(attr(atom, "node"), 5)
  expect_equal(attr(atom, "scale"), 2)
})

test_that("sgwt_atom errors on invalid inputs", {
  g <- graph_ring(10)
  expect_error(sgwt_atom(g, node = 0, scale = 1), "must be between")
  expect_error(sgwt_atom(g, node = 11, scale = 1), "must be between")
  expect_error(sgwt_atom(list(), node = 1, scale = 1), "must be a gsp_graph")
})

test_that("sgwt_node_contribution computes energy", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)

  contrib <- sgwt_node_contribution(W, type = "energy")
  expect_equal(nrow(contrib), 15)
  expect_true(all(contrib >= 0))

  contrib_abs <- sgwt_node_contribution(W, type = "abs")
  expect_true(all(contrib_abs >= 0))

  contrib_raw <- sgwt_node_contribution(W, type = "raw")
  expect_equal(nrow(contrib_raw), 15)
})

test_that("sgwt_node_contribution normalizes correctly", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)

  contrib_norm <- sgwt_node_contribution(W, type = "energy", normalize = TRUE)
  row_sums <- rowSums(contrib_norm)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("sgwt_node_contribution rejects non-sgwt_coeffs", {
  expect_error(sgwt_node_contribution(list()), "must be an sgwt_coeffs")
})

test_that("sgwt_reconstruct_at works", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)

  val <- sgwt_reconstruct_at(W, nodes = 5)
  expect_equal(length(val), 1)

  vals <- sgwt_reconstruct_at(W, nodes = 1:5)
  expect_equal(length(vals), 5)

  all_vals <- sgwt_reconstruct_at(W, nodes = NULL)
  expect_equal(length(all_vals), 15)
})

test_that("sgwt_reconstruct_at errors on invalid nodes", {
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)
  expect_error(sgwt_reconstruct_at(W, nodes = 0), "must be between")
  expect_error(sgwt_reconstruct_at(W, nodes = 11), "must be between")
  expect_error(sgwt_reconstruct_at(list()), "must be an sgwt_coeffs")
})

test_that("sgwt_coeffs_tidy wide format works for single signal", {
  set.seed(42)
  g <- graph_ring(10)
  x <- rnorm(10)
  W <- sgwt(g, x, scales = c(1, 2), K = 15)

  df_wide <- sgwt_coeffs_tidy(W, long = FALSE)
  expect_equal(nrow(df_wide), 10)
  expect_true("node" %in% names(df_wide))
})

test_that("sgwt_coeffs_tidy errors on wide + multi-signal", {
  g <- graph_ring(10)
  X <- matrix(rnorm(10 * 3), nrow = 10)
  W <- sgwt(g, X, scales = c(1, 2), K = 15)
  expect_error(sgwt_coeffs_tidy(W, long = FALSE), "Wide format only")
})

test_that("sgwt_coeffs_tidy rejects non-sgwt_coeffs", {
  expect_error(sgwt_coeffs_tidy(list()), "must be an sgwt_coeffs")
})

test_that("band_energy_df rejects non-sgwt_coeffs", {
  expect_error(band_energy_df(list()), "must be an sgwt_coeffs")
})

test_that("dsgwt_riesz computes directional coefficients", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  dr <- dsgwt_riesz(g, x, scales = c(1, 2), K = 15)

  expect_s3_class(dr, "dsgwt_riesz")
  expect_equal(dim(dr$node_coeffs)[1], 15)
  expect_true(nrow(dr$edge_coeffs) > 0)
})

test_that("print.dsgwt_riesz works", {
  set.seed(42)
  g <- graph_ring(10)
  x <- rnorm(10)
  dr <- dsgwt_riesz(g, x, scales = c(1, 2), K = 15)
  out <- capture.output(print(dr))
  expect_true(any(grepl("Riesz", out)))
})

test_that("dsgwt_riesz_magnitude computes node magnitudes", {
  set.seed(42)
  g <- graph_ring(15)
  x <- rnorm(15)
  dr <- dsgwt_riesz(g, x, scales = c(1, 2), K = 15)
  mag <- dsgwt_riesz_magnitude(dr)

  expect_equal(dim(mag)[1], 15)  # n_nodes
  expect_true(all(mag >= 0))     # magnitudes are non-negative
})

test_that("dsgwt_riesz_magnitude rejects wrong class", {
  expect_error(dsgwt_riesz_magnitude(list()), "must be a dsgwt_riesz")
})

test_that("dsgwt_steer computes steered coefficients", {
  set.seed(42)
  coords <- cbind(cos(seq(0, 2*pi, length.out = 16)[1:15]),
                  sin(seq(0, 2*pi, length.out = 16)[1:15]))
  g <- graph_knn(coords, k = 4, normalized = FALSE)
  x <- rnorm(15)

  ds <- dsgwt_steer(g, x, n_directions = 2, scales = c(1, 2), K = 15)
  expect_s3_class(ds, "dsgwt_steer")
  expect_equal(ds$n_dir, 2)
  expect_equal(dim(ds$coeffs)[1], 15)
})

test_that("print.dsgwt_steer works", {
  set.seed(42)
  coords <- cbind(cos(seq(0, 2*pi, length.out = 16)[1:15]),
                  sin(seq(0, 2*pi, length.out = 16)[1:15]))
  g <- graph_knn(coords, k = 4, normalized = FALSE)
  x <- rnorm(15)
  ds <- dsgwt_steer(g, x, n_directions = 2, scales = c(1, 2), K = 15)
  out <- capture.output(print(ds))
  expect_true(any(grepl("steered", out)))
})

test_that("directional_wavelet_knn steer method works", {
  set.seed(42)
  coords <- cbind(rnorm(20), rnorm(20))
  signal <- rnorm(20)

  result <- directional_wavelet_knn(coords, signal,
                                     method = "steer", k = 4,
                                     n_directions = 2, scales = c(1, 2),
                                     K = 15)
  expect_s3_class(result, "dsgwt_steer")
})

test_that("directional_wavelet_knn riesz method works", {
  set.seed(42)
  coords <- cbind(rnorm(20), rnorm(20))
  signal <- rnorm(20)

  result <- directional_wavelet_knn(coords, signal,
                                     method = "riesz", k = 4,
                                     scales = c(1, 2), K = 15)
  expect_s3_class(result, "dsgwt_riesz")
})
