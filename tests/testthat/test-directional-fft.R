test_that("directional energy aligns with FFT orientation on stripes", {
  skip_if_not_installed("stats")
  # tiny grid to keep FFT simple
  n <- 32
  g <- graph_grid2d(n, n, normalized = FALSE)

  make_stripe <- function(orientation = "horizontal") {
    signal <- rep(0, g$n)
    coords <- g$coords
    if (orientation == "horizontal") {
      signal[coords[, 2] %in% c(15, 16, 17)] <- 1
    } else if (orientation == "vertical") {
      signal[coords[, 1] %in% c(15, 16, 17)] <- 1
    }
    signal
  }

  stripe <- make_stripe("horizontal")
  scales <- auto_wavelet_scales(g, n_scales = 2, lmax = lambda_max(g))
  res <- dsgwt_steer(g, stripe, n_directions = 4, scales = scales,
                     wavelet = "mexican_hat", include_lowpass = TRUE, K = 15,
                     alpha_isotropic = 0.2, p = 2)
  # drop lowpass
  wavelet_bands <- 2:dim(res$coeffs)[2]
  coeffs <- res$coeffs
  if (length(dim(coeffs)) == 3) {
    # single signal stored as 3D: node x band x direction
    coeffs <- array(coeffs, dim = c(dim(coeffs), 1))
  }
  energy_dir <- apply(coeffs[, wavelet_bands, , , drop = FALSE]^2, 3, sum)
  best_dir <- which.max(energy_dir)

  # crude FFT orientation energy (horizontal stripe -> high vertical freq)
  img <- matrix(stripe, n, n, byrow = TRUE)
  F <- abs(stats::fft(img))
  # ignore DC
  kx <- colSums(F)[-1]
  ky <- rowSums(F)[-1]
  fft_vert_dom <- sum(ky) > sum(kx)

  expect_true(fft_vert_dom)
  expect_true(best_dir %in% c(2, 4))  # assuming directions cover axes
})
