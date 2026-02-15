# Minimal directional SGWT demo on a 2D grid (domain-agnostic)
#
# Builds a small grid graph, injects an oriented stripe signal, runs
# dsgwt_steer, and reports which direction has the strongest response.

suppressPackageStartupMessages({
  library(rgsp)
})

run_directional_demo <- function(nx = 32, ny = 32, stripe_angle = "horizontal",
                                 D = 4, n_scales = 2) {
  g <- graph_grid2d(nx, ny, normalized = FALSE)
  g <- graph_set_coords(g, kind = g$coords)

  signal <- rep(0, g$n)
  coords <- g$coords
  if (stripe_angle == "horizontal") {
    signal[coords[, 2] %in% c(15, 16, 17)] <- 1
  } else if (stripe_angle == "vertical") {
    signal[coords[, 1] %in% c(15, 16, 17)] <- 1
  } else {
    # diagonal
    mask <- abs(coords[, 1] - coords[, 2]) <= 1
    signal[mask] <- 1
  }

  scales <- auto_wavelet_scales(g, n_scales = n_scales, lmax = lambda_max(g))
  res <- dsgwt_steer(g, signal,
                     n_directions = D,
                     scales = scales,
                     wavelet = "mexican_hat",
                     include_lowpass = TRUE,
                     K = 20,
                     alpha_isotropic = 0.2,
                     p = 2)

  # Simple energy summary per direction (aggregate over nodes/bands)
  coeffs <- res$coeffs
  # drop lowpass (band 1)
  wavelet_bands <- 2:dim(coeffs)[2]
  energy_dir <- apply(coeffs[, wavelet_bands, , drop = FALSE]^2, 3, sum)
  best_dir <- which.max(energy_dir)
  cat("Stripe:", stripe_angle, " -> strongest direction index:", best_dir, "\n")
  invisible(list(energy = energy_dir, best_dir = best_dir, orientations = res$orientations))
}

if (interactive()) {
  run_directional_demo(stripe_angle = "horizontal", D = 4)
  run_directional_demo(stripe_angle = "vertical", D = 4)
}
