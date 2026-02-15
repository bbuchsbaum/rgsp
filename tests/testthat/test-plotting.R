# Tests for plotting functions

# ==============================================================================
# Layout Tests
# ==============================================================================

test_that("graph_set_coords with circle layout", {
  g <- graph_ring(10)
  g <- graph_set_coords(g, kind = "circle")

  expect_true(!is.null(g$coords))

  expect_equal(nrow(g$coords), 10)
  expect_equal(ncol(g$coords), 2)

  # Circle layout should have points on unit circle
  radii <- sqrt(rowSums(g$coords^2))
  expect_equal(radii, rep(1, 10), tolerance = 1e-10)
})

test_that("graph_set_coords with spring layout", {
  g <- graph_grid2d(3, 3)
  g <- graph_set_coords(g, kind = "spring", seed = 42)

  expect_true(!is.null(g$coords))
  expect_equal(nrow(g$coords), 9)
  expect_equal(ncol(g$coords), 2)
})

test_that("graph_set_coords with spring layout is reproducible with seed", {
  g <- graph_grid2d(3, 3)

  g1 <- graph_set_coords(g, kind = "spring", seed = 123)
  g2 <- graph_set_coords(g, kind = "spring", seed = 123)

  expect_equal(g1$coords, g2$coords)
})

test_that("graph_set_coords with spectral layout", {
  g <- graph_ring(8)
  g <- graph_set_coords(g, kind = "spectral")

  expect_true(!is.null(g$coords))
  expect_equal(nrow(g$coords), 8)
  expect_equal(ncol(g$coords), 2)
})

test_that("graph_set_coords with grid layout", {
  g <- graph_grid2d(4, 5)
  g <- graph_set_coords(g, kind = "grid")

  expect_true(!is.null(g$coords))
  expect_equal(nrow(g$coords), 20)
  expect_equal(ncol(g$coords), 2)
})

test_that("graph_set_coords with random layout", {
  g <- graph_ring(10)
  g <- graph_set_coords(g, kind = "random", seed = 42)

  expect_true(!is.null(g$coords))
  expect_equal(nrow(g$coords), 10)
  expect_equal(ncol(g$coords), 2)

  # Random should be bounded in [0, 1]
  expect_true(all(g$coords >= 0 & g$coords <= 1))
})

test_that("graph_set_coords with line layout", {
  g <- graph_path(10)
  g <- graph_set_coords(g, kind = "line")

  expect_true(!is.null(g$coords))
  expect_equal(nrow(g$coords), 10)
  expect_equal(ncol(g$coords), 2)

  # Line layout should have y = 0 and x equally spaced
  expect_equal(g$coords[, 2], rep(0, 10))
})

test_that("graph_set_coords with 3D spring layout", {
  g <- graph_grid2d(3, 3)
  g <- graph_set_coords(g, kind = "spring", dim = 3, seed = 42)

  expect_true(!is.null(g$coords))
  expect_equal(nrow(g$coords), 9)
  expect_equal(ncol(g$coords), 3)
})

test_that("graph_set_coords with 3D random layout", {
  g <- graph_ring(10)
  g <- graph_set_coords(g, kind = "random", dim = 3, seed = 42)

  expect_true(!is.null(g$coords))
  expect_equal(nrow(g$coords), 10)
  expect_equal(ncol(g$coords), 3)
})

test_that("graph_set_coords errors on invalid kind", {
  g <- graph_ring(10)
  expect_error(graph_set_coords(g, kind = "invalid"))
})

# ==============================================================================
# Plot Function Tests (require ggplot2)
# ==============================================================================

skip_if_no_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    testthat::skip("ggplot2 not installed")
  }
}

test_that("plot_graph returns ggplot object", {
  skip_if_no_ggplot2()

  g <- graph_ring(10)
  g <- graph_set_coords(g, kind = "circle")

  p <- plot_graph(g)
  expect_s3_class(p, "ggplot")
})

test_that("plot_graph with signal", {
  skip_if_no_ggplot2()

  g <- graph_ring(10)
  g <- graph_set_coords(g, kind = "circle")
  signal <- sin(seq(0, 2*pi, length.out = 10))

  p <- plot_graph(g, signal = signal)
  expect_s3_class(p, "ggplot")
})

test_that("plot_graph auto-generates layout if missing", {
  skip_if_no_ggplot2()

  g <- graph_ring(10)
  # No coords set
  expect_null(g$coords)

  p <- plot_graph(g)
  expect_s3_class(p, "ggplot")
})

test_that("plot_signal returns ggplot object", {
  skip_if_no_ggplot2()

  g <- graph_ring(10)
  g <- graph_set_coords(g, kind = "circle")
  signal <- rnorm(10)

  p <- plot_signal(g, signal)
  expect_s3_class(p, "ggplot")
})

test_that("plot_signal with custom vertex size", {
  skip_if_no_ggplot2()

  g <- graph_ring(10)
  g <- graph_set_coords(g, kind = "circle")
  signal <- abs(rnorm(10))

  p <- plot_signal(g, signal, vertex_size = 5)
  expect_s3_class(p, "ggplot")
})

test_that("plot_filter_response with single kernel", {
  skip_if_no_ggplot2()

  p <- plot_filter_response(kernel_heat(1), lmax = 4)
  expect_s3_class(p, "ggplot")
})

test_that("plot_filter_response with multiple kernels", {
  skip_if_no_ggplot2()

  kernels <- list(
    Heat_1 = kernel_heat(1),
    Heat_5 = kernel_heat(5),
    Mexican = kernel_mexican_hat(1)
  )

  p <- plot_filter_response(kernels, lmax = 4)
  expect_s3_class(p, "ggplot")
})

test_that("plot_wavelet_coeffs returns ggplot object", {
  skip_if_no_ggplot2()

  g <- graph_ring(12)
  signal <- rnorm(12)
  W <- sgwt(g, signal, scales = 3)

  p <- plot_wavelet_coeffs(W)
  expect_s3_class(p, "ggplot")
})

test_that("plot_wavelet_coeffs with specific signal_idx", {
  skip_if_no_ggplot2()

  g <- graph_ring(12)
  signal <- matrix(rnorm(24), nrow = 12, ncol = 2)
  W <- sgwt(g, signal, scales = 3)

  p <- plot_wavelet_coeffs(W, signal_idx = 2)
  expect_s3_class(p, "ggplot")
})

test_that("plot_band_energy returns ggplot object", {
  skip_if_no_ggplot2()

  g <- graph_ring(12)
  signal <- rnorm(12)
  W <- sgwt(g, signal, scales = 3)

  p <- plot_band_energy(W)
  expect_s3_class(p, "ggplot")
})

test_that("plot_spectrogram returns ggplot object", {
  skip_if_no_ggplot2()

  g <- graph_ring(12)
  signal <- rnorm(12)

  p <- plot_spectrogram(g, signal)
  expect_s3_class(p, "ggplot")
})

test_that("plot_spectrogram with custom scales", {
  skip_if_no_ggplot2()

  g <- graph_ring(12)
  signal <- rnorm(12)

  p <- plot_spectrogram(g, signal, scales = c(0.5, 1, 2, 4))
  expect_s3_class(p, "ggplot")
})

# ==============================================================================
# Edge cases
# ==============================================================================

test_that("layout works on small graphs", {
  g <- graph_path(2)
  g <- graph_set_coords(g, kind = "spring", seed = 1)
  expect_equal(nrow(g$coords), 2)
})

test_that("spectral layout handles disconnected-like graphs gracefully", {
  # Very sparse graph
  g <- graph_path(5)
  g <- graph_set_coords(g, kind = "spectral")
  expect_equal(nrow(g$coords), 5)
})

# ==============================================================================
# Style Preset Tests
# ==============================================================================

test_that("theme_gsp returns valid themes", {
  skip_if_no_ggplot2()

  for (preset in c("default", "minimal", "dark", "publication")) {
    theme <- theme_gsp(preset)
    expect_s3_class(theme, "theme")
  }
})

test_that("theme_gsp errors on invalid preset", {
  skip_if_no_ggplot2()
  expect_error(theme_gsp("invalid"))
})

test_that("scale_signal returns valid scales", {
  skip_if_no_ggplot2()

  for (palette in c("viridis", "spectral", "diverging", "sequential", "heat")) {
    scale <- scale_signal(palette)
    expect_s3_class(scale, "Scale")
  }
})

test_that("scale_signal errors on invalid palette", {
  skip_if_no_ggplot2()
  expect_error(scale_signal("invalid"))
})

test_that("scale_fill_gsp returns valid scales", {
  skip_if_no_ggplot2()

  for (palette in c("viridis", "spectral", "magma", "inferno", "plasma")) {
    scale <- scale_fill_gsp(palette)
    expect_s3_class(scale, "Scale")
  }
})

test_that("scale_fill_gsp errors on invalid palette", {
  skip_if_no_ggplot2()
  expect_error(scale_fill_gsp("invalid"))
})

test_that("legend_position returns theme element", {
  skip_if_no_ggplot2()

  for (pos in c("right", "left", "top", "bottom", "none")) {
    lp <- legend_position(pos)
    expect_s3_class(lp, "theme")
  }
})

test_that("graph_title returns labs object", {
  skip_if_no_ggplot2()

  labs_obj <- graph_title("Main Title", "Subtitle", "Caption")
  # ggplot2 >= 3.5.0 returns ggplot2::labels class
  expect_true(inherits(labs_obj, "labels") || inherits(labs_obj, "ggplot2::labels"))
})

test_that("theme_gsp can be added to plot", {
  skip_if_no_ggplot2()

  g <- graph_ring(10)
  g <- graph_set_coords(g, "circle")
  p <- plot_graph(g) + theme_gsp("minimal")
  expect_s3_class(p, "ggplot")
})

test_that("scale_signal can be added to plot", {
  skip_if_no_ggplot2()

  g <- graph_ring(10)
  g <- graph_set_coords(g, "circle")
  signal <- sin(seq(0, 2*pi, length.out = 10))

  # Note: scale_signal sets colour scale but plot_signal uses colour aesthetic
  p <- plot_signal(g, signal) + scale_signal("diverging")
  expect_s3_class(p, "ggplot")
})
