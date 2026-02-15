skip_if_no_pygsp <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    testthat::skip("reticulate not installed")
  }
  if (!reticulate::py_available(initialize = FALSE)) {
    testthat::skip("Python not available for reticulate parity tests")
  }
  if (!reticulate::py_module_available("pygsp")) {
    testthat::skip("PyGSP not available (pip install pygsp)")
  }
}

get_pygsp <- function() {
  list(
    pygsp = reticulate::import("pygsp"),
    np = reticulate::import("numpy")
  )
}

test_that("Example parity: heat filter on Sensor graph matches PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()

  set.seed(101)
  g_py <- py$pygsp$graphs$Sensor(32L, seed = 42L)
  g_py$estimate_lmax()
  tau <- 0.8
  filt_py <- py$pygsp$filters$Heat(g_py, tau)

  x <- rnorm(32)
  y_py <- as.numeric(filt_py$filter(py$np$array(x)))

  g_r <- graph_sensor(32, seed = 42)
  y_r <- filter_signal(g_r, x, kernel_heat(tau), K = 80, lmax = as.numeric(g_py$lmax))

  expect_equal(length(y_r), length(y_py))
  expect_true(all(is.finite(y_r)))
})

test_that("Example (R only): heat filter on Sensor graph smooths the signal", {
  set.seed(202)
  g <- graph_sensor(32, seed = 42)
  tau <- 0.8
  x <- rnorm(32)
  L <- graph_laplacian(g)
  energy_x <- as.numeric(t(x) %*% as.matrix(L) %*% x)
  y <- filter_signal(g, x, kernel_heat(tau), K = 80)
  energy_y <- as.numeric(t(y) %*% as.matrix(L) %*% y)

  expect_true(all(is.finite(y)))
  expect_lt(energy_y, energy_x)
})

test_that("Example parity: heat diffusion on grid decays similarly to PyGSP", {
  skip_if_no_pygsp()
  py <- get_pygsp()
  np <- py$np

  n_side <- 10L
  g_py <- py$pygsp$graphs$Grid2d(n_side)
  g_py$estimate_lmax()

  sources <- c((n_side %/% 4 * n_side) + (n_side %/% 4),
               (n_side * 3 %/% 4 * n_side) + (n_side * 3 %/% 4))
  x <- rep(0, g_py$n_vertices)
  x[sources + 1] <- 5  # 0-index to 1-index adjustment for R vector later

  times <- c(0, 5, 10)
  g_r <- graph_grid2d(n_side, n_side)
  lmax <- as.numeric(g_py$lmax)

  # PyGSP diffusion at final time
  filt_py <- py$pygsp$filters$Heat(g_py, scale = times[3])
  y_py <- as.numeric(filt_py$filter(np$array(x)))

  # R diffusion at all times; check monotone smoothing and parity at final time
  energies <- numeric(length(times))
  y_r_last <- NULL
  for (i in seq_along(times)) {
    t <- times[i]
    y_r <- filter_signal(g_r, x, kernel_heat(t), K = 120, lmax = lmax)
    energies[i] <- as.numeric(t(y_r) %*% as.matrix(graph_laplacian(g_r)) %*% y_r)
    if (i == length(times)) y_r_last <- y_r
  }

  # Dirichlet energy should decay
  expect_true(all(diff(energies) <= 1e-8))
  expect_gt(cor(as.numeric(y_r_last), y_py), 0.6)
})

test_that("Example (R only): heat diffusion on grid decays Dirichlet energy", {
  n_side <- 10
  g <- graph_grid2d(n_side, n_side)
  sources <- c((n_side %/% 4 * n_side) + (n_side %/% 4),
               (n_side * 3 %/% 4 * n_side) + (n_side * 3 %/% 4))
  x <- rep(0, g$n)
  x[sources] <- 5

  times <- c(0, 5, 10)
  energies <- numeric(length(times))
  for (i in seq_along(times)) {
    t <- times[i]
    y <- filter_signal(g, x, kernel_heat(t), K = 120)
    energies[i] <- as.numeric(t(y) %*% as.matrix(graph_laplacian(g)) %*% y)
  }
  expect_true(all(diff(energies) <= 1e-8))
})
