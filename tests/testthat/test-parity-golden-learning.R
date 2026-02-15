# Golden-value parity tests: Tikhonov learning against PyGSP

test_that("Golden Tikhonov output matches PyGSP", {
  ref <- skip_if_no_golden("golden_learning.json")
  entry <- ref$tikhonov_ring20

  g <- graph_ring(20)
  x <- entry$signal
  tau <- entry$tau

  y_r <- tikhonov_smooth(g, x, tau = tau)

  expect_equal(as.numeric(y_r), entry$output, tolerance = tol_tikhonov)
})

test_that("Golden Tikhonov output is smoother than input", {
  ref <- skip_if_no_golden("golden_learning.json")
  entry <- ref$tikhonov_ring20

  g <- graph_ring(20)
  x <- entry$signal
  tau <- entry$tau

  y_r <- tikhonov_smooth(g, x, tau = tau)
  L <- as.matrix(graph_laplacian(g))

  qf_input <- as.numeric(t(x) %*% L %*% x)
  qf_output <- as.numeric(t(y_r) %*% L %*% y_r)

  expect_true(qf_output < qf_input)
})

test_that("Golden Tikhonov Laplacian quadratic form matches PyGSP", {
  ref <- skip_if_no_golden("golden_learning.json")
  entry <- ref$tikhonov_ring20

  g <- graph_ring(20)
  x <- entry$signal
  tau <- entry$tau

  y_r <- tikhonov_smooth(g, x, tau = tau)
  L <- as.matrix(graph_laplacian(g))

  qf_output_r <- as.numeric(t(y_r) %*% L %*% y_r)

  # Quadratic forms should be in the same ballpark
  expect_equal(qf_output_r, entry$laplacian_qf_output, tolerance = 0.01)
})

test_that("Golden regression_tikhonov wrapper matches tikhonov_smooth", {
  ref <- skip_if_no_golden("golden_learning.json")
  entry <- ref$tikhonov_ring20

  g <- graph_ring(20)
  x <- entry$signal
  tau <- entry$tau

  y_smooth <- tikhonov_smooth(g, x, tau = tau)
  y_regress <- regression_tikhonov(g, x, tau = tau)

  expect_equal(as.numeric(y_regress), as.numeric(y_smooth), tolerance = 1e-10)
})
