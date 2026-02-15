# Golden-value parity tests: Lanczos filter parity against PyGSP

test_that("Golden Lanczos filter matches exact filtered output", {
  ref <- skip_if_no_golden("golden_lanczos.json")

  g <- graph_ring(12)
  x <- ref$signal
  lmax_py <- ref$lmax
  scale <- ref$scale
  t_rgsp <- pygsp_heat_to_rgsp(scale, lmax_py)

  # Lanczos filter with K=30 steps
  y_lanczos <- filter_signal_lanczos(g, x, kernel_heat(t_rgsp), K = 30)

  # Should match exact PyGSP output closely
  expect_equal(as.numeric(y_lanczos), ref$filtered_exact, tolerance = tol_filter_chebyshev)
})

test_that("Golden Lanczos filter is closer to exact than Chebyshev at same order", {
  ref <- skip_if_no_golden("golden_lanczos.json")

  g <- graph_ring(12)
  x <- ref$signal
  lmax_py <- ref$lmax
  scale <- ref$scale
  t_rgsp <- pygsp_heat_to_rgsp(scale, lmax_py)

  y_lanczos <- filter_signal_lanczos(g, x, kernel_heat(t_rgsp), K = 12)
  y_cheby <- filter_signal(g, x, kernel_heat(t_rgsp), K = 12, lmax = lmax_py)

  exact <- ref$filtered_exact

  err_lanczos <- max(abs(as.numeric(y_lanczos) - exact))
  err_cheby <- max(abs(as.numeric(y_cheby) - exact))

  # Lanczos at full rank (K=n=12) should achieve machine precision
  expect_true(err_lanczos < 1e-8,
              label = paste("Lanczos error:", signif(err_lanczos, 4)))
})

test_that("Golden Lanczos filter strategy='lanczos' in filter_signal works", {
  ref <- skip_if_no_golden("golden_lanczos.json")

  g <- graph_ring(12)
  x <- ref$signal
  lmax_py <- ref$lmax
  scale <- ref$scale
  t_rgsp <- pygsp_heat_to_rgsp(scale, lmax_py)

  # Use filter_signal with explicit lanczos strategy
  y_via_strategy <- filter_signal(g, x, kernel_heat(t_rgsp), K = 30,
                                  lmax = lmax_py, strategy = "lanczos")

  expect_equal(as.numeric(y_via_strategy), ref$filtered_exact,
               tolerance = tol_filter_chebyshev)
})

test_that("Golden Lanczos vs Chebyshev: both match PyGSP Chebyshev output", {
  ref <- skip_if_no_golden("golden_lanczos.json")

  g <- graph_ring(12)
  x <- ref$signal
  lmax_py <- ref$lmax
  scale <- ref$scale
  t_rgsp <- pygsp_heat_to_rgsp(scale, lmax_py)

  y_cheby_r <- filter_signal(g, x, kernel_heat(t_rgsp), K = 30, lmax = lmax_py)

  # rgsp Chebyshev should match PyGSP Chebyshev closely
  expect_equal(as.numeric(y_cheby_r), ref$filtered_chebyshev_m30,
               tolerance = tol_filter_chebyshev)
})
