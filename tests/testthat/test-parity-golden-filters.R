# Golden-value parity tests: filter kernels and filtered output against PyGSP

test_that("Golden heat filter on Ring(12): kernel at eigenvalues", {
  ref <- skip_if_no_golden("golden_filters.json")
  entry <- ref$heat_ring12

  g <- graph_ring(12)
  eig <- sort(graph_eigenpairs(g)$values)
  lmax_py <- entry$lmax
  scale <- entry$scale

  # PyGSP: exp(-scale * lambda / lmax)
  # rgsp:  exp(-t * lambda) with t = scale / lmax
  t_rgsp <- pygsp_heat_to_rgsp(scale, lmax_py)
  kern <- kernel_heat(t_rgsp)
  kernel_r <- kern(eig)

  expect_equal(kernel_r, entry$kernel_at_eigenvalues, tolerance = tol_filter_exact)
})

test_that("Golden heat filter on Ring(12): filtered output matches PyGSP", {
  ref <- skip_if_no_golden("golden_filters.json")
  entry <- ref$heat_ring12

  g <- graph_ring(12)
  x <- entry$signal
  lmax_py <- entry$lmax
  scale <- entry$scale
  t_rgsp <- pygsp_heat_to_rgsp(scale, lmax_py)

  # Use exact method via eigendecomposition for best comparison
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  kern_vals <- kernel_heat(t_rgsp)(eig$values)
  y_exact <- as.numeric(U %*% diag(kern_vals) %*% t(U) %*% x)

  expect_equal(y_exact, entry$filtered_output, tolerance = tol_filter_exact)
})

test_that("Golden heat filter on Ring(12): Chebyshev filtered output matches PyGSP", {
  ref <- skip_if_no_golden("golden_filters.json")
  entry <- ref$heat_ring12

  g <- graph_ring(12)
  x <- entry$signal
  lmax_py <- entry$lmax
  scale <- entry$scale
  t_rgsp <- pygsp_heat_to_rgsp(scale, lmax_py)

  y_cheby <- filter_signal(g, x, kernel_heat(t_rgsp), K = 50, lmax = lmax_py)

  # Chebyshev at K=50 should be very close to exact
  expect_equal(as.numeric(y_cheby), entry$filtered_output, tolerance = tol_filter_chebyshev)
})

test_that("Golden heat filter on Grid2d(4,4): filtered output matches PyGSP", {
  ref <- skip_if_no_golden("golden_filters.json")
  entry <- ref$heat_grid2d_4_4

  g <- graph_grid2d(4, 4)
  x <- entry$signal
  lmax_py <- entry$lmax
  scale <- entry$scale
  t_rgsp <- pygsp_heat_to_rgsp(scale, lmax_py)

  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  kern_vals <- kernel_heat(t_rgsp)(eig$values)
  y_exact <- as.numeric(U %*% diag(kern_vals) %*% t(U) %*% x)

  expect_equal(y_exact, entry$filtered_output, tolerance = tol_filter_exact)
})

test_that("Golden Chebyshev filtered output matches PyGSP", {
  ref <- skip_if_no_golden("golden_chebyshev.json")

  g <- graph_ring(12)
  x <- ref$signal
  lmax_py <- ref$lmax
  scale <- ref$scale
  t_rgsp <- pygsp_heat_to_rgsp(scale, lmax_py)

  # Compare with Chebyshev at same order
  y_r <- filter_signal(g, x, kernel_heat(t_rgsp), K = 30, lmax = lmax_py)

  expect_equal(as.numeric(y_r), ref$filtered_chebyshev_m30, tolerance = tol_filter_chebyshev)

  # Exact filtered output should match closely
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  kern_vals <- kernel_heat(t_rgsp)(eig$values)
  y_exact_r <- as.numeric(U %*% diag(kern_vals) %*% t(U) %*% x)

  expect_equal(y_exact_r, ref$filtered_exact, tolerance = tol_filter_exact)
})
