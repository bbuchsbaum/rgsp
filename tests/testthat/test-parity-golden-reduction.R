# Golden-value parity tests: Kron reduction against PyGSP

test_that("Golden Kron reduction: node count matches PyGSP", {
  ref <- skip_if_no_golden("golden_reduction.json")
  entry <- ref$kron_ring12

  g <- graph_ring(12)
  keep <- entry$keep_1indexed  # 1-indexed for R
  g_red <- kron_reduction(g, keep)

  expect_equal(g_red$n, entry$reduced_n)
})

test_that("Golden Kron reduction: eigenvalues match PyGSP", {
  ref <- skip_if_no_golden("golden_reduction.json")
  entry <- ref$kron_ring12

  g <- graph_ring(12)
  keep <- entry$keep_1indexed
  g_red <- kron_reduction(g, keep)

  eig_r <- sort(graph_eigenpairs(g_red)$values)
  eig_py <- entry$reduced_eigenvalues

  expect_equal(eig_r, eig_py, tolerance = tol_kron)
})

test_that("Golden Kron reduction: reduced Laplacian is PSD", {
  ref <- skip_if_no_golden("golden_reduction.json")
  entry <- ref$kron_ring12

  g <- graph_ring(12)
  keep <- entry$keep_1indexed
  g_red <- kron_reduction(g, keep)

  L_red <- graph_laplacian(g_red, normalized = FALSE)
  eig_vals <- eigen(as.matrix(L_red), symmetric = TRUE, only.values = TRUE)$values

  expect_true(min(eig_vals) >= -1e-8)
})

test_that("Golden Kron reduction: constants in null space", {
  ref <- skip_if_no_golden("golden_reduction.json")
  entry <- ref$kron_ring12

  g <- graph_ring(12)
  keep <- entry$keep_1indexed
  g_red <- kron_reduction(g, keep)

  L_red <- graph_laplacian(g_red, normalized = FALSE)
  ones <- rep(1, g_red$n)

  expect_equal(as.numeric(L_red %*% ones), rep(0, g_red$n), tolerance = 1e-10)
})
