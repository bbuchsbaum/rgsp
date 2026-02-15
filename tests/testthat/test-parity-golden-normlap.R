# Golden-value parity tests: normalized Laplacian eigenvalues against PyGSP

test_that("Golden normalized Laplacian eigenvalues: Ring(10)", {
  ref <- skip_if_no_golden("golden_normalized_laplacian.json")
  entry <- ref$Ring_10

  g <- graph_ring(10)
  L_norm <- graph_laplacian(g, normalized = TRUE)
  eigs_r <- sort(eigen(as.matrix(L_norm), symmetric = TRUE, only.values = TRUE)$values)

  expect_equal(eigs_r, entry$eigenvalues_normalized, tolerance = tol_eigenvalues)
})

test_that("Golden normalized Laplacian eigenvalues: Grid2d(3,4)", {
  ref <- skip_if_no_golden("golden_normalized_laplacian.json")
  entry <- ref$Grid2d_3_4

  g <- graph_grid2d(3, 4)
  L_norm <- graph_laplacian(g, normalized = TRUE)
  eigs_r <- sort(eigen(as.matrix(L_norm), symmetric = TRUE, only.values = TRUE)$values)

  expect_equal(eigs_r, entry$eigenvalues_normalized, tolerance = tol_eigenvalues)
})

test_that("Golden normalized Laplacian eigenvalues are in [0, 2]", {
  ref <- skip_if_no_golden("golden_normalized_laplacian.json")

  for (name in c("Ring_10", "Grid2d_3_4")) {
    entry <- ref[[name]]
    eigs <- entry$eigenvalues_normalized
    expect_true(all(eigs >= -1e-10),
                label = paste(name, "min eigenvalue >= 0"))
    expect_true(all(eigs <= 2 + 1e-10),
                label = paste(name, "max eigenvalue <= 2"))
  }
})

test_that("Normalized and combinatorial Laplacian eigenvalues differ", {
  ref_norm <- skip_if_no_golden("golden_normalized_laplacian.json")
  ref_graphs <- skip_if_no_golden("golden_graphs.json")

  # Ring(10) is regular so normalized eigenvalues = combinatorial / degree
  g <- graph_ring(10)
  eigs_comb <- ref_graphs$Ring_10$eigenvalues
  eigs_norm <- ref_norm$Ring_10$eigenvalues_normalized

  # For 2-regular ring: L_norm = L/2, so eigenvalues should be halved
  expect_equal(eigs_norm, eigs_comb / 2, tolerance = tol_eigenvalues)
})
