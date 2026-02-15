# Golden-value parity tests: graph eigenvalues against PyGSP
# Runs without Python -- compares against pre-generated JSON reference values

test_that("Golden eigenvalues match PyGSP for deterministic graphs", {
  ref <- skip_if_no_golden("golden_graphs.json")

  for (name in names(golden_graph_constructors)) {
    if (is.null(ref[[name]])) next
    g <- golden_graph_constructors[[name]]()
    eig <- sort(graph_eigenpairs(g)$values)
    ref_eig <- ref[[name]]$eigenvalues

    expect_equal(eig, ref_eig, tolerance = tol_eigenvalues,
                 label = paste("Eigenvalues:", name))
  }
})

test_that("Golden node counts match PyGSP for deterministic graphs", {
  ref <- skip_if_no_golden("golden_graphs.json")

  for (name in names(golden_graph_constructors)) {
    if (is.null(ref[[name]])) next
    g <- golden_graph_constructors[[name]]()
    expect_equal(g$n, ref[[name]]$n,
                 label = paste("Node count:", name))
  }
})

test_that("Golden lmax_exact matches PyGSP for deterministic graphs", {
  ref <- skip_if_no_golden("golden_graphs.json")

  for (name in names(golden_graph_constructors)) {
    if (is.null(ref[[name]])) next
    g <- golden_graph_constructors[[name]]()
    lmax_r <- lambda_max(g)
    lmax_exact_py <- ref[[name]]$lmax_exact

    expect_equal(lmax_r, lmax_exact_py, tolerance = tol_eigenvalues,
                 label = paste("lmax_exact:", name))
  }
})

test_that("Golden GFT abs-sorted coefficients match PyGSP", {
  ref <- skip_if_no_golden("golden_gft.json")

  g <- graph_ring(8)
  x <- ref$signal
  gft_r <- gft(g, x)

  expect_abs_sorted_equal(gft_r, ref$gft_abs_sorted, tolerance = tol_gft)
})

test_that("Golden GFT round-trip reconstructs signal", {
  ref <- skip_if_no_golden("golden_gft.json")

  g <- graph_ring(8)
  x <- ref$signal
  gft_r <- gft(g, x)
  rec <- igft(g, gft_r)

  expect_equal(as.numeric(rec), x, tolerance = tol_gft)
})
