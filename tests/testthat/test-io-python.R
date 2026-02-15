# NetworkX / graph-tool I/O tests (optional; skip if Python deps missing)

skip_if_no_networkx <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    testthat::skip("reticulate not installed")
  }
  if (!reticulate::py_available(initialize = FALSE)) {
    testthat::skip("Python not available for reticulate tests")
  }
  if (!reticulate::py_module_available("networkx")) {
    testthat::skip("networkx not available")
  }
  if (!reticulate::py_module_available("scipy.sparse")) {
    testthat::skip("scipy not available")
  }
}

test_that("NetworkX round-trip preserves adjacency", {
  skip_if_no_networkx()

  g <- graph_ring(20)
  g$coords <- cbind(seq_len(g$n), rep(0, g$n))

  g_nx <- graph_to_networkx(g, include_coords = TRUE)
  g2 <- graph_from_networkx(g_nx, normalized = g$normalized, include_coords = TRUE)

  expect_s3_class(g2, "gsp_graph")
  expect_equal(Matrix::drop0(g2$adjacency), Matrix::drop0(g$adjacency))
  expect_equal(as.matrix(g2$coords), as.matrix(g$coords))
})

