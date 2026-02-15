test_that("graph_sparsify_spectral returns connected sparser graph", {
  g <- graph_sensor(200, k = 10, seed = 1)
  e_in <- length(rgsp:::.edge_list_undirected(g$adjacency)$w)

  gs <- graph_sparsify_spectral(g, epsilon = 0.5, seed = 1, n_probes = 24)
  expect_s3_class(gs, "gsp_graph")
  e_out <- length(rgsp:::.edge_list_undirected(gs$adjacency)$w)

  expect_true(rgsp:::.graph_is_connected(gs$adjacency))
  expect_lt(e_out, e_in)

  # Filtering should still work
  x <- rnorm(g$n)
  y <- filter_signal(gs, x, kernel_heat(1), K = 20)
  expect_equal(length(y), g$n)
})

