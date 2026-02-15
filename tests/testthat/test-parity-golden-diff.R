# Golden-value parity tests: differential operators against PyGSP
#
# NOTE: PyGSP D is N x E, rgsp D is E x N (transposed convention).
# We compare edge count, node count, and the div(grad(x))=Lx identity.

test_that("Golden gradient dimensions match PyGSP for Ring(10)", {
  ref <- skip_if_no_golden("golden_differential.json")
  entry <- ref$ring10

  g <- graph_ring(10)
  D <- compute_differential_operator(g)$B

  # rgsp D is E x N; PyGSP D is N x E
  expect_equal(nrow(D), entry$n_edges)   # rows = edges in rgsp
  expect_equal(ncol(D), entry$n_nodes)   # cols = nodes in rgsp
})

test_that("Golden gradient dimensions match PyGSP for Grid2d(3,3)", {
  ref <- skip_if_no_golden("golden_differential.json")
  entry <- ref$grid2d_3_3

  g <- graph_grid2d(3, 3)
  D <- compute_differential_operator(g)$B

  # rgsp D is E x N; PyGSP D is N x E (transposed)
  expect_equal(nrow(D), entry$n_edges)   # rows = edges in rgsp
  expect_equal(ncol(D), entry$n_nodes)   # cols = nodes in rgsp
})

test_that("Golden div(grad(x)) = L*x identity for Grid2d(3,3)", {
  ref <- skip_if_no_golden("golden_differential.json")
  entry <- ref$grid2d_3_3

  g <- graph_grid2d(3, 3)
  x <- entry$signal

  # rgsp: div(grad(x)) = -L*x  (sign convention)
  grad_r <- grad(g, x)
  div_grad_r <- div(g, grad_r)
  L <- graph_laplacian(g)
  Lx <- as.numeric(L %*% x)

  expect_equal(as.numeric(div_grad_r), -Lx, tolerance = tol_differential)

  # Also verify Lx matches PyGSP
  expect_equal(Lx, entry$Lx, tolerance = tol_differential)
})

test_that("Golden gradient output has correct length for Ring(10)", {
  ref <- skip_if_no_golden("golden_differential.json")
  entry <- ref$ring10

  g <- graph_ring(10)
  x <- entry$signal
  grad_r <- grad(g, x)

  expect_equal(length(grad_r), entry$n_edges)
})
