test_that("grad matches incidence product", {
  set.seed(1)
  g <- graph_ring(6)
  B <- graph_incidence(g)$B
  x <- rnorm(6)
  g_grad <- grad(g, x)
  expect_equal(as.numeric(g_grad), as.numeric(B %*% x))
})

test_that("div matches -t(B) product", {
  g <- graph_ring(5)
  inc <- graph_incidence(g)
  edge_sig <- rnorm(nrow(inc$B))
  g_div <- div(g, edge_sig)
  expect_equal(as.numeric(g_div), as.numeric(-Matrix::t(inc$B) %*% edge_sig))
})

test_that("compute_differential_operator returns working closures", {
  g <- graph_ring(4)
  diffops <- compute_differential_operator(g)
  x <- 1:4
  grad_x <- diffops$grad(x)
  div_grad <- diffops$div(grad_x)
  # div(grad(x)) should sum to ~0 on undirected regular graphs
  expect_true(abs(sum(div_grad)) < 1e-8)
})
