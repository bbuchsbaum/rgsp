set.seed(123)

test_that("laplacian is PSD and annihilates constants (small graphs)", {
  graphs <- list(
    graph_ring(6),
    graph_grid2d(2, 3),
    graph_erdos_renyi(6, p = 0.4, seed = 2)
  )
  for (g in graphs) {
    L <- graph_laplacian(g, normalized = FALSE)
    eig <- eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values
    expect_true(min(eig) >= -1e-8)
    ones <- rep(1, g$n)
    expect_equal(as.numeric(L %*% ones), rep(0, g$n), tolerance = 1e-8)
  }
})

test_that("degree sum equals twice edge count for undirected graphs", {
  gs <- list(graph_ring(5), graph_grid2d(2, 2), graph_complete(4))
  for (g in gs) {
    deg_sum <- sum(graph_degree(g))
    edge_count <- length(g$adjacency@x) / 2
    expect_equal(deg_sum, 2 * edge_count)
  }
})

test_that("gft and igft are inverses", {
  g <- graph_grid2d(2, 2)
  x <- rnorm(g$n)
  X <- gft(g, x)
  xr <- igft(g, X)
  expect_equal(drop(xr), x, tolerance = 1e-8)
})

test_that("heat propagation is contractive in l2", {
  g <- graph_ring(8)
  x <- rnorm(8)
  y <- heat_propagate(g, x, t = 0.3)
  expect_true(sum(y^2) <= sum(x^2) + 1e-8)
})

test_that("lambda_max scales with graph_rescale", {
  g <- graph_ring(8)
  l1 <- lambda_max(g)
  g2 <- graph_rescale(g, 3)
  l2 <- lambda_max(g2)
  expect_true(l2 > l1 * 2.9)
})

test_that("random graph seeds produce identical adjacency", {
  g1 <- graph_erdos_renyi(10, p = 0.2, seed = 11)
  g2 <- graph_erdos_renyi(10, p = 0.2, seed = 11)
  expect_equal(g1$adjacency, g2$adjacency)
})

test_that("fuzz: small ER graphs remain PSD", {
  for (i in 1:5) {
    g <- graph_erdos_renyi(6, p = runif(1, 0.1, 0.6), seed = i)
    L <- graph_laplacian(g, normalized = FALSE)
    eig_min <- min(eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values)
    expect_true(eig_min >= -1e-8)
  }
})

test_that("kernel functions stay bounded between 0 and 1 where expected", {
  lam <- seq(0, 4, length.out = 20)
  expect_true(all(kernel_half_cosine(1, 2)(lam) <= 1 + 1e-12))
  expect_true(all(kernel_rectangle(0.5, 2.5)(lam) %in% c(0, 1)))
  expect_true(all(kernel_wave(0.3)(lam) <= 1))
})

test_that("spectrogram energies are nonnegative", {
  g <- graph_ring(7)
  x <- rnorm(7)
  S <- compute_spectrogram(g, x, scales = c(0.4, 0.8), K = 10)
  expect_true(all(S >= 0))
})

test_that("kron_reduction preserves PSD", {
  g <- graph_grid2d(3, 3)
  g2 <- kron_reduction(g, keep = c(1, 3, 5, 7, 9))
  eig <- eigen(as.matrix(graph_laplacian(g2, normalized = FALSE)), symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(eig) >= -1e-8)
})
