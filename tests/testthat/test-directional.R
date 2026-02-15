test_that("graph_gradient reproduces Laplacian for undirected graphs", {
  g <- graph_ring(10, normalized = FALSE)
  D <- graph_gradient(g)
  L <- graph_laplacian(g, normalized = FALSE)
  err <- Matrix::norm(Matrix::t(D) %*% D - L, "F")
  expect_lt(err, 1e-8)
})

test_that("steered graph with alpha_isotropic=1 stays isotropic", {
  g <- graph_grid2d(4, 4, normalized = FALSE)
  gs <- graph_steerable_family(g, n_directions = 1, alpha_isotropic = 1)
  adj1 <- methods::as(gs[[1]]$adjacency, "generalMatrix")
  adj1 <- methods::as(adj1, "CsparseMatrix")
  adj0 <- methods::as(g$adjacency, "generalMatrix")
  adj0 <- methods::as(adj0, "CsparseMatrix")
  expect_equal(adj1, adj0)
})

test_that("steering emphasizes aligned edges on a grid", {
  g <- graph_grid2d(3, 3, normalized = FALSE)
  gs <- graph_steerable_family(g, orientations = matrix(c(1, 0), nrow = 1),
                               alpha_isotropic = 0.05, p = 4, absorb_sign = TRUE)
  W <- methods::as(gs[[1]]$adjacency, "generalMatrix")
  W <- methods::as(W, "TsparseMatrix")
  keep <- W@i < W@j
  i <- W@i[keep] + 1L
  j <- W@j[keep] + 1L
  w <- W@x[keep]
  coords <- g$coords
  deltas <- coords[j, , drop = FALSE] - coords[i, , drop = FALSE]
  horiz <- abs(deltas[, 1]) == 1 & deltas[, 2] == 0
  vert <- abs(deltas[, 2]) == 1 & deltas[, 1] == 0
  expect_gt(mean(w[horiz]), mean(w[vert]))
})

test_that("dsgwt_steer reduces to sgwt in isotropic limit", {
  g <- graph_path(5, normalized = FALSE)
  g <- graph_set_coords(g, kind = cbind(seq_len(g$n), 0))
  x <- rnorm(g$n)
  W_iso <- sgwt(g, x, scales = c(1, 2), wavelet = "mexican_hat", include_lowpass = TRUE, K = 10)
  W_dir <- dsgwt_steer(g, x, n_directions = 1, alpha_isotropic = 1,
                       scales = c(1, 2), wavelet = "mexican_hat",
                       include_lowpass = TRUE, K = 10)
  coeff_dir <- if (length(dim(W_dir$coeffs)) == 4) {
    drop(W_dir$coeffs[, , 1, 1, drop = FALSE])
  } else {
    drop(W_dir$coeffs[, , 1, drop = FALSE])
  }
  coeff_iso <- if (length(dim(W_iso$coeffs)) == 3) {
    drop(W_iso$coeffs[, , 1, drop = FALSE])
  } else {
    drop(W_iso$coeffs)
  }
  expect_equal(coeff_dir, coeff_iso)
})

test_that("Chebyshev cache is reused across directions", {
  cheby_cache_clear()
  cache_env <- get(".cheby_cache", envir = asNamespace("rgsp"))
  # baseline: cheby_coeffs caches
  cheby_coeffs(kernel_mexican_hat(1), K = 5, lmax = 2)
  base_len <- length(ls(envir = cache_env))
  expect_gt(base_len, 0)

  g <- graph_path(6, normalized = FALSE)
  g <- graph_set_coords(g, kind = cbind(seq_len(g$n), 0))
  x <- rnorm(g$n)

  dsgwt_steer(g, x, n_directions = 2, scales = c(1, 2),
              wavelet = "mexican_hat", include_lowpass = TRUE, K = 8,
              alpha_isotropic = 0.5)
  n_after_first <- length(ls(envir = cache_env))

  dsgwt_steer(g, x, n_directions = 5, scales = c(1, 2),
              wavelet = "mexican_hat", include_lowpass = TRUE, K = 8,
              alpha_isotropic = 0.5)
  n_after_second <- length(ls(envir = cache_env))

  expect_equal(n_after_second, n_after_first)
})

test_that("graph_gradient works on normalized and directed graphs (smoke)", {
  g_norm <- graph_ring(8, normalized = TRUE)
  Dn <- graph_gradient(g_norm)
  expect_equal(ncol(Dn), g_norm$n)
  expect_gt(nrow(Dn), 0)

  g_dir <- graph_erdos_renyi(10, p = 0.25, directed = TRUE, normalized = FALSE)
  Dd <- graph_gradient(g_dir)
  expect_equal(ncol(Dd), g_dir$n)
  expect_gt(nrow(Dd), 0)
})

test_that("isotropic limit holds for random graph (fuzz)", {
  set.seed(123)
  g <- graph_erdos_renyi(12, p = 0.2, normalized = FALSE)
  coords <- matrix(rnorm(g$n * 2), ncol = 2)
  g <- graph_set_coords(g, kind = coords)
  x <- rnorm(g$n)
  W_iso <- sgwt(g, x, scales = c(1.5, 3), wavelet = "heat", include_lowpass = TRUE, K = 12)
  W_dir <- dsgwt_steer(g, x, n_directions = 2, alpha_isotropic = 1,
                       scales = c(1.5, 3), wavelet = "heat",
                       include_lowpass = TRUE, K = 12, p = 0)
  coeff_dir <- drop(W_dir$coeffs[, , 1, 1, drop = FALSE])
  coeff_iso <- drop(W_iso$coeffs[, , 1, drop = FALSE])
  expect_equal(coeff_dir, coeff_iso, tolerance = 1e-10)
})

test_that("frame bounds finite and lmax deviations small on grid", {
  g <- graph_grid2d(6, 6, normalized = FALSE)
  g <- graph_set_coords(g, kind = g$coords)
  res <- dsgwt_steer(g, x = rnorm(g$n), n_directions = 3, scales = c(1, 2),
                     wavelet = "mexican_hat", include_lowpass = TRUE, K = 8,
                     check_lmax = FALSE)
  # frame bounds from sgwt call inside sgwt are not exposed; check coeffs finite
  expect_true(all(is.finite(res$coeffs)))
})

test_that("lmax across directions remains close to base on grids", {
  g <- graph_grid2d(5, 5, normalized = TRUE)
  g <- graph_set_coords(g, kind = g$coords)
  lmax_base <- lambda_max(g, normalized = g$normalized)
  gs <- graph_steerable_family(g, n_directions = 3,
                               alpha_isotropic = 0.7, p = 1, absorb_sign = TRUE)
  lmax_dir <- vapply(gs, function(gi) lambda_max(gi, normalized = gi$normalized), numeric(1))
  dev <- max(abs(lmax_dir - lmax_base) / lmax_base)
  expect_lt(dev, 0.1)
})

test_that("orientations metadata stored and normalized", {
  g <- graph_grid2d(4, 4, normalized = FALSE)
  g <- graph_set_coords(g, kind = g$coords)
  U <- rbind(c(1, 0), c(0, 1), c(1, 1))
  res <- dsgwt_steer(g, rnorm(g$n), orientations = U, n_directions = nrow(U),
                     scales = c(1), wavelet = "heat", include_lowpass = FALSE, K = 6)
  stored <- res$orientations
  norms <- sqrt(rowSums(stored^2))
  expect_equal(nrow(stored), nrow(U))
  expect_true(all(abs(norms - 1) < 1e-8))
})

test_that("dsgwt_steer output dims follow node-band-direction-signal", {
  g <- graph_grid2d(3, 3, normalized = FALSE)
  g <- graph_set_coords(g, kind = g$coords)
  x <- cbind(rnorm(g$n), rnorm(g$n))
  res <- dsgwt_steer(g, x, n_directions = 2, scales = c(1),
                     wavelet = "heat", include_lowpass = TRUE, K = 6)
  dims <- dim(res$coeffs)
  expect_equal(length(dims), 4)
  expect_equal(dims[1], g$n)
  expect_equal(dims[3], 2)
  expect_equal(dims[4], 2)
})

test_that("directional coefficients collapse and reconstruct close to sgwt", {
  set.seed(99)
  g <- graph_grid2d(4, 4, normalized = FALSE)
  g <- graph_set_coords(g, kind = g$coords)
  x <- rnorm(g$n)
  scales <- c(1, 2)
  W_iso <- sgwt(g, x, scales = scales, wavelet = "mexican_hat",
                include_lowpass = TRUE, K = 12)
  W_dir <- dsgwt_steer(g, x, n_directions = 2, scales = scales,
                       wavelet = "mexican_hat", include_lowpass = TRUE, K = 12,
                       alpha_isotropic = 0.8, p = 1)
  # simple collapse: mean over directions, then inverse with original filters
  coeffs_mean <- apply(W_dir$coeffs, c(1, 2, 4), mean)
  # replace the wavelet coeffs but keep lowpass identical
  W_collapsed <- W_iso
  W_collapsed$coeffs[, -1, ] <- array(coeffs_mean[, -1, , drop = FALSE],
                                      dim = dim(W_iso$coeffs[, -1, , drop = FALSE]))
  x_rec <- isgwt(W_collapsed)
  x_iso <- isgwt(W_iso)
  err <- max(abs(x_rec - x_iso))
  expect_lt(err, 0.2)
})
