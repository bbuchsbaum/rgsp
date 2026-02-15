#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(reticulate)
})

here <- function(...) file.path(getwd(), ...)

out_dir <- here("inst", "extdata", "pointclouds")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scipy_io <- import("scipy.io", delay_load = FALSE)

as_dgC <- function(sp) {
  if (inherits(sp, "dgCMatrix")) {
    return(Matrix::drop0(sp))
  }
  if (inherits(sp, "Matrix")) {
    return(Matrix::drop0(methods::as(sp, "dgCMatrix")))
  }
  coo <- sp$tocoo()
  n <- as.integer(sp$shape[[1]])
  i <- as.integer(coo$row) + 1L
  j <- as.integer(coo$col) + 1L
  x <- as.numeric(coo$data)
  Matrix::drop0(Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n)))
}

py_radius_graph <- local({
  # Build a radius (epsilon-NN) graph in Python (KDTree) with PyGSP's
  # similarity weights exp(-d^2 / sigma) and average symmetrization.
  reticulate::py_run_string(
    paste(
      "import numpy as np",
      "from scipy import sparse, spatial",
      "from scipy.spatial import distance",
      "",
      "def _rgsp_radius_graph(X, epsilon, sigma=None):",
      "    X = np.asarray(X, dtype=float)",
      "    kdt = spatial.KDTree(X)",
      "    NN = [kdt.query_ball_point(point, r=float(epsilon)) for point in X]",
      "    D = []",
      "    for i, neighbors in enumerate(NN):",
      "        if len(neighbors) > 0:",
      "            D.append([distance.minkowski(X[i], X[j], p=2) for j in neighbors])",
      "        else:",
      "            D.append([])",
      "",
      "    if sigma is None:",
      "        all_dist = []",
      "        for i, neighbors in enumerate(NN):",
      "            for idx, j in enumerate(neighbors):",
      "                if j != i:",
      "                    all_dist.append(D[i][idx])",
      "        if len(all_dist) == 0:",
      "            raise ValueError('No neighbors found')",
      "        sigma = float(np.mean(all_dist))",
      "",
      "    rows = []",
      "    cols = []",
      "    vals = []",
      "    for i, neighbors in enumerate(NN):",
      "        for idx, j in enumerate(neighbors):",
      "            if j == i:",
      "                continue",
      "            rows.append(i)",
      "            cols.append(j)",
      "            d = D[i][idx]",
      "            vals.append(np.exp(-(d * d) / float(sigma)))",
      "",
      "    W = sparse.csc_matrix((vals, (rows, cols)), shape=(X.shape[0], X.shape[0]))",
      "    W = (W + W.T) / 2.0",
      "    return W, sigma",
      sep = "\n"
    )
  )

  function(X, epsilon, sigma = NULL) {
    reticulate::py$`_rgsp_radius_graph`(X, epsilon, sigma)
  }
})

write_rds <- function(name, obj) {
  path <- file.path(out_dir, paste0(name, ".rds"))
  saveRDS(obj, path)
  message("wrote: ", path)
}

load_mat <- function(rel_path) {
  scipy_io$loadmat(here(rel_path))
}

# Logo -----------------------------------------------------------------------

logo <- load_mat("pygsp/pygsp/data/pointclouds/logogsp.mat")
W <- as_dgC(logo$W)
coords <- as.matrix(logo$coords)

write_rds("logogsp", list(
  W = W,
  coords = coords,
  info = list(
    idx_g = as.integer(logo$idx_g),
    idx_s = as.integer(logo$idx_s),
    idx_p = as.integer(logo$idx_p)
  ),
  plotting = list(limits = c(0, 640, -400, 0))
))

# Minnesota ------------------------------------------------------------------

mn <- load_mat("pygsp/pygsp/data/pointclouds/minnesota.mat")
A <- as_dgC(mn$A)
xy <- as.matrix(mn$xy)
labels <- as.character(mn$labels)

write_rds("minnesota", list(
  A = A,
  xy = xy,
  labels = labels,
  plotting = list(limits = c(-98, -89, 43, 50), vertex_size = 40)
))

# Airfoil --------------------------------------------------------------------

air <- load_mat("pygsp/pygsp/data/pointclouds/airfoil.mat")

coords <- cbind(as.numeric(air$x), as.numeric(air$y))
i_inds <- as.integer(air$i_inds)
j_inds <- as.integer(air$j_inds)

A <- Matrix::sparseMatrix(i = i_inds, j = j_inds, x = 1, dims = c(4253, 4253))
W <- (A + Matrix::t(A)) / 2
W <- Matrix::drop0(W)

write_rds("airfoil", list(
  W = W,
  coords = coords,
  plotting = list(
    vertex_size = 30,
    limits = c(-1e-4, 1.01 * max(coords[, 1]), -1e-4, 1.01 * max(coords[, 2]))
  )
))

# Bunny ----------------------------------------------------------------------

bunny <- load_mat("pygsp/pygsp/data/pointclouds/bunny.mat")
X <- as.matrix(bunny$bunny)

res <- py_radius_graph(X, epsilon = 0.02, sigma = NULL)
W <- as_dgC(res[[1]])

write_rds("bunny", list(
  W = W,
  coords = X,
  plotting = list(vertex_size = 10, elevation = -90, azimuth = 90, distance = 8),
  params = list(epsilon = 0.02, sigma = as.numeric(res[[2]]))
))

# Two Moons ------------------------------------------------------------------

tm <- load_mat("pygsp/pygsp/data/pointclouds/two_moons.mat")
features <- as.matrix(tm$features)

write_rds("two_moons", list(
  features = features
))

# DavidSensorNet --------------------------------------------------------------

d64 <- load_mat("pygsp/pygsp/data/pointclouds/david64.mat")
write_rds("david64", list(
  N = as.integer(d64$N[1, 1]),
  W = as_dgC(d64$W),
  coords = as.matrix(d64$coords),
  plotting = list(limits = c(0, 1, 0, 1))
))

d500 <- load_mat("pygsp/pygsp/data/pointclouds/david500.mat")
write_rds("david500", list(
  N = as.integer(d500$N[1, 1]),
  W = as_dgC(d500$W),
  coords = as.matrix(d500$coords),
  plotting = list(limits = c(0, 1, 0, 1))
))
