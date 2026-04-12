# Graph constructors and helpers

#' Create a graph object from adjacency
#' @param adjacency sparse adjacency (dgCMatrix)
#' @param coords optional node coordinates (matrix/data.frame)
#' @param normalized logical; if TRUE compute normalized Laplacian
#' @param directed logical; if TRUE, Laplacian uses directed form (random-walk)
#' @export
new_graph <- function(adjacency, coords = NULL, normalized = FALSE, directed = FALSE) {
  if (methods::is(adjacency, "lMatrix")) {
    # avoid deprecated lgC -> dgC coercion warnings
    adjacency <- Matrix::drop0(adjacency + 0)
  } else if (methods::is(adjacency, "dsCMatrix")) {
    # recommended path: symmetric -> general -> Csparse
    adjacency <- methods::as(adjacency, "generalMatrix")
    adjacency <- methods::as(adjacency, "CsparseMatrix")
  } else if (!inherits(adjacency, "dgCMatrix")) {
    adjacency <- methods::as(adjacency, "dgCMatrix")
  }
  if (nrow(adjacency) != ncol(adjacency)) {
    stop("Adjacency must be square")
  }
  # enforce symmetry within tolerance
  if (!directed && !isTRUE(Matrix::isSymmetric(adjacency, tol = 1e-9))) {
    stop("Adjacency must be symmetric unless directed=TRUE")
  }

  deg <- degree_cpp(adjacency)
  lap <- laplacian_cpp(adjacency, normalized = normalized, directed = directed)
  cache_env <- new.env(parent = emptyenv())
  cache_env$degree <- deg
  cache_env[[paste0("lap_", normalized)]] <- lap

  g <- list(
    adjacency = adjacency,
    laplacian = lap,
    degree = deg,
    n = as.integer(nrow(adjacency)),
    coords = coords,
    normalized = normalized,
    directed = directed,
    cache = cache_env
  )
  class(g) <- "gsp_graph"
  g
}

#' Ring (or path) graph
#' @param n number of nodes
#' @param periodic if FALSE, returns a path graph
#' @param normalized logical; compute normalized Laplacian
#' @export
graph_ring <- function(n, periodic = TRUE, normalized = FALSE) {
  adj <- adjacency_ring_cpp(as.integer(n), periodic = periodic)
  g <- new_graph(adj, normalized = normalized, directed = FALSE)
  g$graph_type <- if (isTRUE(periodic)) "ring" else "path"
  g$graph_params <- list(periodic = isTRUE(periodic), n = as.integer(n))
  g
}

#' Path graph (alias)
#' @param n number of nodes
#' @param normalized logical; compute normalized Laplacian
#' @return a gsp_graph object
#' @export
graph_path <- function(n, normalized = FALSE) {
  g <- graph_ring(n, periodic = FALSE, normalized = normalized)
  g$graph_type <- "path"
  g$graph_params <- list(periodic = FALSE, n = as.integer(n))
  g
}

#' 3D grid graph (optionally torus)
#' @param nx number of nodes in x
#' @param ny number of nodes in y
#' @param nz number of nodes in z
#' @param periodic wrap edges to make a 3D torus
#' @param normalized logical; compute normalized Laplacian
#' @export
graph_grid3d <- function(nx, ny, nz, periodic = FALSE, normalized = FALSE) {
  adj <- adjacency_grid3d_cpp(as.integer(nx), as.integer(ny), as.integer(nz), periodic = periodic)
  coords <- as.matrix(expand.grid(x = seq_len(nx), y = seq_len(ny), z = seq_len(nz)))
  new_graph(adj, coords = coords, normalized = normalized, directed = FALSE)
}

#' 2D grid graph (optionally torus)
#' @param nrow number of rows
#' @param ncol number of columns
#' @param periodic wrap edges to make a torus
#' @param normalized logical; compute normalized Laplacian
#' @export
graph_grid2d <- function(nrow, ncol, periodic = FALSE, normalized = FALSE) {
  adj <- adjacency_grid2d_cpp(as.integer(nrow), as.integer(ncol), periodic = periodic)
  coords <- as.matrix(expand.grid(x = seq_len(ncol), y = seq_len(nrow)))
  g <- new_graph(adj, coords = coords, normalized = normalized, directed = FALSE)
  g$graph_type <- "grid2d"
  g$graph_params <- list(nrow = as.integer(nrow), ncol = as.integer(ncol), periodic = isTRUE(periodic))
  g
}

#' Random geometric graph in 2D (unit square)
#' @param n nodes
#' @param radius connect nodes within this distance
#' @param normalized compute normalized Laplacian
#' @param seed optional RNG seed for reproducibility
#' @return a gsp_graph object
#' @export
graph_random_geometric <- function(n, radius, normalized = FALSE, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  coords <- cbind(runif(n), runif(n))
  distmat <- as.matrix(dist(coords))
  adj <- Matrix::Matrix(distmat <= radius & distmat > 0, sparse = TRUE)
  new_graph(adj, coords = coords, normalized = normalized, directed = FALSE)
}

#' Sensor graph (k-nearest neighbors)
#' @param n nodes
#' @param k neighbors per node
#' @param normalized compute normalized Laplacian
#' @param seed optional RNG seed for reproducibility
#' @return a gsp_graph object
#' @export
graph_sensor <- function(n, k = 6, normalized = FALSE, seed = NULL) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  coords <- cbind(runif(n), runif(n))
  distmat <- as.matrix(dist(coords))
  diag(distmat) <- Inf
  idx <- t(apply(distmat, 1, function(d) order(d)[seq_len(k)]))
  i <- rep(seq_len(n), each = k)
  j <- as.vector(idx)
  x <- rep(1, length(i))
  adj <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
  adj <- (adj + t(adj)) > 0
  Matrix::diag(adj) <- FALSE
  adj <- Matrix::Matrix(adj, sparse = TRUE)
  new_graph(adj, coords = coords, normalized = normalized, directed = FALSE)
}

# Built-in dataset graphs -----------------------------------------------------

.load_pointcloud_rds <- function(name) {
  filename <- paste0(name, ".rds")
  path <- system.file("extdata", "pointclouds", filename, package = "rgsp")
  if (path == "") {
    path <- file.path("inst", "extdata", "pointclouds", filename)
  }
  if (!file.exists(path)) {
    stop("Missing extdata file: ", filename, call. = FALSE)
  }
  readRDS(path)
}

#' GSP Logo graph (PyGSP dataset)
#'
#' @param normalized logical; compute normalized Laplacian
#' @return a `gsp_graph` with `coords` and `info` indices
#' @export
graph_logo <- function(normalized = FALSE) {
  d <- .load_pointcloud_rds("logogsp")
  g <- new_graph(d$W, coords = d$coords, normalized = normalized, directed = FALSE)
  g$info <- d$info
  g$plotting <- d$plotting
  g
}

#' Minnesota road network graph (PyGSP dataset)
#'
#' @param connected logical; if TRUE, match PyGSP's connected/binarized variant
#' @param normalized logical; compute normalized Laplacian
#' @return a `gsp_graph` with `coords` and `labels`
#' @export
graph_minnesota <- function(connected = TRUE, normalized = FALSE) {
  d <- .load_pointcloud_rds("minnesota")
  adj <- d$A

  if (isTRUE(connected)) {
    # Match PyGSP: add missing edge then binarize.
    adj_bin <- adj > 0
    adj_bin[349, 355] <- TRUE
    adj_bin[355, 349] <- TRUE
    Matrix::diag(adj_bin) <- FALSE
    adj <- Matrix::drop0(adj_bin + 0)
  } else {
    Matrix::diag(adj) <- 0
    adj <- Matrix::drop0(adj)
  }

  g <- new_graph(adj, coords = d$xy, normalized = normalized, directed = FALSE)
  g$labels <- d$labels
  g$plotting <- d$plotting
  g
}

#' Airfoil graph (PyGSP dataset)
#'
#' @param normalized logical; compute normalized Laplacian
#' @return a `gsp_graph` with `coords`
#' @export
graph_airfoil <- function(normalized = FALSE) {
  d <- .load_pointcloud_rds("airfoil")
  g <- new_graph(d$W, coords = d$coords, normalized = normalized, directed = FALSE)
  g$plotting <- d$plotting
  g
}

#' Stanford Bunny graph (PyGSP dataset)
#'
#' @param normalized logical; compute normalized Laplacian
#' @return a `gsp_graph` with 3D `coords`
#' @export
graph_bunny <- function(normalized = FALSE) {
  d <- .load_pointcloud_rds("bunny")
  g <- new_graph(d$W, coords = d$coords, normalized = normalized, directed = FALSE)
  g$plotting <- d$plotting
  g$params <- d$params
  g
}

.two_moons_arc <- function(N, sigmad, distance, number) {
  phi <- runif(N) * pi
  rb <- sigmad * rnorm(N)
  ab <- runif(N) * 2 * pi
  bx <- rb * cos(ab)
  by <- rb * sin(ab)

  if (number == 1L) {
    x <- cos(phi) + bx + 0.5
    y <- -sin(phi) + by - (distance - 1) / 2
  } else {
    x <- cos(phi) + bx - 0.5
    y <- sin(phi) + by + (distance - 1) / 2
  }

  cbind(x, y)
}

#' Two Moons graph (PyGSP dataset)
#'
#' @param moontype `"standard"` (from dataset) or `"synthesized"`
#' @param dim point-cloud dimensionality for `moontype="standard"` (default 2)
#' @param sigmag similarity bandwidth for weights (default 0.05)
#' @param N number of vertices for `moontype="synthesized"` (default 400)
#' @param sigmad noise variance for synthesized moons (default 0.07)
#' @param distance distance between moons for synthesized graph (default 0.5)
#' @param seed optional RNG seed for reproducibility (synthesized)
#' @param normalized logical; compute normalized Laplacian
#' @return a `gsp_graph` with `coords` and `labels`
#' @export
graph_two_moons <- function(moontype = c("standard", "synthesized"),
                            dim = 2,
                            sigmag = 0.05,
                            N = 400,
                            sigmad = 0.07,
                            distance = 0.5,
                            seed = NULL,
                            normalized = FALSE) {
  moontype <- match.arg(moontype)

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }

  if (moontype == "standard") {
    d <- .load_pointcloud_rds("two_moons")
    features <- as.matrix(d$features)
    dim <- as.integer(dim)[1]
    if (dim < 1L || dim > nrow(features)) stop("dim must be in [1, n_features]")
    coords <- t(features[seq_len(dim), , drop = FALSE])
    n_total <- nrow(coords)
    N1 <- n_total %/% 2L
    N2 <- n_total - N1
    labels <- c(rep(0L, N1), rep(1L, N2))
  } else {
    N <- as.integer(N)[1]
    if (N < 2L) stop("N must be >= 2")
    N1 <- N %/% 2L
    N2 <- N - N1
    coords1 <- .two_moons_arc(N1, sigmad, distance, 1L)
    coords2 <- .two_moons_arc(N2, sigmad, distance, 2L)
    coords <- rbind(coords1, coords2)
    labels <- c(rep(0L, N1), rep(1L, N2))
  }

  g <- graph_knn(coords,
                 k = 5,
                 weight = "pygsp",
                 sigma = sigmag,
                 sym = "average",
                 normalized = normalized)
  g$labels <- labels
  g$plotting <- list(vertex_size = 30)
  g
}

#' David sensor network graph (PyGSP dataset)
#'
#' For `N = 64` and `N = 500`, loads the precomputed dataset graphs.
#' Otherwise, generates a random sensor graph using the PyGSP heuristic.
#'
#' @param N number of vertices (default 64)
#' @param seed optional RNG seed
#' @param normalized logical; compute normalized Laplacian
#' @return a `gsp_graph` with `coords`
#' @export
graph_david_sensor_net <- function(N = 64, seed = NULL, normalized = FALSE) {
  N <- as.integer(N)[1]
  if (N < 2L) stop("N must be >= 2")

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }

  if (N %in% c(64L, 500L)) {
    d <- .load_pointcloud_rds(paste0("david", N))
    g <- new_graph(d$W, coords = d$coords, normalized = normalized, directed = FALSE)
    g$plotting <- d$plotting
    return(g)
  }

  coords <- cbind(runif(N), runif(N))

  target_dist_cutoff <- -0.125 * N / 436.075 + 0.2183
  T <- 0.6
  s <- sqrt(-(target_dist_cutoff^2) / (2 * log(T)))

  d <- as.matrix(dist(coords))
  W <- exp(-(d^2) / (2 * s^2))
  W[W < T] <- 0
  diag(W) <- 0

  adj <- Matrix::Matrix(W, sparse = TRUE)
  g <- new_graph(adj, coords = coords, normalized = normalized, directed = FALSE)
  g$plotting <- list(limits = c(0, 1, 0, 1))
  g
}

#' k-nearest-neighbor graph
#'
#' Build an undirected k-NN graph from node coordinates. Supports simple binary
#' weights, inverse-distance weights, or a heat kernel with an optional
#' bandwidth parameter. Uses `FNN::get.knn()` when available and falls back to a
#' distance matrix for very small graphs.
#'
#' @param coords numeric matrix of node coordinates (n x d)
#' @param k number of neighbors per node (clipped to `n-1`)
#' @param weight weighting scheme: `"heat"`, `"binary"`, `"distance"`, or `"pygsp"`
#' @param sigma optional bandwidth for `weight = "heat"`; defaults to the median
#'   k-NN distance
#' @param sym symmetrization mode: `"union"` (default), `"mutual"`, or `"average"`
#' @param normalized logical; if TRUE compute normalized Laplacian
#' @param seed optional RNG seed for reproducibility
#' @param eps small constant to avoid divide-by-zero for distance weights
#'
#' @return `gsp_graph` object
#' @importFrom FNN get.knn
#' @export
graph_knn <- function(coords, k = 6, weight = c("heat", "binary", "distance", "pygsp"),
                      sigma = NULL, sym = c("union", "mutual", "average"),
                      normalized = FALSE, seed = NULL, eps = 1e-8) {
  if (length(weight) > 1) weight <- weight[1]
  if (length(sym) > 1) sym <- sym[1]
  weight <- match.arg(weight)
  sym <- match.arg(sym)
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if (n < 2) stop("coords must have at least 2 rows")

  k <- max(1L, min(as.integer(k), n - 1L))

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }

  # nearest neighbours -------------------------------------------------------
  nn_idx <- NULL
  nn_dist <- NULL
  if (requireNamespace("FNN", quietly = TRUE) && n > 3) {
    res <- FNN::get.knn(coords, k = k)
    nn_idx <- res$nn.index
    nn_dist <- res$nn.dist
  } else {
    distmat <- as.matrix(dist(coords))
    diag(distmat) <- Inf
    nn_idx <- t(apply(distmat, 1, function(d) order(d)[seq_len(k)]))
    nn_dist <- t(apply(distmat, 1, function(d) sort(d)[seq_len(k)]))
  }

  # weights ------------------------------------------------------------------
  if (weight == "pygsp") {
    # PyGSP NNGraph-style similarity weights: exp(-d^2 / sigma)
    sigma <- sigma %||% mean(nn_dist[nn_dist > 0], na.rm = TRUE)
    sigma <- if (is.na(sigma) || sigma == 0) 1 else sigma
    w <- exp(-(nn_dist^2) / sigma)
  } else if (weight == "heat") {
    sigma <- sigma %||% stats::median(nn_dist[nn_dist > 0], na.rm = TRUE)
    sigma <- if (is.na(sigma) || sigma == 0) 1 else sigma
    w <- exp(-(nn_dist^2) / (2 * sigma^2))
  } else if (weight == "distance") {
    w <- 1 / (nn_dist + eps)
  } else {
    w <- matrix(1, nrow = nrow(nn_dist), ncol = ncol(nn_dist))
  }

  i <- rep(seq_len(n), each = k)
  j <- as.vector(nn_idx)
  x <- as.vector(w)

  adj <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n))

  if (sym == "union") {
    adj <- (adj + Matrix::t(adj) + abs(adj - Matrix::t(adj))) / 2
  } else if (sym == "mutual") {
    adj <- (adj + Matrix::t(adj) - abs(adj - Matrix::t(adj))) / 2
  } else {
    # average
    adj <- (adj + Matrix::t(adj)) / 2
  }
  adj <- Matrix::drop0(adj)

  new_graph(adj, coords = coords, normalized = normalized, directed = FALSE)
}

#' @rdname graph_knn
#' @export
knn_graph <- graph_knn

#' Stochastic Block Model graph
#' @param block_sizes integer vector of community sizes
#' @param P matrix of edge probabilities (symmetric)
#' @param normalized compute normalized Laplacian
#' @param seed optional RNG seed for reproducibility
#' @return a gsp_graph object
#' @export
graph_sbm <- function(block_sizes, P, normalized = FALSE, seed = NULL) {
  block_sizes <- as.integer(block_sizes)
  if (!is.matrix(P) || nrow(P) != ncol(P) || nrow(P) != length(block_sizes)) {
    stop("P must be a square matrix with dimension equal to length(block_sizes)")
  }
  if (!isTRUE(all.equal(P, t(P)))) {
    stop("P must be symmetric")
  }
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }

  n <- sum(block_sizes)
  offsets <- cumsum(c(0L, head(block_sizes, -1L)))

  edge_i <- list()
  edge_j <- list()
  idx <- 1L
  append_edges <- function(i_local, j_local, offset_i, offset_j) {
    edge_i[[idx]] <<- i_local + offset_i
    edge_j[[idx]] <<- j_local + offset_j
    idx <<- idx + 1L
  }

  for (a in seq_along(block_sizes)) {
    for (b in a:length(block_sizes)) {
      prob <- P[a, b]
      if (prob <= 0) next

      ka <- block_sizes[a]
      kb <- block_sizes[b]
      off_a <- offsets[a]
      off_b <- offsets[b]

      if (a == b) {
        if (ka < 2) next
        if (prob >= 1) {
          comb <- utils::combn(ka, 2)
          i_local <- comb[1, ]
          j_local <- comb[2, ]
        } else {
          M <- Matrix::rsparsematrix(ka, ka, density = prob, symmetric = TRUE,
                                     rand.x = function(nnz) rep(1, nnz))
          if (ka > 0) diag(M) <- 0
          M <- Matrix::triu(M, 1)
          if (length(M@x) == 0) next
          s <- Matrix::summary(M)
          i_local <- s$i
          j_local <- s$j
        }
        append_edges(i_local, j_local, off_a, off_a)
      } else {
        if (prob >= 1) {
          i_local <- rep.int(seq_len(ka), times = kb)
          j_local <- rep(seq_len(kb), each = ka)
        } else {
          M <- Matrix::rsparsematrix(ka, kb, density = prob,
                                     rand.x = function(nnz) rep(1, nnz))
          if (length(M@x) == 0) next
          s <- Matrix::summary(M)
          i_local <- s$i
          j_local <- s$j
        }
        append_edges(i_local, j_local, off_a, off_b)
      }
    }
  }

  if (length(edge_i) == 0) {
    adj <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
  } else {
    i_idx <- unlist(edge_i, use.names = FALSE)
    j_idx <- unlist(edge_j, use.names = FALSE)
    adj <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1, dims = c(n, n))
  }
  adj <- Matrix::forceSymmetric(adj, uplo = "U")
  new_graph(adj, normalized = normalized, directed = FALSE)
}

#' Erdos-Renyi G(n,p) random graph
#'
#' @param n number of nodes
#' @param p edge probability
#' @param directed logical; if TRUE, edges oriented independently
#' @param seed optional RNG seed
#' @param normalized logical; use normalized Laplacian
#' @export
graph_erdos_renyi <- function(n, p, directed = FALSE, seed = NULL, normalized = FALSE) {
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  if (!directed) {
    comb <- utils::combn(n, 2)
    draw <- rbinom(ncol(comb), 1, p)
    keep <- which(draw == 1)
    i <- c(comb[1, keep], comb[2, keep])
    j <- c(comb[2, keep], comb[1, keep])
  } else {
    m <- n * (n - 1)
    draw <- rbinom(m, 1, p)
    pos <- which(draw == 1) - 1
    i <- pos %% n + 1
    j <- pos %/% n + 1
    keep <- i != j
    i <- i[keep]; j <- j[keep]
  }
  adj <- Matrix::sparseMatrix(i = i, j = j, x = 1, dims = c(n, n))
  new_graph(adj, normalized = normalized, directed = directed)
}

#' Complete graph
#' @param n number of nodes
#' @param directed logical; directed edges (no self-loops)
#' @param normalized logical; normalized Laplacian
#' @export
graph_complete <- function(n, directed = FALSE, normalized = FALSE) {
  adj <- Matrix::Matrix(1, n, n, sparse = TRUE)
  diag(adj) <- 0
  if (!directed) adj <- Matrix::forceSymmetric(adj, uplo = "U")
  new_graph(adj, normalized = normalized, directed = directed)
}

#' Random regular graph (configuration model)
#' @param n number of nodes
#' @param degree degree for each node (n*degree must be even)
#' @param seed optional RNG seed
#' @param normalized logical; normalized Laplacian
#' @export
graph_random_regular <- function(n, degree, seed = NULL, normalized = FALSE) {
  if (degree * n %% 2 != 0) stop("n*degree must be even")
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  stubs <- rep(seq_len(n), each = degree)
  stubs <- sample(stubs)
  i <- stubs[seq(1, length(stubs), by = 2)]
  j <- stubs[seq(2, length(stubs), by = 2)]
  adj <- Matrix::sparseMatrix(i = c(i, j), j = c(j, i), x = 1, dims = c(n, n))
  adj <- Matrix::drop0(adj)
  new_graph(adj, normalized = normalized, directed = FALSE)
}

#' Barabasi-Albert preferential attachment graph
#' @param n number of nodes
#' @param m edges per new node
#' @param seed optional RNG seed
#' @param normalized logical; normalized Laplacian
#' @export
graph_barabasi_albert <- function(n, m = 2, seed = NULL, normalized = FALSE) {
  if (m < 1 || m >= n) stop("m must be between 1 and n-1")
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
    set.seed(seed)
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = .GlobalEnv) else assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  deg <- rep(0L, n)
  edges_i <- integer(0); edges_j <- integer(0)
  # start with m+1 clique
  for (i in 1:(m + 1)) {
    for (j in (i + 1):(m + 1)) {
      edges_i <- c(edges_i, i, j)
      edges_j <- c(edges_j, j, i)
      deg[i] <- deg[i] + 1
      deg[j] <- deg[j] + 1
    }
  }
  for (v in (m + 2):n) {
    probs <- deg[1:(v - 1)]
    probs <- probs / sum(probs)
    targets <- sample.int(v - 1, m, prob = probs, replace = FALSE)
    for (tgt in targets) {
      edges_i <- c(edges_i, v, tgt)
      edges_j <- c(edges_j, tgt, v)
      deg[v] <- deg[v] + 1
      deg[tgt] <- deg[tgt] + 1
    }
  }
  adj <- Matrix::sparseMatrix(i = edges_i, j = edges_j, x = 1, dims = c(n, n))
  adj <- Matrix::forceSymmetric(adj, uplo = "U")
  new_graph(adj, normalized = normalized, directed = FALSE)
}

# Caching helpers -----------------------------------------------------------

#' Clear cached quantities stored in a graph.
#' @param g a gsp_graph object
#' @return the graph with cleared cache
#' @export
graph_cache_clear <- function(g) {
  if (is.environment(g$cache)) {
    rm(list = ls(envir = g$cache, all.names = TRUE), envir = g$cache)
  } else {
    g$cache <- list()
  }
  g
}

#' Get (and cache) degree vector.
#' @param g a gsp_graph object
#' @return numeric vector of node degrees
#' @export
graph_degree <- function(g) {
  if (!is.null(g$degree)) return(g$degree)
  if (!is.null(g$cache$degree)) return(g$cache$degree)

  deg <- degree_cpp(g$adjacency)
  g$cache$degree <- deg
  deg
}

#' Get (and cache) Laplacian.
#' @param g a gsp_graph object
#' @param normalized use normalized Laplacian
#' @return sparse Laplacian matrix
#' @export
graph_laplacian <- function(g, normalized = g$normalized) {
  key <- paste0("lap_", normalized)
  if (identical(normalized, g$normalized) && !is.null(g$laplacian)) {
    return(g$laplacian)
  }
  if (!is.null(g$cache[[key]])) return(g$cache[[key]])

  L <- laplacian_cpp(g$adjacency, normalized = normalized, directed = g$directed)
  g$cache[[key]] <- L
  L
}

#' Rescale edge weights by a factor (returns a new graph)
#' @param g a gsp_graph object
#' @param factor multiplicative factor for edge weights
#' @return a new gsp_graph with rescaled weights
#' @export
graph_rescale <- function(g, factor) {
  adj <- g$adjacency * factor
  g2 <- new_graph(adj, coords = g$coords, normalized = g$normalized, directed = g$directed)
  g2
}

#' Sparsify small weights (in-place threshold)
#' @param g a gsp_graph object
#' @param tol values with abs < tol are dropped
#' @return a new gsp_graph with small weights removed
#' @export
graph_sparsify <- function(g, tol = 1e-12) {
  adj <- g$adjacency
  adj@x[abs(adj@x) < tol] <- 0
  adj <- Matrix::drop0(adj)
  if (!g$directed) {
    adj <- Matrix::forceSymmetric(adj, uplo = "U")
  }
  g2 <- new_graph(adj, coords = g$coords, normalized = g$normalized, directed = g$directed)
  g2
}

#' Print method for gsp_graph
#' @param x a gsp_graph object
#' @param ... additional arguments (ignored)
#' @return invisibly returns the graph
#' @method print gsp_graph
#' @export
print.gsp_graph <- function(x, ...) {
  edges <- length(x$adjacency@x) / 2
  cat(sprintf("<gsp_graph> n = %d, edges = %d, normalized = %s\n",
              x$n, edges, x$normalized))
  if (!is.null(x$coords)) {
    cat(sprintf(" coords: %d x %d\n", nrow(x$coords), ncol(x$coords)))
  }
  invisible(x)
}

# Directional utilities ------------------------------------------------------

#' Oriented incidence matrix of a graph
#'
#' Returns an E x N sparse incidence matrix with a fixed orientation
#' (i<j for undirected graphs, edge direction for directed graphs).
#' Attributes `edge_weights`, `edge_i`, `edge_j` store edge data.
#' @param g gsp_graph
#' @return sparse matrix of size (n_edges x n_nodes)
#' @export
graph_incidence <- function(g) {
  if (!inherits(g, "gsp_graph")) stop("'g' must be a gsp_graph")
  A <- methods::as(g$adjacency, "TsparseMatrix")
  if (g$directed) {
    idx_i <- A@i + 1L
    idx_j <- A@j + 1L
    w <- A@x
  } else {
    keep <- A@i < A@j
    idx_i <- A@i[keep] + 1L
    idx_j <- A@j[keep] + 1L
    w <- A@x[keep]
  }
  m <- length(idx_i)
  B <- Matrix::sparseMatrix(
    i = rep(seq_len(m), 2),
    j = c(idx_i, idx_j),
    x = c(rep(1, m), rep(-1, m)),
    dims = c(m, g$n)
  )
  attr(B, "edge_weights") <- w
  attr(B, "edge_i") <- idx_i
  attr(B, "edge_j") <- idx_j
  B
}

#' Weighted gradient operator on a graph
#'
#' Gradient matrix D such that Dx gives weighted edge differences (sqrt(w)*(x_i - x_j)).
#' @param g gsp_graph
#' @return sparse matrix (n_edges x n_nodes)
#' @export
graph_gradient <- function(g) {
  A <- methods::as(g$adjacency, "TsparseMatrix")
  if (g$directed) {
    idx_i <- A@i + 1L
    idx_j <- A@j + 1L
    w <- A@x
  } else {
    keep <- A@i < A@j
    idx_i <- A@i[keep] + 1L
    idx_j <- A@j[keep] + 1L
    w <- A@x[keep]
  }
  m <- length(idx_i)
  if (m == 0) {
    stop("graph_gradient(): graph has no edges")
  }
  sqrtw <- sqrt(w)
  Matrix::sparseMatrix(
    i = rep(seq_len(m), 2),
    j = c(idx_i, idx_j),
    x = c(sqrtw, -sqrtw),
    dims = c(m, g$n)
  )
}

.unit_vectors_axes <- function(d) {
  diag(d)
}

.unit_vectors_circle <- function(n_dir) {
  theta <- seq(0, 2 * pi, length.out = n_dir + 1)[- (n_dir + 1)]
  cbind(cos(theta), sin(theta))
}

.unit_vectors_fibonacci <- function(n_dir) {
  # Evenly spread points on sphere (rows)
  if (n_dir == 1) return(matrix(c(1, 0, 0), nrow = 1))
  i <- seq_len(n_dir) - 1
  phi <- (1 + sqrt(5)) / 2
  theta <- 2 * pi * i / phi
  y <- 1 - 2 * (i + 0.5) / n_dir
  r <- sqrt(pmax(0, 1 - y^2))
  x <- r * cos(theta)
  z <- r * sin(theta)
  cbind(x, y, z)
}

.normalize_rows <- function(M) {
  norms <- sqrt(rowSums(M^2))
  norms[norms == 0] <- 1
  M / norms
}

#' Build a family of orientation-biased graphs
#'
#' Reweights edges according to their alignment with target orientations,
#' returning one graph per orientation. Useful for directional wavelets.
#' @param g base gsp_graph with coords
#' @param n_directions number of canonical directions (used when `orientations` is NULL)
#' @param orientations optional matrix of unit vectors (rows = directions) overriding n_directions
#' @param alpha_isotropic blend factor with isotropic weights (0=fully directional, 1=isotropic)
#' @param p power applied to alignment (larger -> sharper directionality)
#' @param absorb_sign if TRUE, +u and -u are equivalent (abs dot product)
#' @param min_degree minimum degree to enforce (warns if violated)
#' @param min_weight weights below this are dropped
#' @param check_degrees logical; enforce min_degree and warn on isolation
#' @return list of gsp_graph objects, one per orientation
#' @export
graph_steerable_family <- function(g,
                                   n_directions = 3L,
                                   orientations = NULL,
                                   alpha_isotropic = 0.1,
                                   p = 4L,
                                   absorb_sign = TRUE,
                                   min_degree = 3L,
                                   min_weight = 1e-6,
                                   check_degrees = TRUE) {
  if (is.null(g$coords)) stop("graph_steerable_family(): graph must have coords")
  coords <- as.matrix(g$coords)
  d <- ncol(coords)

  U <- orientations
  if (is.null(U)) {
    if (d == 2) {
      U <- .unit_vectors_circle(n_directions)
    } else if (d == 3) {
      if (n_directions <= 3) {
        U <- .unit_vectors_axes(3)[seq_len(n_directions), , drop = FALSE]
      } else {
        U <- .unit_vectors_fibonacci(n_directions)
      }
    } else {
      U <- .unit_vectors_axes(d)
    }
  }
  U <- as.matrix(U)
  if (ncol(U) != d && nrow(U) == d) {
    # allow transposed input
    U <- t(U)
  }
  if (ncol(U) != d) stop("orientation vectors must have same dimension as coords")
  U <- .normalize_rows(U)
  n_dir <- nrow(U)

  A <- methods::as(g$adjacency, "TsparseMatrix")
  if (g$directed) {
    idx_i <- A@i + 1L
    idx_j <- A@j + 1L
    w0 <- A@x
  } else {
    keep <- A@i < A@j
    idx_i <- A@i[keep] + 1L
    idx_j <- A@j[keep] + 1L
    w0 <- A@x[keep]
  }

  vecs <- coords[idx_j, , drop = FALSE] - coords[idx_i, , drop = FALSE]
  dists <- sqrt(rowSums(vecs^2))
  dists[dists == 0] <- 1
  vecs <- vecs / dists

  graphs <- vector("list", n_dir)

  for (d_idx in seq_len(n_dir)) {
    u <- U[d_idx, , drop = TRUE]
    alignment <- as.numeric(vecs %*% u)
    if (absorb_sign) alignment <- abs(alignment)
    modulation <- alignment^p
    new_w <- w0 * ((1 - alpha_isotropic) * modulation + alpha_isotropic)

    # Optional degree rescue
    if (check_degrees && min_degree > 0) {
      deg_guess <- numeric(g$n)
      # approximate degrees before dropping tiny weights
      for (k in seq_along(idx_i)) {
        deg_guess[idx_i[k]] <- deg_guess[idx_i[k]] + new_w[k]
        deg_guess[idx_j[k]] <- deg_guess[idx_j[k]] + new_w[k]
      }
      if (any(deg_guess < min_degree)) {
        beta <- 0.25
        new_w <- pmax(new_w, beta * w0)
      }
    }

    keep_w <- new_w > min_weight
    A_i <- idx_i[keep_w]
    A_j <- idx_j[keep_w]
    w_new <- new_w[keep_w]

    if (g$directed) {
      W_theta <- Matrix::sparseMatrix(i = A_i, j = A_j, x = w_new,
                                      dims = dim(g$adjacency))
    } else {
      W_theta <- Matrix::sparseMatrix(
        i = c(A_i, A_j),
        j = c(A_j, A_i),
        x = c(w_new, w_new),
        dims = dim(g$adjacency)
      )
      W_theta <- Matrix::forceSymmetric(W_theta, uplo = "U")
    }

    if (check_degrees && min_degree > 0) {
      deg <- degree_cpp(W_theta)
      if (any(deg == 0)) {
        warning("graph_steerable_family(): isolated nodes in direction ", d_idx,
                "; consider lowering p or increasing alpha_isotropic")
      }
    }

    g_theta <- new_graph(
      Matrix::drop0(W_theta),
      coords = g$coords,
      normalized = g$normalized,
      directed = g$directed
    )
    g_theta$orientation <- u
    graphs[[d_idx]] <- g_theta
  }

  graphs
}
