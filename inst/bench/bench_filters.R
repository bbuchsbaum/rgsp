# Benchmark Chebyshev filtering vs PyGSP (optional)
#
# Usage:
#   Rscript inst/bench/bench_filters.R
#
# Requirements:
#   - R packages: bench, reticulate, Matrix
#   - Python with pygsp installed and on PATH/RETICULATE_PYTHON

suppressPackageStartupMessages({
  library(bench)
  library(Matrix)
  library(reticulate)
  library(rgsp)
})

pygsp <- NULL
try({
  pygsp <- import("pygsp", delay_load = TRUE)
}, silent = TRUE)

make_ring <- function(n) graph_ring(n)
make_grid <- function(n) graph_grid2d(n, n, periodic = TRUE)

rand_signal <- function(n, m = 4) matrix(rnorm(n * m), ncol = m)

bench_one <- function(gfun, n, kernel, K = 50, reps = 5) {
  g <- gfun(n)
  x <- rand_signal(g$n, 4)
  lmax <- lambda_max(g)
  b <- mark(
    rgsp = filter_signal(g, x, kernel, K = K, lmax = lmax),
    iterations = reps,
    check = FALSE
  )
  if (!is.null(pygsp)) {
    Gpy <- if (identical(gfun, make_ring)) {
      pygsp$graphs$Ring(n)
    } else {
      pygsp$graphs$Grid2d(n, n, periodic = TRUE)
    }
    fk <- function(l) exp(-0.5 * l)
    fp <- pygsp$filters$Filter(Gpy, fk)
    xpy <- x
    bpy <- system_time({
      res <- fp$filter(xpy)
    })
    b$pygsp <- bpy[["real"]]
  }
  b
}

main <- function() {
  ks <- c(50)
  sizes <- c(128, 256)
  kernel <- kernel_heat(0.5)
  results <- list()
  for (n in sizes) {
    results[[paste0("ring_", n)]] <- bench_one(make_ring, n, kernel, K = ks)
    results[[paste0("grid_", n)]] <- bench_one(make_grid, n, kernel, K = ks)
  }
  print(results)
}

if (identical(sys.nframe(), 0L)) {
  main()
}
