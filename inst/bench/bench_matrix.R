# Benchmark matrix of configurations for Chebyshev filtering
#
# Usage:
#   Rscript inst/bench/bench_matrix.R
#
# Outputs timing summary to console; extend to write CSV if desired.

suppressPackageStartupMessages({
  library(bench)
  library(Matrix)
  library(rgsp)
})

run_case <- function(graph_fun, n, K, t) {
  g <- graph_fun(n)
  x <- matrix(rnorm(g$n * 4), ncol = 4)
  lmax <- lambda_max(g)
  kernel <- kernel_heat(t)
  mark(
    rgsp = filter_signal(g, x, kernel, K = K, lmax = lmax),
    iterations = 3,
    check = FALSE
  )
}

main <- function() {
  cases <- expand.grid(
    graph = c("ring", "grid"),
    n = c(128, 256, 512),
    K = c(30, 50, 80),
    t = c(0.1, 0.5)
  )
  res <- vector("list", nrow(cases))
  for (i in seq_len(nrow(cases))) {
    row <- cases[i, ]
    gfun <- if (row$graph == "ring") function(n) graph_ring(n) else function(n) graph_grid2d(sqrt(n), sqrt(n), periodic = TRUE)
    res[[i]] <- cbind(row, summary = list(run_case(gfun, row$n, row$K, row$t)))
  }
  print(res)
}

if (identical(sys.nframe(), 0L)) main()
