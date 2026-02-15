# Thread scaling benchmark for Chebyshev filtering
#
# Usage:
#   Rscript inst/bench/bench_threads.R

suppressPackageStartupMessages({
  library(bench)
  library(Matrix)
  library(rgsp)
})

set_threads <- function(n) {
  # use MKL/BLAS threads via RhpcBLASctl if available
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(n)
    RhpcBLASctl::omp_set_num_threads(n)
  }
  options(rgsp.threads = n)
}

bench_threads <- function(n_threads, g, x, kernel, K, lmax) {
  set_threads(n_threads)
  mark(
    filter_signal(g, x, kernel, K = K, lmax = lmax),
    iterations = 3,
    check = FALSE
  )
}

main <- function() {
  g <- graph_ring(1024)
  x <- matrix(rnorm(g$n * 4), ncol = 4)
  lmax <- lambda_max(g)
  kernel <- kernel_heat(0.2)
  thread_counts <- c(1, 2, 4, 8)
  res <- lapply(thread_counts, function(nt) {
    data.frame(
      threads = nt,
      bench_threads(nt, g, x, kernel, K = 50, lmax = lmax)
    )
  })
  print(res)
}

if (identical(sys.nframe(), 0L)) main()
