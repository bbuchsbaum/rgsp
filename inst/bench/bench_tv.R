# Benchmark TV denoising vs Tikhonov smoothing

suppressPackageStartupMessages({
  library(bench)
  library(rgsp)
})

bench_tv <- function(n = 1024, lambda_tv = 1, tau = 0.5, sigma = 0.5) {
  g <- graph_ring(n)
  y <- rnorm(n)
  lmax <- lambda_max(g)
  # Tikhonov baseline
  b_tikh <- mark(
    tikhonov = tikhonov_smooth(g, y, tau = 0.1),
    iterations = 3, check = FALSE
  )
  # TV denoise (short maxit for smoke)
  b_tv <- mark(
    tv = tv_denoise(g, y, lambda = lambda_tv, maxit = 80, tau = tau, sigma = sigma),
    iterations = 3, check = FALSE
  )
  list(tikhonov = b_tikh, tv = b_tv)
}

main <- function() {
  res <- bench_tv()
  print(res)
}

if (identical(sys.nframe(), 0L)) main()
