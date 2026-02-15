# Incidence-based isotropic TV benchmark (rgsp)
#
# Usage:
#   Rscript inst/bench/bench_incidence_tv.R

suppressPackageStartupMessages({
  library(bench)
  library(Matrix)
  library(rgsp)
})

run_case <- function(n = 2048, lambda = 1, maxit = 80) {
  g <- graph_ring(n)
  y <- sin(seq(0, 2*pi, length.out = n)) + rnorm(n, sd = 0.2)
  res <- mark(
    tv = tv_denoise(g, y, lambda = lambda, maxit = maxit, tau = 0.5, sigma = 0.5),
    tikh = tikhonov_smooth(g, y, tau = 0.2),
    iterations = 3,
    check = FALSE
  )
  res
}

main <- function() {
  res <- run_case()
  print(res)
  out_path <- Sys.getenv("BENCH_OUT", NA)
  if (!is.na(out_path) && nzchar(out_path)) {
    df <- as.data.frame(res)
    df$name <- rownames(df)
    write.csv(df, out_path, row.names = FALSE)
    message("Wrote benchmark to ", out_path)
  }
}

if (identical(sys.nframe(), 0L)) main()
