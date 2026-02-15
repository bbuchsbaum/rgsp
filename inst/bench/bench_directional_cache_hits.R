# Benchmark cache hits/misses for directional SGWT
# Not run by default; intended for local inspection.

suppressPackageStartupMessages({
  library(rgsp)
})

cache_counts <- function() length(ls(envir = get(".cheby_cache", envir = asNamespace("rgsp"))))

run_cache_bench <- function(nx = 20, ny = 20, D1 = 3, D2 = 8, n_scales = 2) {
  g <- graph_grid2d(nx, ny, normalized = FALSE)
  g <- graph_set_coords(g, kind = g$coords)
  x <- matrix(rnorm(g$n * 3), nrow = g$n)  # 3 signals
  scales <- auto_wavelet_scales(g, n_scales = n_scales, lmax = lambda_max(g))

  cheby_cache_clear()
  cat("cache before:", cache_counts(), "\n")
  dsgwt_steer(g, x, n_directions = D1, scales = scales,
              wavelet = "mexican_hat", include_lowpass = TRUE, K = 12,
              alpha_isotropic = 0.2, p = 2)
  c1 <- cache_counts()
  cat("cache after D1:", c1, "\n")

  dsgwt_steer(g, x, n_directions = D2, scales = scales,
              wavelet = "mexican_hat", include_lowpass = TRUE, K = 12,
              alpha_isotropic = 0.2, p = 2)
  c2 <- cache_counts()
  cat("cache after D2:", c2, " (should match D1)\n")
}

if (interactive()) {
  options(rgsp.threads = 1L)
  run_cache_bench()
}
