# Directional SGWT smoke benchmark (small, non-CRAN)
#
# Purpose: quick sanity on runtime/memory scaling with #directions.
# Keep sizes tiny to avoid CI timeouts; adjust locally for larger graphs.

suppressPackageStartupMessages({
  library(rgsp)
})

bench_one <- function(nx = 30, ny = 30, D = 4, n_scales = 2, wavelet = "mexican_hat",
                      memory_saving = FALSE) {
  cat(sprintf("\nGrid %dx%d, directions=%d, scales=%d, wavelet=%s\n",
              nx, ny, D, n_scales, wavelet))
  g <- graph_grid2d(nx, ny, normalized = FALSE)
  g <- graph_set_coords(g, kind = g$coords)
  x <- matrix(rnorm(g$n * 5), nrow = g$n)  # 5 signals
  scales <- auto_wavelet_scales(g, n_scales = n_scales, lmax = lambda_max(g))

  t_dir <- system.time({
    dsgwt_steer(g, x,
                n_directions = D,
                scales = scales,
                wavelet = wavelet,
                include_lowpass = TRUE,
                K = 20,
                alpha_isotropic = 0.3,
                p = 2,
                memory_saving = memory_saving)
  })
  cat("  elapsed (sec):", round(t_dir["elapsed"], 3), "\n")
}

if (interactive()) {
  options(rgsp.threads = 1L)
  bench_one(20, 20, D = 3)
  bench_one(30, 30, D = 6)
}
