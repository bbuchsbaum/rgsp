# Compare rgsp TV denoise with PyGSP heat/TV (if available)
#
# Usage:
#   Rscript inst/bench/bench_tv_parity_py.R

suppressPackageStartupMessages({
  library(rgsp)
  library(Matrix)
  library(reticulate)
})

py_available <- FALSE
pygsp <- NULL
try({
  pygsp <- import("pygsp", delay_load = TRUE)
  py_available <- TRUE
}, silent = TRUE)

if (!py_available) {
  message("PyGSP not available; skipping parity check.")
  quit(status = 0)
}

set.seed(1)
n <- 256
g <- graph_ring(n)
y <- sin(seq(0, 2 * pi, length.out = n)) + rnorm(n, sd = 0.1)

# rgsp TV
x_r <- tv_denoise(g, y, lambda = 1, maxit = 150)

# PyGSP: approximate by heat filter as proxy (PyGSP lacks direct TV prox API)
Gpy <- pygsp$graphs$Ring(n)
filt <- pygsp$filters$Heat(Gpy, tau = 1.0)
x_py <- filt$filter(y)

err <- sqrt(mean((x_r - x_py)^2))
cat(sprintf("Parity RMSE (rgsp TV vs PyGSP heat proxy): %.4f\n", err))

out_path <- Sys.getenv("BENCH_OUT", NA)
if (!is.na(out_path) && nzchar(out_path)) {
  write.csv(data.frame(rmse = err), out_path, row.names = FALSE)
  message("Wrote parity result to ", out_path)
}
