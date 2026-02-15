# Benchmark rgsp vs PyGSP for key ops (GFT/IGFT, heat filtering)
#
# Usage:
#   Rscript inst/bench/bench_pygsp_compare.R
#   BENCH_OUT=bench_pygsp.csv Rscript inst/bench/bench_pygsp_compare.R
#
# Requires: reticulate + PyGSP installed; skips gracefully otherwise.

suppressPackageStartupMessages({
  library(rgsp)
  library(bench)
  library(reticulate)
})

py_available <- requireNamespace("reticulate", quietly = TRUE) &&
  reticulate::py_available(initialize = TRUE) &&
  reticulate::py_module_available("pygsp")

if (!py_available) {
  message("PyGSP not available; skipping.")
  quit(status = 0)
}

py <- reticulate::import("pygsp")
np <- reticulate::import("numpy")

set.seed(123)
n_vec <- c(1024, 4096, 8192)   # sizes for heat filtering sweep
K_vec <- c(20, 40, 80)
n_signals <- 8

# Helpers -----------------------------------------------------------------

bench_op <- function(fun, warmup = 2, reps = 8) {
  # warm-up
  for (i in seq_len(warmup)) fun()
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    gc()
    times[i] <- system.time(fun())[["elapsed"]]
  }
  1000 * median(times)  # ms
}

rows <- list()
for (n in n_vec) {
  for (K in K_vec) {
    X_r <- matrix(rnorm(n * n_signals), nrow = n)
    X_py <- np$array(X_r)

    # rgsp objects (heat filter)
    g_r <- graph_ring(n)
    lmax_r <- 4  # analytic for ring

    # PyGSP objects (heat filter)
    g_py <- py$graphs$Ring(as.integer(n))
    g_py$estimate_lmax()
    heat_py <- py$filters$Heat(g_py, 1.0)

    t_rg <- bench_op(function() {
      filter_signal(g_r, X_r, kernel_heat(1.0), K = K, lmax = lmax_r)
    })
    t_py <- bench_op(function() {
      heat_py$filter(X_py)
    })

    rows[[length(rows) + 1]] <- data.frame(
      n = n, K = K, impl = "rgsp", median_ms = t_rg
    )
    rows[[length(rows) + 1]] <- data.frame(
      n = n, K = K, impl = "pygsp", median_ms = t_py
    )
  }
}

df <- do.call(rbind, rows)

print(df)

print(df)

out <- Sys.getenv("BENCH_OUT", NA_character_)
if (!is.na(out) && nzchar(out)) {
  write.csv(df, out, row.names = FALSE)
  message("Wrote results to ", out)
}
