# Run all rgsp benchmarks
#
# Usage:
#   Rscript inst/bench/bench_all.R [--quick]
#
# Options:
#   --quick  Run quick benchmarks with smaller sizes (for CI smoke tests)

suppressPackageStartupMessages({
  library(rgsp)
})

args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args

cat(strrep("=", 70), "\n")
cat("rgsp Comprehensive Benchmark Suite\n")
cat(sprintf("Mode: %s\n", if (quick_mode) "QUICK (CI smoke test)" else "FULL"))
cat(sprintf("Date: %s\n", Sys.time()))
cat(strrep("=", 70), "\n\n")

# Check for PyGSP
pygsp_available <- FALSE
if (requireNamespace("reticulate", quietly = TRUE)) {
  pygsp_available <- tryCatch({
    reticulate::import("pygsp")
    TRUE
  }, error = function(e) FALSE)
}
cat(sprintf("PyGSP available: %s\n\n", pygsp_available))

# ==============================================================================
# Quick Benchmarks (for CI)
# ==============================================================================

run_quick_benchmarks <- function() {
  cat("Running quick benchmarks...\n\n")

  results <- list()

  # 1. GFT/IGFT smoke test
  cat("1. GFT/IGFT...\n")
  g <- graph_ring(128)
  x <- rnorm(128)

  t_eigen <- system.time(graph_eigenpairs(g))["elapsed"]
  t_gft <- system.time(for (i in 1:10) gft(g, x))["elapsed"] / 10
  t_igft <- system.time({
    coeffs <- gft(g, x)
    for (i in 1:10) igft(g, coeffs)
  })["elapsed"] / 10

  results$gft <- data.frame(
    operation = c("eigenpairs", "gft", "igft"),
    n = 128,
    time_ms = c(t_eigen, t_gft, t_igft) * 1000
  )
  print(results$gft)
  cat("\n")

  # 2. Chebyshev filtering smoke test
  cat("2. Chebyshev filtering...\n")
  cheby_cache_clear()
  t_filter <- system.time({
    for (i in 1:10) filter_signal(g, x, kernel_heat(1), K = 30)
  })["elapsed"] / 10

  results$chebyshev <- data.frame(
    operation = "filter_signal",
    n = 128,
    K = 30,
    time_ms = t_filter * 1000
  )
  print(results$chebyshev)
  cat("\n")

  # 3. SGWT smoke test
  cat("3. SGWT...\n")
  cheby_cache_clear()
  t_sgwt <- system.time({
    for (i in 1:5) sgwt(g, x, scales = 3)
  })["elapsed"] / 5

  W <- sgwt(g, x, scales = 3)
  t_isgwt <- system.time({
    for (i in 1:5) isgwt(W)
  })["elapsed"] / 5

  results$sgwt <- data.frame(
    operation = c("sgwt", "isgwt"),
    n = 128,
    scales = 3,
    time_ms = c(t_sgwt, t_isgwt) * 1000
  )
  print(results$sgwt)
  cat("\n")

  # 4. Graph construction smoke test
  cat("4. Graph construction...\n")
  t_ring <- system.time(for (i in 1:100) graph_ring(100))["elapsed"] / 100
  t_grid <- system.time(for (i in 1:100) graph_grid2d(10, 10))["elapsed"] / 100
  t_sensor <- system.time(for (i in 1:10) graph_sensor(100, seed = i))["elapsed"] / 10

  results$graphs <- data.frame(
    graph = c("ring", "grid2d", "sensor"),
    n = 100,
    time_ms = c(t_ring, t_grid, t_sensor) * 1000
  )
  print(results$graphs)
  cat("\n")

  # Summary
  cat(strrep("=", 70), "\n")
  cat("Quick Benchmark Summary: PASS\n")
  cat(strrep("=", 70), "\n")

  invisible(results)
}

# ==============================================================================
# Full Benchmarks
# ==============================================================================

run_full_benchmarks <- function() {
  cat("Running full benchmarks...\n\n")

  bench_dir <- system.file("bench", package = "rgsp")
  if (bench_dir == "") bench_dir <- "inst/bench"

  results <- list()

  # Run individual benchmark scripts
  scripts <- c("bench_gft.R", "bench_chebyshev.R", "bench_wavelets.R")

  for (script in scripts) {
    script_path <- file.path(bench_dir, script)
    if (file.exists(script_path)) {
      cat(strrep("-", 70), "\n")
      cat(sprintf("Running %s\n", script))
      cat(strrep("-", 70), "\n")
      tryCatch({
        source(script_path, local = TRUE)
        results[[script]] <- "completed"
      }, error = function(e) {
        cat("Error:", conditionMessage(e), "\n")
        results[[script]] <- paste("error:", conditionMessage(e))
      })
      cat("\n")
    } else {
      cat(sprintf("Script not found: %s\n", script_path))
    }
  }

  invisible(results)
}

# ==============================================================================
# Main
# ==============================================================================

if (quick_mode) {
  run_quick_benchmarks()
} else {
  run_full_benchmarks()
}

cat("\nBenchmark run completed.\n")
