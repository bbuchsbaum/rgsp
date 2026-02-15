# Benchmark GFT/IGFT operations vs PyGSP
#
# Usage:
#   Rscript inst/bench/bench_gft.R
#
# Requirements:
#   - R packages: bench, reticulate, Matrix
#   - Python with pygsp installed (optional, for comparison)

suppressPackageStartupMessages({
  library(bench)
  library(Matrix)
  library(rgsp)
})

# Try to load PyGSP for comparison
pygsp_available <- FALSE
pygsp <- NULL
np <- NULL

if (requireNamespace("reticulate", quietly = TRUE)) {
  tryCatch({
    pygsp <- reticulate::import("pygsp", delay_load = FALSE)
    np <- reticulate::import("numpy", delay_load = FALSE)
    pygsp_available <- TRUE
    message("PyGSP available for comparison benchmarks")
  }, error = function(e) {
    message("PyGSP not available; running R-only benchmarks")
  })
}

# ==============================================================================
# Benchmark Functions
# ==============================================================================

bench_gft <- function(g, signals, reps = 5) {
  # Ensure eigenpairs are cached
  graph_eigenpairs(g)

  mark(
    gft = gft(g, signals),
    iterations = reps,
    check = FALSE
  )
}

bench_igft <- function(g, coeffs, reps = 5) {
  mark(
    igft = igft(g, coeffs),
    iterations = reps,
    check = FALSE
  )
}

bench_gft_roundtrip <- function(g, signals, reps = 5) {
  # Ensure eigenpairs are cached
  graph_eigenpairs(g)

  mark(
    roundtrip = {
      coeffs <- gft(g, signals)
      igft(g, coeffs)
    },
    iterations = reps,
    check = FALSE
  )
}

bench_gft_vs_pygsp <- function(n, graph_type = "ring", n_signals = 1, reps = 5) {
  # Create graph
  g_r <- switch(graph_type,
    ring = graph_ring(n),
    path = graph_path(n),
    grid = graph_grid2d(sqrt(n), sqrt(n)),
    stop("Unknown graph type")
  )

  # Generate signals
  set.seed(42)
  signals <- matrix(rnorm(n * n_signals), nrow = n, ncol = n_signals)

  # R benchmark (with eigenpair computation)
  t_eigen_r <- system.time({
    graph_eigenpairs(g_r)
  })["elapsed"]

  # R GFT benchmark (eigenpairs cached)
  b_r <- bench_gft(g_r, signals, reps)

  result <- data.frame(
    graph = graph_type,
    n = n,
    n_signals = n_signals,
    rgsp_eigen_s = t_eigen_r,
    rgsp_gft_median_ms = as.numeric(median(b_r$median)) * 1000
  )

  # PyGSP comparison if available
if (pygsp_available) {
    g_py <- switch(graph_type,
      ring = pygsp$graphs$Ring(as.integer(n)),
      path = pygsp$graphs$Path(as.integer(n)),
      grid = pygsp$graphs$Grid2d(as.integer(sqrt(n)), as.integer(sqrt(n))),
      stop("Unknown graph type")
    )

    # PyGSP eigenpair computation
    t_eigen_py <- system.time({
      g_py$compute_fourier_basis()
    })["elapsed"]

    # PyGSP GFT
    signals_py <- np$array(signals)
    t_gft_py <- bench::mark(
      pygsp_gft = g_py$gft(signals_py),
      iterations = reps,
      check = FALSE
    )

    result$pygsp_eigen_s <- t_eigen_py
    result$pygsp_gft_median_ms <- as.numeric(median(t_gft_py$median)) * 1000
    result$speedup_eigen <- t_eigen_py / t_eigen_r
    result$speedup_gft <- result$pygsp_gft_median_ms / result$rgsp_gft_median_ms
  }

  result
}

# ==============================================================================
# Main Benchmark Suite
# ==============================================================================

main <- function() {
  cat("=" %s+% strrep("=", 70) %s+% "\n")
  cat("GFT/IGFT Benchmark Suite\n")
  cat("=" %s+% strrep("=", 70) %s+% "\n\n")

  # Graph sizes to test
  sizes <- c(64, 128, 256, 512, 1024)
  graph_types <- c("ring", "path")
  signal_counts <- c(1, 4, 16)

  results <- list()

  for (gt in graph_types) {
    for (n in sizes) {
      for (ns in signal_counts) {
        cat(sprintf("Benchmarking: %s graph, n=%d, signals=%d\n", gt, n, ns))
        tryCatch({
          res <- bench_gft_vs_pygsp(n, gt, ns, reps = 5)
          results[[length(results) + 1]] <- res
        }, error = function(e) {
          cat("  Error:", conditionMessage(e), "\n")
        })
      }
    }
  }

  # Combine results
  df <- do.call(rbind, results)

  cat("\n")
  cat("=" %s+% strrep("=", 70) %s+% "\n")
  cat("Results Summary\n")
  cat("=" %s+% strrep("=", 70) %s+% "\n\n")

  print(df, row.names = FALSE)

  # Save results
  outfile <- "inst/bench/results_gft.csv"
  if (file.exists(dirname(outfile))) {
    write.csv(df, outfile, row.names = FALSE)
    cat("\nResults saved to:", outfile, "\n")
  }

  invisible(df)
}

# String concatenation helper
`%s+%` <- function(a, b) paste0(a, b)

if (identical(environment(), globalenv())) {
  main()
}
