# Benchmark Chebyshev filtering vs PyGSP
#
# Usage:
#   Rscript inst/bench/bench_chebyshev.R
#
# Compares rgsp Chebyshev polynomial filtering against PyGSP for various:
# - Graph sizes
# - Filter orders (K)
# - Number of signals
# - Different kernels

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

bench_chebyshev_r <- function(g, signals, kernel, K, lmax, reps = 5) {
  # Clear cache to ensure fair timing
  cheby_cache_clear()

  mark(
    chebyshev = filter_signal(g, signals, kernel, K = K, lmax = lmax),
    iterations = reps,
    check = FALSE
  )
}

bench_chebyshev_vs_pygsp <- function(n, kernel_name = "heat", K = 30,
                                      n_signals = 1, reps = 5) {
  # Create ring graph (simple, well-understood structure)
  g_r <- graph_ring(n)
  lmax_r <- lambda_max(g_r)

  # Generate signals
  set.seed(42)
  signals <- matrix(rnorm(n * n_signals), nrow = n, ncol = n_signals)

  # Select kernel
  kernel_r <- switch(kernel_name,
    heat = kernel_heat(1.0),
    mexican_hat = kernel_mexican_hat(1.0),
    exponential = kernel_exponential(1.0),
    stop("Unknown kernel")
  )

  # R benchmark
  b_r <- bench_chebyshev_r(g_r, signals, kernel_r, K, lmax_r, reps)

  result <- data.frame(
    n = n,
    K = K,
    n_signals = n_signals,
    kernel = kernel_name,
    rgsp_median_ms = as.numeric(median(b_r$median)) * 1000,
    rgsp_mem_alloc = as.numeric(b_r$mem_alloc)
  )

  # PyGSP comparison if available
  if (pygsp_available) {
    g_py <- pygsp$graphs$Ring(as.integer(n))
    g_py$estimate_lmax()

    # Create PyGSP filter
    filt_py <- switch(kernel_name,
      heat = pygsp$filters$Heat(g_py, tau = 1.0),
      mexican_hat = pygsp$filters$MexicanHat(g_py, Nf = 1L),
      exponential = pygsp$filters$Expwin(g_py),
      stop("Unknown kernel for PyGSP")
    )

    signals_py <- np$array(signals)

    # PyGSP benchmark
    b_py <- mark(
      pygsp = filt_py$filter(signals_py, method = "chebyshev", order = as.integer(K)),
      iterations = reps,
      check = FALSE
    )

    result$pygsp_median_ms <- as.numeric(median(b_py$median)) * 1000
    result$speedup <- result$pygsp_median_ms / result$rgsp_median_ms
  }

  result
}

# ==============================================================================
# Scaling Benchmarks
# ==============================================================================

bench_scaling_n <- function(K = 30, n_signals = 4) {
  sizes <- c(64, 128, 256, 512, 1024, 2048)
  results <- lapply(sizes, function(n) {
    cat(sprintf("  n = %d\n", n))
    bench_chebyshev_vs_pygsp(n, "heat", K, n_signals, reps = 5)
  })
  do.call(rbind, results)
}

bench_scaling_K <- function(n = 512, n_signals = 4) {
  orders <- c(10, 20, 30, 50, 100)
  results <- lapply(orders, function(K) {
    cat(sprintf("  K = %d\n", K))
    bench_chebyshev_vs_pygsp(n, "heat", K, n_signals, reps = 5)
  })
  do.call(rbind, results)
}

bench_scaling_signals <- function(n = 512, K = 30) {
  signal_counts <- c(1, 4, 16, 64)
  results <- lapply(signal_counts, function(ns) {
    cat(sprintf("  signals = %d\n", ns))
    bench_chebyshev_vs_pygsp(n, "heat", K, ns, reps = 5)
  })
  do.call(rbind, results)
}

# ==============================================================================
# Main Benchmark Suite
# ==============================================================================

main <- function() {
  cat(strrep("=", 70), "\n")
  cat("Chebyshev Filtering Benchmark Suite\n")
  cat(strrep("=", 70), "\n\n")

  all_results <- list()

  # 1. Scaling with graph size
  cat("1. Scaling with graph size (K=30, signals=4)\n")
  cat(strrep("-", 50), "\n")
  res_n <- bench_scaling_n()
  all_results$scaling_n <- res_n
  print(res_n[, c("n", "K", "rgsp_median_ms",
                  if ("pygsp_median_ms" %in% names(res_n)) "pygsp_median_ms" else NULL,
                  if ("speedup" %in% names(res_n)) "speedup" else NULL)],
        row.names = FALSE)
  cat("\n")

  # 2. Scaling with filter order
  cat("2. Scaling with filter order K (n=512, signals=4)\n")
  cat(strrep("-", 50), "\n")
  res_K <- bench_scaling_K()
  all_results$scaling_K <- res_K
  print(res_K[, c("n", "K", "rgsp_median_ms",
                  if ("pygsp_median_ms" %in% names(res_K)) "pygsp_median_ms" else NULL,
                  if ("speedup" %in% names(res_K)) "speedup" else NULL)],
        row.names = FALSE)
  cat("\n")

  # 3. Scaling with number of signals
  cat("3. Scaling with number of signals (n=512, K=30)\n")
  cat(strrep("-", 50), "\n")
  res_sig <- bench_scaling_signals()
  all_results$scaling_signals <- res_sig
  print(res_sig[, c("n", "n_signals", "rgsp_median_ms",
                    if ("pygsp_median_ms" %in% names(res_sig)) "pygsp_median_ms" else NULL,
                    if ("speedup" %in% names(res_sig)) "speedup" else NULL)],
        row.names = FALSE)
  cat("\n")

  # 4. Different kernels
  cat("4. Different kernels (n=512, K=30, signals=4)\n")
  cat(strrep("-", 50), "\n")
  kernels <- c("heat", "exponential")
  res_kern <- lapply(kernels, function(k) {
    cat(sprintf("  kernel = %s\n", k))
    tryCatch(
      bench_chebyshev_vs_pygsp(512, k, 30, 4, reps = 5),
      error = function(e) NULL
    )
  })
  res_kern <- do.call(rbind, Filter(Negate(is.null), res_kern))
  if (!is.null(res_kern) && nrow(res_kern) > 0) {
    all_results$kernels <- res_kern
    print(res_kern[, c("kernel", "rgsp_median_ms",
                       if ("pygsp_median_ms" %in% names(res_kern)) "pygsp_median_ms" else NULL,
                       if ("speedup" %in% names(res_kern)) "speedup" else NULL)],
          row.names = FALSE)
  }

  # Summary
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("Summary\n")
  cat(strrep("=", 70), "\n")

  if (pygsp_available && "speedup" %in% names(res_n)) {
    avg_speedup <- mean(c(res_n$speedup, res_K$speedup, res_sig$speedup), na.rm = TRUE)
    cat(sprintf("Average speedup vs PyGSP: %.2fx\n", avg_speedup))
  }

  # Save results
  outfile <- "inst/bench/results_chebyshev.rds"
  if (file.exists(dirname(outfile))) {
    saveRDS(all_results, outfile)
    cat("Results saved to:", outfile, "\n")
  }

  invisible(all_results)
}

if (identical(environment(), globalenv())) {
  main()
}
