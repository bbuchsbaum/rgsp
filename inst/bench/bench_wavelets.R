# Benchmark Wavelet Transforms (SGWT) vs PyGSP
#
# Usage:
#   Rscript inst/bench/bench_wavelets.R
#
# Compares rgsp spectral graph wavelet transform against PyGSP

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

bench_sgwt_r <- function(g, signal, scales = 4, reps = 5) {
  # Clear caches
  cheby_cache_clear()

  mark(
    sgwt = sgwt(g, signal, scales = scales),
    iterations = reps,
    check = FALSE
  )
}

bench_isgwt_r <- function(W, reps = 5) {
  mark(
    isgwt = isgwt(W),
    iterations = reps,
    check = FALSE
  )
}

bench_wavelet_banks <- function(n = 256, reps = 5) {
  g <- graph_ring(n)
  signal <- rnorm(n)

  results <- list()

  # Mexican hat wavelets
  cat("  Mexican hat bank...\n")
  cheby_cache_clear()
  b_mh <- mark(
    mexican_hat = {
      bank <- mexican_hat_wavelet_bank(g, n_scales = 4)
      filter_bank_apply(g, signal, bank)
    },
    iterations = reps,
    check = FALSE
  )
  results$mexican_hat <- data.frame(
    wavelet = "mexican_hat",
    median_ms = as.numeric(median(b_mh$median)) * 1000
  )

  # Meyer wavelets
  cat("  Meyer bank...\n")
  cheby_cache_clear()
  b_meyer <- mark(
    meyer = {
      bank <- meyer_wavelet_bank(g, n_scales = 4)
      filter_bank_apply(g, signal, bank)
    },
    iterations = reps,
    check = FALSE
  )
  results$meyer <- data.frame(
    wavelet = "meyer",
    median_ms = as.numeric(median(b_meyer$median)) * 1000
  )

  # Heat wavelets
  cat("  Heat bank...\n")
  cheby_cache_clear()
  b_heat <- mark(
    heat = {
      bank <- wavelet_heat_bank(g, n_scales = 4)
      filter_bank_apply(g, signal, bank)
    },
    iterations = reps,
    check = FALSE
  )
  results$heat <- data.frame(
    wavelet = "heat",
    median_ms = as.numeric(median(b_heat$median)) * 1000
  )

  # Tight frame (regular)
  cat("  Tight frame (regular)...\n")
  cheby_cache_clear()
  b_tight <- mark(
    tight_regular = {
      bank <- tight_frame_regular(g)
      filter_bank_apply(g, signal, bank)
    },
    iterations = reps,
    check = FALSE
  )
  results$tight_regular <- data.frame(
    wavelet = "tight_regular",
    median_ms = as.numeric(median(b_tight$median)) * 1000
  )

  do.call(rbind, results)
}

bench_sgwt_vs_pygsp <- function(n, scales = 4, reps = 5) {
  g_r <- graph_ring(n)
  set.seed(42)
  signal <- rnorm(n)

  # R benchmark - SGWT
  b_r <- bench_sgwt_r(g_r, signal, scales, reps)
  W <- sgwt(g_r, signal, scales = scales)

  # R benchmark - ISGWT
  b_r_inv <- bench_isgwt_r(W, reps)

  result <- data.frame(
    n = n,
    scales = scales,
    rgsp_sgwt_ms = as.numeric(median(b_r$median)) * 1000,
    rgsp_isgwt_ms = as.numeric(median(b_r_inv$median)) * 1000
  )

  # PyGSP comparison if available
  if (pygsp_available) {
    g_py <- pygsp$graphs$Ring(as.integer(n))
    g_py$estimate_lmax()

    # PyGSP Mexican hat filter bank (similar to SGWT)
    filt_py <- pygsp$filters$MexicanHat(g_py, Nf = as.integer(scales))

    signal_py <- np$array(signal)

    # Forward transform
    b_py <- mark(
      pygsp_analysis = filt_py$filter(signal_py),
      iterations = reps,
      check = FALSE
    )

    # Get coefficients for synthesis
    coeffs_py <- filt_py$filter(signal_py)

    # Inverse transform
    b_py_inv <- mark(
      pygsp_synthesis = filt_py$synthesize(coeffs_py),
      iterations = reps,
      check = FALSE
    )

    result$pygsp_analysis_ms <- as.numeric(median(b_py$median)) * 1000
    result$pygsp_synthesis_ms <- as.numeric(median(b_py_inv$median)) * 1000
    result$speedup_analysis <- result$pygsp_analysis_ms / result$rgsp_sgwt_ms
    result$speedup_synthesis <- result$pygsp_synthesis_ms / result$rgsp_isgwt_ms
  }

  result
}

# ==============================================================================
# Main Benchmark Suite
# ==============================================================================

main <- function() {
  cat(strrep("=", 70), "\n")
  cat("Wavelet Transform Benchmark Suite\n")
  cat(strrep("=", 70), "\n\n")

  all_results <- list()

  # 1. SGWT scaling with graph size
  cat("1. SGWT scaling with graph size (scales=4)\n")
  cat(strrep("-", 50), "\n")
  sizes <- c(64, 128, 256, 512, 1024)
  res_n <- lapply(sizes, function(n) {
    cat(sprintf("  n = %d\n", n))
    bench_sgwt_vs_pygsp(n, scales = 4, reps = 5)
  })
  res_n <- do.call(rbind, res_n)
  all_results$scaling_n <- res_n
  print(res_n, row.names = FALSE)
  cat("\n")

  # 2. SGWT scaling with number of scales
  cat("2. SGWT scaling with number of scales (n=256)\n")
  cat(strrep("-", 50), "\n")
  scale_counts <- c(2, 3, 4, 5, 6)
  res_scales <- lapply(scale_counts, function(s) {
    cat(sprintf("  scales = %d\n", s))
    bench_sgwt_vs_pygsp(256, scales = s, reps = 5)
  })
  res_scales <- do.call(rbind, res_scales)
  all_results$scaling_scales <- res_scales
  print(res_scales, row.names = FALSE)
  cat("\n")

  # 3. Different wavelet banks
  cat("3. Different wavelet banks (n=256)\n")
  cat(strrep("-", 50), "\n")
  res_banks <- bench_wavelet_banks(256, reps = 5)
  all_results$wavelet_banks <- res_banks
  print(res_banks, row.names = FALSE)
  cat("\n")

  # Summary
  cat(strrep("=", 70), "\n")
  cat("Summary\n")
  cat(strrep("=", 70), "\n")

  if (pygsp_available && "speedup_analysis" %in% names(res_n)) {
    avg_speedup <- mean(c(res_n$speedup_analysis, res_scales$speedup_analysis), na.rm = TRUE)
    cat(sprintf("Average speedup vs PyGSP (analysis): %.2fx\n", avg_speedup))
  }

  # Save results
  outfile <- "inst/bench/results_wavelets.rds"
  if (file.exists(dirname(outfile))) {
    saveRDS(all_results, outfile)
    cat("Results saved to:", outfile, "\n")
  }

  invisible(all_results)
}

if (identical(environment(), globalenv())) {
  main()
}
