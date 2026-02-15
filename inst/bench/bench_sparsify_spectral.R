# Benchmark spectral sparsification speedups
#
# Usage:
#   Rscript inst/bench/bench_sparsify_spectral.R [--quick]
#
# Measures:
# - edge count reduction from graph_sparsify_spectral()
# - downstream filter_signal() speedup (Chebyshev)

suppressPackageStartupMessages({
  library(bench)
  library(Matrix)
  library(rgsp)
})

args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args

edge_count <- function(A) length(Matrix::summary(Matrix::triu(A, 1))$x)

set.seed(1)
n <- if (quick_mode) 1000 else 2000
k <- if (quick_mode) 20 else 30
epsilon <- 0.5
K <- 50
reps <- if (quick_mode) 3 else 5

coords <- cbind(runif(n), runif(n))
g <- graph_knn(coords, k = k, weight = "binary", normalized = FALSE, seed = 1)

cat(strrep("=", 70), "\n")
cat("Spectral sparsification benchmark\n")
cat(sprintf("Mode: %s\n", if (quick_mode) "QUICK" else "FULL"))
cat(sprintf("n=%d, k=%d, epsilon=%.2f, K=%d, reps=%d\n", n, k, epsilon, K, reps))
cat(strrep("=", 70), "\n\n")

e_in <- edge_count(g$adjacency)
cat(sprintf("edges_in:  %d\n", e_in))

sp_t <- system.time({
  gs <- graph_sparsify_spectral(g, epsilon = epsilon, seed = 1, n_probes = 24)
})

e_out <- edge_count(gs$adjacency)
cat(sprintf("edges_out: %d\n", e_out))
cat(sprintf("edge reduction: %.2fx\n", e_in / max(1, e_out)))
cat(sprintf("sparsify elapsed: %.3f sec\n\n", sp_t[["elapsed"]]))

# Warm caches so we time ~O(nnz) parts.
x <- rnorm(n)
kern <- kernel_heat(1)
lmax_in <- lambda_max(g)
lmax_out <- lambda_max(gs)
invisible(filter_signal(g, x, kern, K = K, lmax = lmax_in))
invisible(filter_signal(gs, x, kern, K = K, lmax = lmax_out))

b <- mark(
  original = filter_signal(g, x, kern, K = K, lmax = lmax_in),
  sparsified = filter_signal(gs, x, kern, K = K, lmax = lmax_out),
  iterations = reps,
  check = FALSE
)

print(b)
cat("\n")
med_ms <- setNames(as.numeric(b$median) * 1000, b$expression)
cat(sprintf("median_ms original:  %.3f\n", med_ms[["original"]]))
cat(sprintf("median_ms sparsified: %.3f\n", med_ms[["sparsified"]]))
cat(sprintf("speedup: %.2fx\n", med_ms[["original"]] / med_ms[["sparsified"]]))

