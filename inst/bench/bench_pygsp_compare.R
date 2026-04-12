# Benchmark local rgsp working tree vs vendored PyGSP on matched operations.
#
# Usage:
#   Rscript inst/bench/bench_pygsp_compare.R
#   BENCH_OUT=inst/bench/results_pygsp_compare.csv Rscript inst/bench/bench_pygsp_compare.R

suppressPackageStartupMessages({
  library(pkgload)
  library(reticulate)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
pkg_root <- dirname(dirname(dirname(script_path)))

pkgload::load_all(pkg_root, quiet = TRUE)
source(file.path(pkg_root, "inst", "bench", "pygsp_setup.R"))

py <- init_pygsp_bench()
if (is.null(py)) {
  message("PyGSP not available; skipping.")
  quit(status = 0)
}

message("Using Python: ", py$python)
if (!is.null(py$pygsp_root)) {
  message("Using vendored PyGSP at: ", py$pygsp_root)
}

set.seed(123)
gft_sizes <- c(512L, 2048L)
filter_cases <- list(
  list(graph = "ring", n = 512L, exact_ref = TRUE),
  list(graph = "ring", n = 2048L, exact_ref = TRUE),
  list(graph = "ring", n = 8192L, exact_ref = FALSE),
  list(graph = "path", n = 512L, exact_ref = TRUE),
  list(graph = "path", n = 2048L, exact_ref = FALSE),
  list(graph = "grid2d", n1 = 32L, n2 = 32L, exact_ref = TRUE),
  list(graph = "grid2d", n1 = 64L, n2 = 64L, exact_ref = FALSE),
  list(graph = "grid2d", n1 = 96L, n2 = 96L, exact_ref = FALSE),
  list(graph = "sensor", n = 512L, k = 6L, seed = 42L, exact_ref = FALSE),
  list(graph = "sensor", n = 2048L, k = 6L, seed = 42L, exact_ref = FALSE),
  list(graph = "random_geometric", n = 512L, radius = 0.09, seed = 7L, exact_ref = TRUE),
  list(graph = "sbm", block_sizes = c(512L, 512L), P = matrix(c(0.02, 0.002, 0.002, 0.02), 2, 2), seed = 11L, exact_ref = FALSE)
)
signal_counts <- c(1L, 8L)
heat_scale <- 1.0
heat_order <- 30L

graph_label <- function(case) {
  if (identical(case$graph, "grid2d")) {
    return(sprintf("grid2d_%dx%d", case$n1, case$n2))
  }
  if (identical(case$graph, "sensor")) {
    return(sprintf("sensor_%d_k%d", case$n, case$k))
  }
  if (identical(case$graph, "random_geometric")) {
    return(sprintf("random_geometric_%d_r%.3f", case$n, case$radius))
  }
  if (identical(case$graph, "sbm")) {
    return(sprintf("sbm_%d_blocks%d", sum(case$block_sizes), length(case$block_sizes)))
  }
  sprintf("%s_%d", case$graph, case$n)
}

graph_order <- function(case) {
  if (identical(case$graph, "grid2d")) {
    return(case$n1 * case$n2)
  }
  if (identical(case$graph, "sbm")) {
    return(sum(case$block_sizes))
  }
  case$n
}

make_graph_pair <- function(case) {
  if (identical(case$graph, "ring")) {
    return(list(
      r = graph_ring(case$n),
      py = py$graphs$Ring(as.integer(case$n))
    ))
  }
  if (identical(case$graph, "path")) {
    return(list(
      r = graph_path(case$n),
      py = py$graphs$Path(as.integer(case$n))
    ))
  }
  if (identical(case$graph, "grid2d")) {
    return(list(
      r = graph_grid2d(case$n1, case$n2, periodic = FALSE),
      py = py$graphs$Grid2d(as.integer(case$n1), as.integer(case$n2))
    ))
  }
  if (identical(case$graph, "sensor")) {
    g_r <- graph_sensor(case$n, k = case$k, seed = case$seed)
    return(list(
      r = g_r,
      py = pygsp_graph_from_adjacency(py, g_r$adjacency)
    ))
  }
  if (identical(case$graph, "random_geometric")) {
    g_r <- graph_random_geometric(case$n, radius = case$radius, seed = case$seed)
    return(list(
      r = g_r,
      py = pygsp_graph_from_adjacency(py, g_r$adjacency)
    ))
  }
  if (identical(case$graph, "sbm")) {
    g_r <- graph_sbm(case$block_sizes, case$P, seed = case$seed)
    return(list(
      r = g_r,
      py = pygsp_graph_from_adjacency(py, g_r$adjacency)
    ))
  }
  stop("Unsupported graph case: ", case$graph, call. = FALSE)
}

bench_gft_case <- function(n, n_signals) {
  pair <- make_graph_pair(list(graph = "ring", n = n))
  X_r <- matrix(rnorm(n * n_signals), nrow = n)
  X_py <- py$np$array(X_r)

  graph_eigenpairs(pair$r)
  pair$py$compute_fourier_basis()

  coeffs_r <- gft(pair$r, X_r)
  coeffs_py <- pair$py$gft(X_py)
  rg_gft <- bench_op_details(function() gft(pair$r, X_r))
  rg_igft <- bench_op_details(function() igft(pair$r, coeffs_r))
  py_gft <- bench_op_details(function() pair$py$gft(X_py))
  py_igft <- bench_op_details(function() pair$py$igft(coeffs_py))

  data.frame(
    graph = "ring",
    op = c("gft", "igft"),
    n = n,
    n_signals = n_signals,
    rgsp_ms = c(rg_gft$median_ms, rg_igft$median_ms),
    rgsp_mad_ms = c(rg_gft$mad_ms, rg_igft$mad_ms),
    pygsp_ms = c(py_gft$median_ms, py_igft$median_ms),
    pygsp_mad_ms = c(py_gft$mad_ms, py_igft$mad_ms),
    rel_err = 0,
    rel_err_vs_pygsp = 0,
    lmax_py = NA_real_,
    stringsAsFactors = FALSE
  )
}

bench_heat_case <- function(case, n_signals) {
  pair <- make_graph_pair(case)
  n <- graph_order(case)
  X_r <- matrix(rnorm(n * n_signals), nrow = n)
  X_py <- py$np$array(X_r)

  pair$py$estimate_lmax()
  lmax_py <- as.numeric(pair$py$lmax)
  t_rgsp <- heat_scale / lmax_py
  heat_py <- py$filters$Heat(pair$py, scale = heat_scale)

  rg_cheby <- function() {
    filter_signal(
      pair$r, X_r, kernel_heat(t_rgsp),
      K = heat_order, lmax = lmax_py, strategy = "cheby"
    )
  }
  rg_lanczos <- function() {
    filter_signal(
      pair$r, X_r, kernel_heat(t_rgsp),
      K = heat_order, lmax = lmax_py, strategy = "lanczos"
    )
  }
  rg_exact <- function() {
    filter_signal(
      pair$r, X_r, kernel_heat(t_rgsp),
      K = heat_order, lmax = lmax_py, strategy = "exact"
    )
  }
  rg_auto <- function() {
    filter_signal(
      pair$r, X_r, kernel_heat(t_rgsp),
      K = heat_order, lmax = lmax_py, strategy = "auto"
    )
  }
  py_cheby <- function() {
    heat_py$filter(X_py, method = "chebyshev", order = as.integer(heat_order))
  }

  y_rg_cheby <- rg_cheby()
  y_rg_lanczos <- rg_lanczos()
  y_rg_exact <- rg_exact()
  y_rg_auto <- rg_auto()
  y_py_cheby <- reticulate::py_to_r(py_cheby())
  rg_cheby_bench <- bench_op_details(rg_cheby)
  rg_lanczos_bench <- bench_op_details(rg_lanczos)
  rg_exact_bench <- bench_op_details(rg_exact)
  rg_auto_bench <- bench_op_details(rg_auto)
  py_cheby_bench <- bench_op_details(py_cheby)

  if (isTRUE(case$exact_ref)) {
    pair$py$compute_fourier_basis()
    y_exact <- reticulate::py_to_r(heat_py$filter(X_py, method = "exact"))
    err_rg_cheby <- rel_l2(y_rg_cheby, y_exact)
    err_rg_lanczos <- rel_l2(y_rg_lanczos, y_exact)
    err_rg_exact <- rel_l2(y_rg_exact, y_exact)
    err_rg_auto <- rel_l2(y_rg_auto, y_exact)
    err_py_cheby <- rel_l2(y_py_cheby, y_exact)
  } else {
    err_rg_cheby <- NA_real_
    err_rg_lanczos <- NA_real_
    err_rg_exact <- NA_real_
    err_rg_auto <- NA_real_
    err_py_cheby <- NA_real_
  }

  data.frame(
    graph = graph_label(case),
    op = c("heat_cheby", "heat_lanczos", "heat_exact", "heat_auto", "pygsp_heat_cheby"),
    n = n,
    n_signals = n_signals,
    rgsp_ms = c(
      rg_cheby_bench$median_ms,
      rg_lanczos_bench$median_ms,
      rg_exact_bench$median_ms,
      rg_auto_bench$median_ms,
      NA_real_
    ),
    rgsp_mad_ms = c(
      rg_cheby_bench$mad_ms,
      rg_lanczos_bench$mad_ms,
      rg_exact_bench$mad_ms,
      rg_auto_bench$mad_ms,
      NA_real_
    ),
    pygsp_ms = c(
      py_cheby_bench$median_ms,
      py_cheby_bench$median_ms,
      py_cheby_bench$median_ms,
      py_cheby_bench$median_ms,
      py_cheby_bench$median_ms
    ),
    pygsp_mad_ms = c(
      py_cheby_bench$mad_ms,
      py_cheby_bench$mad_ms,
      py_cheby_bench$mad_ms,
      py_cheby_bench$mad_ms,
      py_cheby_bench$mad_ms
    ),
    rel_err = c(err_rg_cheby, err_rg_lanczos, err_rg_exact, err_rg_auto, err_py_cheby),
    rel_err_vs_pygsp = c(
      rel_l2(y_rg_cheby, y_py_cheby),
      rel_l2(y_rg_lanczos, y_py_cheby),
      rel_l2(y_rg_exact, y_py_cheby),
      rel_l2(y_rg_auto, y_py_cheby),
      0
    ),
    lmax_py = lmax_py,
    stringsAsFactors = FALSE
  )
}

rows <- list()

for (n in gft_sizes) {
  for (n_signals in signal_counts) {
    rows[[length(rows) + 1L]] <- bench_gft_case(n, n_signals)
  }
}

for (case in filter_cases) {
  for (n_signals in signal_counts) {
    rows[[length(rows) + 1L]] <- bench_heat_case(case, n_signals)
  }
}

df <- do.call(rbind, rows)
df$speedup <- with(df, pygsp_ms / rgsp_ms)
df$stability_ratio <- with(df, ifelse(is.na(rgsp_mad_ms) | is.na(pygsp_mad_ms), NA_real_, (rgsp_mad_ms + pygsp_mad_ms) / pmax(rgsp_ms + pygsp_ms, 1e-9)))

print(df, row.names = FALSE)

out <- Sys.getenv("BENCH_OUT", NA_character_)
if (!is.na(out) && nzchar(out)) {
  write.csv(df, out, row.names = FALSE)
  message("Wrote results to ", out)
}
