find_vendored_pygsp <- function() {
  candidates <- c("pygsp", file.path("..", "pygsp"), file.path("..", "..", "pygsp"))
  for (candidate in candidates) {
    init_file <- file.path(candidate, "pygsp", "__init__.py")
    if (file.exists(init_file)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  NULL
}

prepend_python_path <- function(path) {
  current <- Sys.getenv("PYTHONPATH", unset = "")
  parts <- Filter(nzchar, strsplit(current, .Platform$path.sep, fixed = TRUE)[[1]])
  if (!path %in% parts) {
    Sys.setenv(PYTHONPATH = paste(c(path, parts), collapse = .Platform$path.sep))
  }
}

init_pygsp_bench <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(NULL)
  }

  py_bin <- Sys.getenv("RETICULATE_PYTHON", unset = Sys.which("python"))
  if (!nzchar(py_bin)) {
    return(NULL)
  }

  pygsp_root <- find_vendored_pygsp()
  if (!is.null(pygsp_root)) {
    prepend_python_path(pygsp_root)
  }

  Sys.setenv(RETICULATE_PYTHON = py_bin)
  reticulate::use_python(py_bin, required = FALSE)

  if (!reticulate::py_available(initialize = TRUE)) {
    return(NULL)
  }

  if (!is.null(pygsp_root)) {
    reticulate::py_run_string(
      sprintf(
        "import sys\np = %s\nif p not in sys.path:\n    sys.path.insert(0, p)\n",
        dQuote(pygsp_root)
      )
    )
  }

  if (!reticulate::py_module_available("pygsp")) {
    return(NULL)
  }

  list(
    python = py_bin,
    pygsp_root = pygsp_root,
    graphs = reticulate::import("pygsp.graphs", delay_load = FALSE),
    filters = reticulate::import("pygsp.filters", delay_load = FALSE),
    np = reticulate::import("numpy", delay_load = FALSE),
    scipy_sparse = reticulate::import("scipy.sparse", delay_load = FALSE)
  )
}

as_benchmark_dgC <- function(adjacency) {
  adjacency <- Matrix::drop0(adjacency)
  if (methods::is(adjacency, "dsCMatrix")) {
    adjacency <- methods::as(adjacency, "generalMatrix")
    adjacency <- methods::as(adjacency, "CsparseMatrix")
  } else if (!inherits(adjacency, "dgCMatrix")) {
    adjacency <- methods::as(adjacency, "dgCMatrix")
  }
  adjacency
}

pygsp_graph_from_adjacency <- function(py, adjacency) {
  adjacency <- as_benchmark_dgC(adjacency)
  stopifnot(inherits(adjacency, "dgCMatrix"))
  pymain <- reticulate::import_main()
  pymain$scipy_sparse <- py$scipy_sparse
  pymain$x <- as.numeric(adjacency@x)
  pymain$i <- as.integer(adjacency@i)
  pymain$p <- as.integer(adjacency@p)
  pymain$n0 <- as.integer(nrow(adjacency))
  pymain$n1 <- as.integer(ncol(adjacency))
  reticulate::py_run_string("A_py = scipy_sparse.csc_matrix((x, i, p), shape=(n0, n1))\n")
  py$graphs$Graph(pymain$A_py)
}

bench_op <- function(fun, warmup = 2L, reps = 7L) {
  for (i in seq_len(warmup)) fun()
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    gc()
    times[[i]] <- system.time(fun())[["elapsed"]]
  }
  1000 * median(times)
}

bench_op_details <- function(fun, warmup = 2L, reps = 7L) {
  for (i in seq_len(warmup)) fun()
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    gc()
    times[[i]] <- system.time(fun())[["elapsed"]]
  }
  times_ms <- 1000 * times
  list(
    median_ms = stats::median(times_ms),
    mad_ms = stats::mad(times_ms, constant = 1),
    min_ms = min(times_ms),
    max_ms = max(times_ms),
    samples_ms = times_ms
  )
}

rel_l2 <- function(x, y) {
  sqrt(sum((x - y) ^ 2)) / sqrt(sum(y ^ 2))
}
