# Ensure package code is available for tests
# During R CMD check, the package is already installed
library(rgsp)

# If reticulate is present, prefer the vendored ./pygsp tree with the system
# Python before eager initialization. This keeps live parity tests runnable in
# clean environments where PyGSP is not installed site-wide.
if (requireNamespace("reticulate", quietly = TRUE)) {
  pygsp_candidates <- c("pygsp", file.path("..", "pygsp"), file.path("..", "..", "pygsp"))
  vendored_pygsp <- ""
  for (candidate in pygsp_candidates) {
    candidate_path <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (dir.exists(candidate_path)) {
      vendored_pygsp <- candidate_path
      break
    }
  }
  if (dir.exists(vendored_pygsp)) {
    if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
      py_bin <- Sys.which("python")
      if (nzchar(py_bin)) {
        Sys.setenv(RETICULATE_PYTHON = py_bin)
        reticulate::use_python(py_bin, required = FALSE)
      }
    }

    py_path <- Sys.getenv("PYTHONPATH")
    py_parts <- Filter(nzchar, strsplit(py_path, .Platform$path.sep, fixed = TRUE)[[1]])
    if (!(vendored_pygsp %in% py_parts)) {
      Sys.setenv(PYTHONPATH = paste(c(vendored_pygsp, py_parts), collapse = .Platform$path.sep))
    }
  }

  try(reticulate::py_available(initialize = TRUE), silent = TRUE)
}
