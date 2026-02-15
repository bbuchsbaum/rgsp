# Ensure package code is available for tests
# During R CMD check, the package is already installed
library(rgsp)

# If reticulate is present, eagerly initialize Python so py_available(initialize = FALSE)
# reflects the configured interpreter (helps parity tests that skip on py_available(FALSE)).
if (requireNamespace("reticulate", quietly = TRUE)) {
  try(reticulate::py_available(initialize = TRUE), silent = TRUE)
}
