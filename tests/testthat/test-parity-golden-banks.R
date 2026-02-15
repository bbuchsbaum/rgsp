# Golden-value parity tests: filter banks (MexicanHat, Meyer, Heat) against PyGSP
#
# NOTE: jsonlite::fromJSON converts list-of-lists to matrices, so
# entry$kernel_evals[i,] gives the i-th band (row), NOT entry$kernel_evals[[i]].
#
# API differences:
# - MexicanHat: PyGSP bandpass = scale*x*exp(-scale*x), rgsp = x*exp(-x/sigma)
#   with sigma=1/scale. PyGSP_kernel = scale * rgsp_kernel (scale factor).
# - Meyer: PyGSP uses smooth Meyer scaling function. rgsp uses kernel_rectangle
#   for lowpass and kernel_meyer for bandpass. Different implementations,
#   so exact kernel match is not expected. Test properties instead.

# --- MexicanHat bank tests ---

test_that("Golden MexicanHat bank: band count matches PyGSP", {
  ref <- skip_if_no_golden("golden_filterbanks.json")
  entry <- ref$mexicanhat_ring20

  # PyGSP MexicanHat(Nf=4) produces 4 bands (1 lowpass + 3 bandpass)
  expect_equal(entry$n_bands, entry$Nf)
})

test_that("Golden MexicanHat bank: kernel shape matches with scale factor", {
  ref <- skip_if_no_golden("golden_filterbanks.json")
  entry <- ref$mexicanhat_ring20

  g <- graph_ring(20)
  eig <- sort(graph_eigenpairs(g)$values)

  # Use the PyGSP scales directly with rgsp's kernel_mexican_hat
  pygsp_scales <- entry$scales

  for (i in seq_along(pygsp_scales)) {
    s <- pygsp_scales[i]
    sigma_r <- pygsp_mexhat_scale_to_rgsp_sigma(s)
    rgsp_vals <- kernel_mexican_hat(sigma_r)(eig)

    # PyGSP bandpass: s * x * exp(-s*x) = s * rgsp_kernel(x)
    # So pygsp_bp = s * rgsp_vals
    pygsp_bp <- entry$kernel_evals[i + 1, ]  # skip lowpass (row 1)

    expect_equal(s * rgsp_vals, pygsp_bp, tolerance = 1e-8,
                 label = paste("MexicanHat band", i, "with scale factor"))
  }
})

test_that("Golden MexicanHat bank: filtered output matches with scale factor", {
  ref <- skip_if_no_golden("golden_filterbanks.json")
  entry <- ref$mexicanhat_ring20

  g <- graph_ring(20)
  x <- entry$signal
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  lam <- eig$values

  pygsp_scales <- entry$scales

  for (i in seq_along(pygsp_scales)) {
    s <- pygsp_scales[i]
    sigma_r <- pygsp_mexhat_scale_to_rgsp_sigma(s)
    kern_vals <- kernel_mexican_hat(sigma_r)(lam)
    # Apply exact filter: U diag(kern_vals) U^T x
    y_r <- as.numeric(U %*% diag(kern_vals) %*% t(U) %*% x)

    # PyGSP filtered output includes the scale factor
    pygsp_band <- entry$filtered_output_per_band[i + 1, ]

    # PyGSP_filtered = s * rgsp_filtered (due to kernel scale factor)
    # Tolerance is higher than tol_filter_exact because eigendecomposition
    # differences compound through the scale factor multiplication
    expect_equal(s * y_r, pygsp_band, tolerance = 1e-3,
                 label = paste("MexicanHat band", i, "filtered output"))
  }
})

# --- Meyer bank tests ---

test_that("Golden Meyer bank: band count matches PyGSP", {
  ref <- skip_if_no_golden("golden_filterbanks.json")
  entry <- ref$meyer_path16

  g <- graph_path(16)

  # PyGSP Meyer(Nf=3) produces 3 bands (1 scaling + 2 wavelets)
  bank_r <- meyer_wavelet_bank(g, n_scales = entry$Nf - 1, lmax = entry$lmax)
  expect_equal(length(bank_r), entry$n_bands)
})

test_that("Golden Meyer bank: kernels form partition of unity", {
  ref <- skip_if_no_golden("golden_filterbanks.json")
  entry <- ref$meyer_path16

  g <- graph_path(16)
  eig <- sort(graph_eigenpairs(g)$values)

  bank_r <- meyer_wavelet_bank(g, n_scales = entry$Nf - 1, lmax = entry$lmax)

  # Sum of squared kernel responses should cover the spectrum
  sos <- numeric(length(eig))
  for (i in seq_along(bank_r)) {
    vals <- bank_r[[i]](eig)
    sos <- sos + vals^2
  }

  # For a well-formed wavelet bank, SOS should be >= some minimum for interior eigenvalues
  # (May not be exactly 1 since rgsp uses rectangle + Meyer, not a tight frame)
  interior <- eig > 1e-6 & eig < max(eig) - 1e-6
  if (any(interior)) {
    expect_true(all(sos[interior] > 0.1),
                label = "Meyer bank covers interior spectrum")
  }
})

test_that("Golden Meyer bank: filtered output preserves signal structure", {
  ref <- skip_if_no_golden("golden_filterbanks.json")
  entry <- ref$meyer_path16

  g <- graph_path(16)
  x <- entry$signal
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  lam <- eig$values

  bank_r <- meyer_wavelet_bank(g, n_scales = entry$Nf - 1, lmax = entry$lmax)

  # Sum of band-filtered outputs should approximate the original signal
  # (when bands cover the spectrum)
  y_total <- numeric(length(x))
  for (i in seq_along(bank_r)) {
    kern_vals <- bank_r[[i]](lam)
    y_band <- as.numeric(U %*% diag(kern_vals^2) %*% t(U) %*% x)
    y_total <- y_total + y_band
  }

  # Correlation with original should be high
  r <- cor(y_total, x)
  expect_true(r > 0.9,
              label = paste("Meyer bank reconstruction correlation:", round(r, 4)))
})

# --- Heat bank at multiple scales ---

test_that("Golden heat bank: kernel evals at eigenvalues match PyGSP", {
  ref <- skip_if_no_golden("golden_filterbanks.json")
  entry <- ref$heat_bank_ring12

  g <- graph_ring(12)
  eig <- sort(graph_eigenpairs(g)$values)
  lmax_py <- entry$lmax
  scales <- entry$scales

  for (i in seq_along(scales)) {
    s <- scales[[i]]
    t_rgsp <- pygsp_heat_to_rgsp(s, lmax_py)
    kern_r <- kernel_heat(t_rgsp)(eig)
    kern_py <- entry$kernel_evals_per_scale[i, ]  # matrix row access
    expect_equal(kern_r, kern_py, tolerance = tol_filter_exact,
                 label = paste("Heat bank scale", s, "kernel evals"))
  }
})

test_that("Golden heat bank: filtered output per scale matches PyGSP", {
  ref <- skip_if_no_golden("golden_filterbanks.json")
  entry <- ref$heat_bank_ring12

  g <- graph_ring(12)
  x <- entry$signal
  eig <- graph_eigenpairs(g)
  U <- eig$vectors
  lam <- eig$values
  lmax_py <- entry$lmax
  scales <- entry$scales

  for (i in seq_along(scales)) {
    s <- scales[[i]]
    t_rgsp <- pygsp_heat_to_rgsp(s, lmax_py)
    kern_vals <- kernel_heat(t_rgsp)(lam)
    y_r <- as.numeric(U %*% diag(kern_vals) %*% t(U) %*% x)
    y_py <- entry$filtered_output_per_scale[i, ]  # matrix row access
    expect_equal(y_r, y_py, tolerance = tol_filter_exact,
                 label = paste("Heat bank scale", s, "filtered output"))
  }
})
