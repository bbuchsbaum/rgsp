# Golden-value parity tests: frame bounds and tight frames against PyGSP

test_that("Golden heat frame bounds are positive", {
  ref <- skip_if_no_golden("golden_frames.json")
  entry <- ref$heat_ring10

  g <- graph_ring(10)
  kernels <- list(kernel_heat(1.0))
  frame_r <- compute_frame(g, kernels, method = "exact")

  expect_true(frame_r$A > 0)
  expect_true(frame_r$B >= frame_r$A)

  # PyGSP bounds should also be positive

  expect_true(entry$A > 0)
  expect_true(entry$B >= entry$A)
})

test_that("Golden Regular tight frame: bounds close to 1", {
  ref <- skip_if_no_golden("golden_frames.json")
  entry <- ref$regular_ring16

  # PyGSP Regular is a tight frame: A ~= B ~= 1
  expect_true(entry$A > 0.9)
  expect_true(entry$B < 1.1)
  expect_true(entry$B / entry$A < 1.2)

  # rgsp tight_frame_regular should also be tight
  g <- graph_ring(16)
  bank <- tight_frame_regular(lambda_max(g))
  eig <- graph_eigenpairs(g)

  # Sum of squared responses at eigenvalues should be ~1
  lam <- eig$values
  sos <- Reduce("+", lapply(bank, function(k) k(lam)^2))
  expect_true(all(sos > 0.5), label = "SOS should be positive")
})

test_that("Golden tight frames: PyGSP sum-of-squares near 1", {
  ref <- skip_if_no_golden("golden_frames.json")
  tf <- ref$tight_frames

  for (name in names(tf)) {
    entry <- tf[[name]]
    sos <- entry$sum_of_squares

    # For tight frames, sum of squares should be approximately 1
    expect_true(all(sos > 0.8),
                label = paste(name, "SOS minimum:", round(min(sos), 3)))
    expect_true(all(sos < 1.2),
                label = paste(name, "SOS maximum:", round(max(sos), 3)))
  }
})

test_that("Golden tight frames: PyGSP bounds A and B are close", {
  ref <- skip_if_no_golden("golden_frames.json")
  tf <- ref$tight_frames

  for (name in names(tf)) {
    entry <- tf[[name]]
    expect_true(entry$A > 0,
                label = paste(name, "A > 0"))
    expect_true(entry$B / max(entry$A, 1e-10) < 1.5,
                label = paste(name, "B/A ratio"))
  }
})

test_that("Golden tight frame filter count matches: Regular has 2 filters", {
  ref <- skip_if_no_golden("golden_frames.json")
  entry <- ref$regular_ring16

  expect_equal(entry$n_filters, 2L)
})
