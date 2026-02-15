# Golden-value parity tests: classification Tikhonov against PyGSP

test_that("Golden classification_tikhonov: output logits match PyGSP", {
  ref <- skip_if_no_golden("golden_learning.json")
  entry <- ref$classification_ring20

  g <- graph_ring(20)
  labels <- entry$labels + 1L  # PyGSP 0-indexed -> R 1-indexed
  tau <- entry$tau

  scores <- classification_tikhonov(g, labels, tau = tau)

  # rgsp returns n x K matrix, PyGSP returns K columns as separate arrays
  expect_equal(ncol(scores), entry$n_classes)
  expect_equal(nrow(scores), 20L)

  for (i in seq_len(entry$n_classes)) {
    expect_equal(as.numeric(scores[, i]), entry$output_logits[i, ],
                 tolerance = tol_tikhonov,
                 label = paste("class", i, "logits"))
  }
})

test_that("Golden classification_tikhonov: predicted classes match PyGSP", {
  ref <- skip_if_no_golden("golden_learning.json")
  entry <- ref$classification_ring20

  g <- graph_ring(20)
  labels <- entry$labels + 1L
  tau <- entry$tau

  scores <- classification_tikhonov(g, labels, tau = tau)
  predicted <- apply(scores, 1, which.max)

  # PyGSP predicted_classes are 0-indexed; convert to 1-indexed
  expected <- entry$predicted_classes + 1L

  expect_equal(predicted, expected)
})

test_that("Golden classification_tikhonov: rows sum to ~1 (partition of unity)", {
  ref <- skip_if_no_golden("golden_learning.json")
  entry <- ref$classification_ring20

  g <- graph_ring(20)
  labels <- entry$labels + 1L
  tau <- entry$tau

  scores <- classification_tikhonov(g, labels, tau = tau)
  row_sums <- rowSums(scores)

  expect_equal(row_sums, rep(1, 20), tolerance = 0.01)
})
