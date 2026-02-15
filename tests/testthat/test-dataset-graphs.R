test_that("graph_logo loads the PyGSP Logo dataset graph", {
  g <- graph_logo()
  expect_s3_class(g, "gsp_graph")
  expect_true(is.matrix(g$coords))
  expect_equal(nrow(g$coords), g$n)
  expect_true(is.list(g$info))
  expect_true(is.list(g$plotting))
})

test_that("graph_minnesota loads the PyGSP Minnesota dataset graph", {
  g <- graph_minnesota(connected = TRUE)
  expect_s3_class(g, "gsp_graph")
  expect_true(is.matrix(g$coords))
  expect_equal(nrow(g$coords), g$n)
  expect_equal(length(g$labels), g$n)

  expect_true(rgsp:::.graph_is_connected(g$adjacency))

  # Connected variant is binarized.
  x <- g$adjacency@x
  expect_true(all(abs(x - 1) < 1e-12))
})

test_that("graph_airfoil loads the PyGSP Airfoil dataset graph", {
  g <- graph_airfoil()
  expect_s3_class(g, "gsp_graph")
  expect_equal(g$n, 4253L)
  expect_true(is.matrix(g$coords))
  expect_equal(nrow(g$coords), g$n)
  expect_equal(ncol(g$coords), 2L)
})

test_that("graph_bunny loads the PyGSP Bunny dataset graph", {
  g <- graph_bunny()
  expect_s3_class(g, "gsp_graph")
  expect_true(is.matrix(g$coords))
  expect_equal(nrow(g$coords), g$n)
  expect_equal(ncol(g$coords), 3L)
  expect_true(is.list(g$params))
  expect_equal(g$params$epsilon, 0.02, tolerance = 1e-12)
})

test_that("graph_two_moons standard produces labels and coords", {
  g <- graph_two_moons("standard", dim = 2, sigmag = 0.05)
  expect_s3_class(g, "gsp_graph")
  expect_true(is.matrix(g$coords))
  expect_equal(nrow(g$coords), g$n)
  expect_equal(ncol(g$coords), 2L)
  expect_equal(length(g$labels), g$n)
  expect_equal(sort(unique(g$labels)), c(0L, 1L))
})

test_that("graph_two_moons synthesized works for smaller N", {
  g <- graph_two_moons("synthesized", N = 200, seed = 1)
  expect_s3_class(g, "gsp_graph")
  expect_equal(g$n, 200L)
  expect_equal(length(g$labels), 200L)
})

test_that("graph_david_sensor_net loads and generates graphs", {
  g64 <- graph_david_sensor_net(64)
  expect_equal(g64$n, 64L)
  expect_true(is.matrix(g64$coords))

  g80 <- graph_david_sensor_net(80, seed = 1)
  expect_equal(g80$n, 80L)
  expect_true(is.matrix(g80$coords))
})
