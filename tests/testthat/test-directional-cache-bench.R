test_that("cheby cache hit/miss counts stable across repeated steer calls", {
  g <- graph_grid2d(5, 5, normalized = FALSE)
  g <- graph_set_coords(g, kind = g$coords)
  x <- rnorm(g$n)
  scales <- c(1, 2)

  cheby_cache_clear()
  cache_env <- get(".cheby_cache", envir = asNamespace("rgsp"))
  count_keys <- function() length(ls(envir = cache_env))

  dsgwt_steer(g, x, n_directions = 2, scales = scales,
              wavelet = "mexican_hat", include_lowpass = TRUE, K = 6,
              alpha_isotropic = 0.2, p = 2)
  first <- count_keys()

  dsgwt_steer(g, x, n_directions = 4, scales = scales,
              wavelet = "mexican_hat", include_lowpass = TRUE, K = 6,
              alpha_isotropic = 0.2, p = 2)
  second <- count_keys()
  expect_equal(second, first)
})
