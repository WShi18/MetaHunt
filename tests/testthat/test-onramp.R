test_that("f_hat_from_models stacks lm predictions correctly", {
  set.seed(1)
  make <- function(slope) {
    x <- runif(60)
    data.frame(x = x, y = slope * x + rnorm(60, sd = 0.05))
  }
  models <- lapply(c(-1, 0, 1, 0.5), function(s) stats::lm(y ~ poly(x, 2),
                                                           data = make(s)))
  grid <- data.frame(x = seq(0, 1, length.out = 25))
  F_hat <- f_hat_from_models(models, grid)
  expect_equal(dim(F_hat), c(4L, 25L))
  expect_equal(F_hat[1, ],
               as.numeric(stats::predict(models[[1]], newdata = grid)))
})

test_that("f_hat_from_models honours custom predict_fn", {
  models <- list(
    function(x) sin(x),
    function(x) cos(x),
    function(x) x^2
  )
  grid <- data.frame(x = seq(0, 1, length.out = 10))
  F_hat <- f_hat_from_models(
    models, grid,
    predict_fn = function(model, grid) model(grid$x)
  )
  expect_equal(F_hat[1, ], sin(grid$x))
  expect_equal(F_hat[2, ], cos(grid$x))
  expect_equal(F_hat[3, ], grid$x^2)
})

test_that("f_hat_from_models errors when predict lengths disagree", {
  models <- list(
    function(x) x[1:5],
    function(x) x
  )
  grid <- data.frame(x = 1:10)
  expect_error(
    f_hat_from_models(models, grid,
                      predict_fn = function(m, g) m(g$x)),
    "different-length|returned 5 predictions"
  )
})

test_that("f_hat_from_models errors on NA predictions", {
  models <- list(function(x) c(NA_real_, x[-1]))
  grid <- data.frame(x = 1:5)
  expect_error(
    f_hat_from_models(models, grid,
                      predict_fn = function(m, g) m(g$x)),
    "NA"
  )
})

test_that("build_grid returns reference unchanged when n_grid is NULL or large", {
  ref <- data.frame(a = rnorm(40), b = rnorm(40))
  expect_identical(build_grid(ref), ref)
  expect_identical(build_grid(ref, n_grid = 100), ref)
})

test_that("build_grid samples without replacement and is reproducible", {
  ref <- data.frame(a = rnorm(80), b = rnorm(80))
  g1 <- build_grid(ref, n_grid = 25, seed = 7)
  g2 <- build_grid(ref, n_grid = 25, seed = 7)
  expect_equal(g1, g2)
  expect_equal(nrow(g1), 25L)
  expect_equal(anyDuplicated(rownames(g1)), 0L)
})

test_that("build_grid rejects malformed inputs", {
  expect_error(build_grid("nope"), "matrix or data frame")
  expect_error(build_grid(data.frame(a = double(0))), "empty")
  expect_error(build_grid(data.frame(a = 1:5), n_grid = -1),
               "positive integer")
})
