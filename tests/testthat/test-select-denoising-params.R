test_that("select_denoising_params returns the expected structure", {
  set.seed(1)
  G <- 25; m <- 60; K_true <- 3
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.5), c(-0.4, 1), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  res <- select_denoising_params(F_hat, W, K = 3,
                                 N_grid = c(0.2, 0.5) * log(m),
                                 Delta_grid = c(0.1, 0.2),
                                 n_folds = 3L, seed = 1)

  expect_type(res, "list")
  expect_s3_class(res, "metahunt_denoising_search")
  expect_named(res$grid,
               c("N", "Delta", "cv_error", "cv_se", "n_folds_ok"),
               ignore.order = TRUE)
  expect_equal(nrow(res$grid), 4L)
  expect_named(res$best, c("N", "Delta", "cv_error"))
  expect_true(res$best$cv_error == min(res$grid$cv_error, na.rm = TRUE))
})

test_that("print.metahunt_denoising_search runs and prints a summary", {
  set.seed(11)
  G <- 25; m <- 60
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.5), c(-0.4, 1), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  res <- select_denoising_params(F_hat, W, K = 3,
                                 N_grid = c(0.2, 0.5) * log(m),
                                 Delta_grid = c(0.1, 0.2),
                                 n_folds = 3L, seed = 1)

  expect_output(print(res), "Denoising parameter search")
  expect_invisible(print(res))
})

test_that("select_denoising_params honours custom grids", {
  set.seed(2)
  G <- 20; m <- 50
  F_hat <- matrix(rnorm(m * G), m, G)
  W <- data.frame(w1 = rnorm(m))

  res <- suppressWarnings(
    select_denoising_params(F_hat, W, K = 2,
                            N_grid = c(1, 2),
                            Delta_grid = c(0.5, 1, 2),
                            n_folds = 3L, seed = 1)
  )
  expect_setequal(unique(res$grid$N), c(1, 2))
  expect_setequal(unique(res$grid$Delta), c(0.5, 1, 2))
  expect_equal(nrow(res$grid), 6L)
})

test_that("select_denoising_params is reproducible with a seed", {
  set.seed(3)
  G <- 20; m <- 50
  F_hat <- matrix(rnorm(m * G), m, G)
  W <- data.frame(w1 = rnorm(m))

  a <- select_denoising_params(F_hat, W, K = 2,
                               N_grid = c(0.5, 1),
                               Delta_grid = c(1, 2),
                               n_folds = 3L, seed = 99)
  b <- select_denoising_params(F_hat, W, K = 2,
                               N_grid = c(0.5, 1),
                               Delta_grid = c(1, 2),
                               n_folds = 3L, seed = 99)
  expect_equal(a$grid, b$grid)
  expect_equal(a$best, b$best)
})

test_that("select_denoising_params rejects N/Delta inside dfspa_args", {
  G <- 15; m <- 30
  F_hat <- matrix(rnorm(m * G), m, G)
  W <- data.frame(w1 = rnorm(m))
  expect_error(
    select_denoising_params(F_hat, W, K = 2, n_folds = 2L,
                            dfspa_args = list(N = 1)),
    "must not contain N"
  )
})
