test_that("summary.metahunt returns the right class and fields", {
  set.seed(1)
  G <- 25; m <- 40
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  fit <- metahunt(F_hat, W, K = 3,
                  dfspa_args = list(denoise = FALSE))
  s <- summary(fit)

  expect_s3_class(s, "summary.metahunt")
  expected_fields <- c("m", "G_grid", "K", "weight_method",
                       "predictor_names", "pi_summary",
                       "n_kept", "n_dropped", "denoising")
  expect_true(all(expected_fields %in% names(s)))

  expect_equal(s$m, m)
  expect_equal(s$G_grid, G)
  expect_equal(s$K, 3L)
  expect_equal(s$weight_method, fit$weight_model$method)
  expect_equal(s$predictor_names, fit$weight_model$predictor_names)
  expect_true(is.list(s$denoising))
  expect_true(all(c("N", "Delta") %in% names(s$denoising)))
})

test_that("print.summary.metahunt runs without error", {
  set.seed(2)
  G <- 25; m <- 40
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  fit <- metahunt(F_hat, W, K = 3,
                  dfspa_args = list(denoise = FALSE))
  s <- summary(fit)

  expect_output(print(s), "MetaHunt fit summary")
  expect_output(print(s), "Per-basis pi_hat summary")
})

test_that("pi_summary has shape K x 5 with correct column names", {
  set.seed(3)
  G <- 30; m <- 50
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  fit <- metahunt(F_hat, W, K = 3,
                  dfspa_args = list(denoise = FALSE))
  s <- summary(fit)

  expect_true(is.matrix(s$pi_summary))
  expect_equal(dim(s$pi_summary), c(3L, 5L))
  expect_equal(colnames(s$pi_summary),
               c("min", "mean", "median", "max", "sd"))
})

test_that("n_kept + n_dropped equals m", {
  set.seed(4)
  G <- 25; m <- 40
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  fit <- metahunt(F_hat, W, K = 3,
                  dfspa_args = list(denoise = FALSE))
  s <- summary(fit)

  expect_equal(s$n_kept + s$n_dropped, s$m)
})

test_that("coef.metahunt_weight_model returns a non-empty numeric object", {
  set.seed(5)
  G <- 25; m <- 50
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  fit <- metahunt(F_hat, W, K = 3,
                  dfspa_args = list(denoise = FALSE))
  cf <- coef(fit$weight_model)

  # DirichletReg::coef returns a list of numeric vectors (one per category
  # in the "common" parametrization), or a numeric vector / matrix in
  # other parametrizations.
  expect_true(is.numeric(cf) || is.matrix(cf) || is.list(cf))
  expect_true(length(cf) > 0L)
  expect_true(all(is.finite(unlist(cf))))
})
