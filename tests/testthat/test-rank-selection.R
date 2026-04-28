test_that("reconstruction_error_curve is non-increasing for noiseless rank-K data", {
  set.seed(1)
  G <- 40; K_true <- 3; m <- 30
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  pi_mat <- matrix(stats::rgamma(m * K_true, 0.5), m, K_true)
  pi_mat <- pi_mat / rowSums(pi_mat)
  # include near-pure studies to guarantee basis recovery
  pi_mat[1:3, ] <- diag(K_true)
  F_hat <- pi_mat %*% basis

  res <- reconstruction_error_curve(F_hat, K_range = 2:K_true,
                                    dfspa_args = list(denoise = FALSE))

  expect_s3_class(res, "data.frame")
  expect_named(res, c("K", "error"))
  expect_equal(res$K, 2:K_true)
  expect_true(all(diff(res$error) <= 1e-8))
  # at K = K_true, error should be ~0 in the noiseless case
  expect_lt(res$error[res$K == K_true], 1e-6)
})

test_that("reconstruction_error_curve respects dfspa_args and grid_weights", {
  set.seed(2)
  G <- 20; m <- 15
  F_hat <- matrix(stats::rnorm(m * G), m, G)
  w <- runif(G, min = 0.1)

  res <- reconstruction_error_curve(
    F_hat, K_range = 2:4, grid_weights = w,
    dfspa_args = list(denoise = FALSE)
  )
  expect_true(all(is.finite(res$error)))
  expect_equal(nrow(res), 3L)
})

test_that("cv_error_curve returns the expected structure", {
  set.seed(3)
  G <- 30; m <- 40; K_true <- 3
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(0.8, -0.5), c(-0.3, 0.9), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.02), m, G)

  res <- cv_error_curve(F_hat, W, K_range = 2:4, n_folds = 4L,
                       dfspa_args = list(denoise = FALSE), seed = 42)

  expect_s3_class(res, "data.frame")
  expect_named(res, c("K", "cv_error", "cv_se"))
  expect_equal(res$K, 2:4)
  fe <- attr(res, "fold_errors")
  expect_equal(dim(fe), c(3L, 4L))
  expect_true(all(is.finite(res$cv_error)))
})

test_that("cv_error_curve is reproducible with a seed", {
  set.seed(4)
  G <- 20; m <- 30
  F_hat <- matrix(rnorm(m * G), m, G)
  W <- data.frame(w1 = rnorm(m))

  a <- cv_error_curve(F_hat, W, K_range = 2:3, n_folds = 3L,
                      dfspa_args = list(denoise = FALSE), seed = 99)
  b <- cv_error_curve(F_hat, W, K_range = 2:3, n_folds = 3L,
                      dfspa_args = list(denoise = FALSE), seed = 99)
  expect_equal(a, b)
})

test_that("cv_error_curve surfaces failures as NA with a warning", {
  set.seed(5)
  G <- 20; m <- 10
  F_hat <- matrix(rnorm(m * G), m, G)
  W <- data.frame(w1 = rnorm(m))

  # Force failures by asking for K too large relative to fold training size
  expect_warning(
    res <- cv_error_curve(F_hat, W, K_range = 9:10, n_folds = 2L,
                          dfspa_args = list(denoise = FALSE), seed = 1),
    "failed"
  )
  expect_true(anyNA(res$cv_error) || anyNA(attr(res, "fold_errors")))
})
