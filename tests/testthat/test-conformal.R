make_toy <- function(seed = 1, m = 100, G = 30, K_true = 3, sd_noise = 0.05) {
  set.seed(seed)
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = sd_noise), m, G)
  list(F_hat = F_hat, W = W, basis = basis, pi_true = pi_true, x = x)
}

test_that("split_conformal pointwise: shapes and containment", {
  d <- make_toy()
  W_new <- d$W[1:5, , drop = FALSE]
  res <- split_conformal(d$F_hat, d$W, W_new, K = 3, seed = 1)

  expect_s3_class(res, "metahunt_conformal")
  expect_equal(res$method, "split")
  expect_equal(dim(res$prediction), c(5L, ncol(d$F_hat)))
  expect_equal(dim(res$lower), dim(res$prediction))
  expect_equal(dim(res$upper), dim(res$prediction))
  expect_length(res$quantile, ncol(d$F_hat))
  expect_true(all(res$lower <= res$prediction))
  expect_true(all(res$upper >= res$prediction))
  expect_equal(res$alpha, 0.05)
})

test_that("split_conformal scalar: shapes and containment", {
  d <- make_toy()
  W_new <- d$W[1:7, , drop = FALSE]
  res <- split_conformal(d$F_hat, d$W, W_new, K = 3,
                         wrapper = mean, seed = 1)

  expect_length(res$prediction, 7L)
  expect_length(res$lower, 7L)
  expect_length(res$upper, 7L)
  expect_length(res$quantile, 1L)
  expect_true(all(res$lower <= res$prediction))
  expect_true(all(res$upper >= res$prediction))
})

test_that("split_conformal achieves roughly nominal coverage on held-out data", {
  d <- make_toy(m = 200)
  # leave-out 50 studies as test, use the rest to fit + calibrate
  set.seed(2)
  test_idx <- sample(nrow(d$F_hat), 50)
  train_idx <- setdiff(seq_len(nrow(d$F_hat)), test_idx)

  res <- split_conformal(
    d$F_hat[train_idx, , drop = FALSE], d$W[train_idx, , drop = FALSE],
    d$W[test_idx, , drop = FALSE], K = 3,
    wrapper = mean, alpha = 0.1, seed = 3
  )
  obs_scalar <- apply_wrapper(d$F_hat[test_idx, , drop = FALSE],
                              wrapper = mean)
  covered <- obs_scalar >= res$lower & obs_scalar <= res$upper
  emp_cov <- mean(covered)
  # with n_cal ~ 45 and alpha = 0.1, expect ~90% empirical coverage;
  # be generous (0.75-1.0) since we only have 50 test points.
  expect_gt(emp_cov, 0.75)
})

test_that("split_conformal is reproducible with seed", {
  d <- make_toy(m = 80)
  W_new <- d$W[1:3, , drop = FALSE]
  a <- split_conformal(d$F_hat, d$W, W_new, K = 3, wrapper = mean, seed = 42)
  b <- split_conformal(d$F_hat, d$W, W_new, K = 3, wrapper = mean, seed = 42)
  expect_equal(a$prediction, b$prediction)
  expect_equal(a$lower, b$lower)
  expect_equal(a$upper, b$upper)
})

test_that("split_conformal honours calibration_idx", {
  d <- make_toy(m = 80)
  W_new <- d$W[1:3, , drop = FALSE]
  idx <- 1:20
  res <- split_conformal(d$F_hat, d$W, W_new, K = 3, wrapper = mean,
                         calibration_idx = idx)
  expect_equal(res$n_cal, length(idx))
})

test_that("split_conformal: infinite quantile when n_cal too small for alpha", {
  d <- make_toy(m = 30)
  W_new <- d$W[1:2, , drop = FALSE]
  expect_warning(
    res <- split_conformal(d$F_hat, d$W, W_new, K = 3,
                           alpha = 0.01, cal_frac = 0.1, wrapper = mean,
                           seed = 1),
    "infinite"
  )
  # with n_cal = 3 and alpha = 0.01, k = ceiling(0.99 * 4) = 4 > 3 → Inf
  expect_true(is.infinite(res$quantile))
  expect_true(all(is.infinite(res$lower) | is.infinite(res$upper)))
})

test_that("cross_conformal: shapes and basic correctness", {
  d <- make_toy(m = 60)
  W_new <- d$W[1:4, , drop = FALSE]
  res <- cross_conformal(d$F_hat, d$W, W_new, K = 3,
                         n_folds = 4L, wrapper = mean, seed = 7)

  expect_s3_class(res, "metahunt_conformal")
  expect_equal(res$method, "cross")
  expect_length(res$prediction, 4L)
  expect_true(res$n_cal <= nrow(d$F_hat))
  expect_true(all(res$lower <= res$prediction))
  expect_true(all(res$upper >= res$prediction))
})

test_that("cross_conformal pointwise: shapes match grid", {
  d <- make_toy(m = 60)
  W_new <- d$W[1:3, , drop = FALSE]
  res <- cross_conformal(d$F_hat, d$W, W_new, K = 3,
                         n_folds = 4L, seed = 9)
  expect_equal(dim(res$prediction), c(3L, ncol(d$F_hat)))
  expect_length(res$quantile, ncol(d$F_hat))
})

test_that("conformal input validation catches malformed inputs", {
  d <- make_toy(m = 40)
  expect_error(split_conformal(d$F_hat, d$W, d$W[1:2, ], K = 3, alpha = 1.5),
               "alpha.*\\(0, 1\\)")
  expect_error(split_conformal(d$F_hat, d$W, d$W[1:2, ], K = 3,
                               cal_frac = 1.2),
               "cal_frac.*\\(0, 1\\)")
  expect_error(cross_conformal(d$F_hat, d$W, d$W[1:2, ], K = 3,
                               n_folds = 1),
               "n_folds")
  expect_error(split_conformal(d$F_hat, d$W[-1, ], d$W[1:2, ], K = 3),
               "same number of rows")
})

test_that("print.metahunt_conformal runs without error", {
  d <- make_toy(m = 40)
  suppressWarnings({
    res1 <- split_conformal(d$F_hat, d$W, d$W[1:2, ], K = 3, seed = 1)
    res2 <- split_conformal(d$F_hat, d$W, d$W[1:2, ], K = 3,
                            wrapper = mean, seed = 1)
  })
  expect_output(print(res1), "MetaHunt conformal prediction")
  expect_output(print(res2), "MetaHunt conformal prediction")
})
