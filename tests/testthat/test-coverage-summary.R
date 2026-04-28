make_toy_cs <- function(seed = 1, m = 80, G = 25, sd_noise = 0.05) {
  set.seed(seed)
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = sd_noise), m, G)
  list(F_hat = F_hat, W = W, basis = basis, pi_true = pi_true, x = x)
}

test_that("summary.metahunt_conformal: pointwise object has expected fields", {
  d <- make_toy_cs(m = 80, G = 25)
  W_new <- d$W[1:5, , drop = FALSE]
  res <- split_conformal(d$F_hat, d$W, W_new, K = 3,
                         dfspa_args = list(denoise = FALSE), seed = 1)

  s <- summary(res)
  expect_s3_class(s, "summary.metahunt_conformal")
  expect_true(all(c("n_targets", "G_grid", "n_cal", "alpha", "method",
                    "mean_interval_width", "frac_finite_quantile",
                    "quantile_summary") %in% names(s)))
  expect_equal(s$n_targets, 5L)
  expect_equal(s$G_grid, ncol(d$F_hat))
  expect_equal(s$method, "split")
  expect_true(s$frac_finite_quantile >= 0 && s$frac_finite_quantile <= 1)
  expect_output(print(s), "Summary of MetaHunt conformal prediction")
})

test_that("summary.metahunt_conformal: scalar object has expected fields", {
  d <- make_toy_cs(m = 80, G = 25)
  W_new <- d$W[1:5, , drop = FALSE]
  res <- split_conformal(d$F_hat, d$W, W_new, K = 3, wrapper = mean,
                         dfspa_args = list(denoise = FALSE), seed = 1)

  s <- summary(res)
  expect_s3_class(s, "summary.metahunt_conformal")
  expect_true(all(c("n_targets", "n_cal", "alpha", "method",
                    "mean_interval_width", "quantile",
                    "quantile_finite") %in% names(s)))
  expect_equal(s$n_targets, 5L)
  expect_equal(s$method, "split")
  expect_output(print(s), "Summary of MetaHunt conformal prediction")
})

test_that("coverage(): pointwise object returns correct shapes and ranges", {
  d <- make_toy_cs(m = 80, G = 25)
  test_idx <- 1:8
  train_idx <- setdiff(seq_len(nrow(d$F_hat)), test_idx)
  res <- split_conformal(d$F_hat[train_idx, , drop = FALSE],
                         d$W[train_idx, , drop = FALSE],
                         d$W[test_idx, , drop = FALSE], K = 3,
                         dfspa_args = list(denoise = FALSE), seed = 1)

  cov <- coverage(res, F_obs = d$F_hat[test_idx, , drop = FALSE])
  expect_true(all(c("pointwise", "per_target", "per_grid_point",
                    "overall", "nominal") %in% names(cov)))
  expect_equal(dim(cov$pointwise),
               c(length(test_idx), ncol(d$F_hat)))
  expect_length(cov$per_target, length(test_idx))
  expect_length(cov$per_grid_point, ncol(d$F_hat))
  expect_true(cov$overall >= 0 && cov$overall <= 1)
  expect_true(all(cov$per_target >= 0 & cov$per_target <= 1))
  expect_true(all(cov$per_grid_point >= 0 & cov$per_grid_point <= 1))
  expect_equal(cov$nominal, 1 - res$alpha)
})

test_that("coverage(): scalar object returns correct shapes and ranges", {
  d <- make_toy_cs(m = 80, G = 25)
  test_idx <- 1:8
  train_idx <- setdiff(seq_len(nrow(d$F_hat)), test_idx)
  res <- split_conformal(d$F_hat[train_idx, , drop = FALSE],
                         d$W[train_idx, , drop = FALSE],
                         d$W[test_idx, , drop = FALSE], K = 3,
                         wrapper = mean,
                         dfspa_args = list(denoise = FALSE), seed = 1)

  cov <- coverage(res, F_obs = d$F_hat[test_idx, , drop = FALSE])
  expect_true(all(c("pointwise", "overall", "nominal") %in% names(cov)))
  expect_length(cov$pointwise, length(test_idx))
  expect_true(is.logical(cov$pointwise))
  expect_true(cov$overall >= 0 && cov$overall <= 1)
  expect_equal(cov$nominal, 1 - res$alpha)
})

test_that("coverage(): rejects mismatched F_obs shape", {
  d <- make_toy_cs(m = 80, G = 25)
  test_idx <- 1:5
  train_idx <- setdiff(seq_len(nrow(d$F_hat)), test_idx)
  res_pt <- split_conformal(d$F_hat[train_idx, , drop = FALSE],
                            d$W[train_idx, , drop = FALSE],
                            d$W[test_idx, , drop = FALSE], K = 3,
                            dfspa_args = list(denoise = FALSE), seed = 1)
  # wrong number of rows
  expect_error(
    coverage(res_pt, F_obs = d$F_hat[1:3, , drop = FALSE]),
    "rows"
  )
  # wrong number of columns
  expect_error(
    coverage(res_pt, F_obs = d$F_hat[test_idx, 1:10, drop = FALSE]),
    "columns"
  )

  res_sc <- split_conformal(d$F_hat[train_idx, , drop = FALSE],
                            d$W[train_idx, , drop = FALSE],
                            d$W[test_idx, , drop = FALSE], K = 3,
                            wrapper = mean,
                            dfspa_args = list(denoise = FALSE), seed = 1)
  expect_error(
    coverage(res_sc, F_obs = d$F_hat[1:2, , drop = FALSE]),
    "rows"
  )
})

test_that("coverage(): rejects non-metahunt_conformal input", {
  expect_error(coverage(list(prediction = 1), F_obs = matrix(0, 1, 1)),
               "metahunt_conformal")
  expect_error(coverage(42, F_obs = matrix(0, 1, 1)),
               "metahunt_conformal")
})

test_that("coverage(): empirical coverage is roughly close to nominal", {
  d <- make_toy_cs(m = 200, G = 25)
  set.seed(11)
  test_idx <- sample(nrow(d$F_hat), 50)
  train_idx <- setdiff(seq_len(nrow(d$F_hat)), test_idx)
  res <- split_conformal(d$F_hat[train_idx, , drop = FALSE],
                         d$W[train_idx, , drop = FALSE],
                         d$W[test_idx, , drop = FALSE], K = 3,
                         wrapper = mean, alpha = 0.1,
                         dfspa_args = list(denoise = FALSE), seed = 3)
  cov <- coverage(res, F_obs = d$F_hat[test_idx, , drop = FALSE])
  expect_gt(cov$overall, 0.6)
})
