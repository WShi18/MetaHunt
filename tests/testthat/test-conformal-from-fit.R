make_toy <- function(seed = 1, m = 100, G = 30) {
  set.seed(seed)
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
  list(F_hat = F_hat, W = W)
}

test_that("conformal_from_fit matches split_conformal with the same split", {
  d <- make_toy(m = 100)
  m <- nrow(d$F_hat)
  # deterministic split
  cal_idx <- 1:30
  tr_idx  <- setdiff(seq_len(m), cal_idx)
  W_new <- d$W[1:5, , drop = FALSE]

  # via high-level split_conformal with explicit calibration_idx
  a <- split_conformal(d$F_hat, d$W, W_new, K = 3,
                       calibration_idx = cal_idx, wrapper = mean)

  # via low-level conformal_from_fit with manual fit
  fit <- dfspa(d$F_hat[tr_idx, , drop = FALSE], K = 3)
  pih <- project_to_simplex(d$F_hat[tr_idx, , drop = FALSE], fit$bases)
  wm  <- fit_weight_model(pih, d$W[tr_idx, , drop = FALSE])
  b <- conformal_from_fit(fit, wm,
                          F_cal = d$F_hat[cal_idx, , drop = FALSE],
                          W_cal = d$W[cal_idx, , drop = FALSE],
                          W_new = W_new,
                          wrapper = mean)

  expect_equal(a$prediction, b$prediction)
  expect_equal(a$lower, b$lower)
  expect_equal(a$upper, b$upper)
  expect_equal(a$quantile, b$quantile)
})

test_that("conformal_from_fit pointwise: shapes", {
  d <- make_toy(m = 60)
  tr  <- 1:40; cal <- 41:55; new <- 56:60

  fit <- dfspa(d$F_hat[tr, ], K = 3)
  pih <- project_to_simplex(d$F_hat[tr, ], fit$bases)
  wm  <- fit_weight_model(pih, d$W[tr, , drop = FALSE])

  res <- suppressWarnings(conformal_from_fit(
    fit, wm,
    F_cal = d$F_hat[cal, , drop = FALSE],
    W_cal = d$W[cal, , drop = FALSE],
    W_new = d$W[new, , drop = FALSE]
  ))
  expect_equal(dim(res$prediction), c(length(new), ncol(d$F_hat)))
  expect_length(res$quantile, ncol(d$F_hat))
  expect_true(all(res$lower <= res$prediction))
  expect_true(all(res$upper >= res$prediction))
})

test_that("conformal_from_fit rejects basis/weight-model K mismatch", {
  d <- make_toy(m = 50)
  fit2 <- dfspa(d$F_hat[1:30, ], K = 2)
  fit3 <- dfspa(d$F_hat[1:30, ], K = 3)
  pih3 <- project_to_simplex(d$F_hat[1:30, ], fit3$bases)
  wm3  <- fit_weight_model(pih3, d$W[1:30, , drop = FALSE])
  expect_error(
    conformal_from_fit(fit2, wm3,
                       F_cal = d$F_hat[31:50, ], W_cal = d$W[31:50, , drop = FALSE],
                       W_new = d$W[1:2, , drop = FALSE]),
    "mismatch"
  )
})
