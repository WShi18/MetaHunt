test_that("metahunt fits and predicts end-to-end", {
  set.seed(1)
  G <- 30; m <- 60
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  fit <- metahunt(F_hat, W, K = 3)

  expect_s3_class(fit, "metahunt")
  expect_equal(fit$m, m)
  expect_equal(fit$G_grid, G)
  expect_equal(fit$K, 3L)
  expect_equal(dim(fit$pi_hat), c(m, 3L))

  # predict with no wrapper returns matrix
  W_new <- W[1:5, , drop = FALSE]
  f_pred <- predict(fit, newdata = W_new)
  expect_equal(dim(f_pred), c(5L, G))

  # predict with wrapper returns vector
  s_pred <- predict(fit, newdata = W_new, wrapper = mean)
  expect_length(s_pred, 5L)
  expect_equal(s_pred, rowMeans(f_pred))
})

test_that("metahunt matches manual pipeline", {
  set.seed(2)
  G <- 25; m <- 50
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m))
  eta <- cbind(0.5 * W$w1, -0.3 * W$w1, rep(0, m))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.03), m, G)

  fit <- metahunt(F_hat, W, K = 3)

  manual_fit <- dfspa(F_hat, K = 3)
  manual_pi  <- project_to_simplex(F_hat, manual_fit$bases)
  manual_wm  <- fit_weight_model(manual_pi, W)
  manual_pred <- predict_target(manual_fit, manual_wm, W[1:3, , drop = FALSE])
  wrapped_pred <- predict(fit, newdata = W[1:3, , drop = FALSE])

  expect_equal(wrapped_pred, manual_pred)
})

test_that("metahunt fit integrates cleanly with conformal_from_fit", {
  set.seed(3)
  G <- 25; m <- 80
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m))
  eta <- cbind(0.5 * W$w1, -0.3 * W$w1, rep(0, m))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  tr <- 1:50; cal <- 51:70; new <- 71:80
  fit <- metahunt(F_hat[tr, ], W[tr, , drop = FALSE], K = 3,
                  dfspa_args = list(denoise = FALSE))

  res <- conformal_from_fit(
    fit$dfspa_fit, fit$weight_model,
    F_cal = F_hat[cal, , drop = FALSE],
    W_cal = W[cal, , drop = FALSE],
    W_new = W[new, , drop = FALSE],
    wrapper = mean
  )
  expect_s3_class(res, "metahunt_conformal")
  expect_length(res$prediction, length(new))
})

test_that("metahunt rejects invalid inputs", {
  G <- 20; m <- 30
  F_hat <- matrix(rnorm(m * G), m, G)
  W <- data.frame(w1 = rnorm(m))
  expect_error(metahunt(F_hat[-1, ], W, K = 3), "same number of rows")
  expect_error(metahunt(F_hat, W, K = 0), "positive integer")
  expect_error(metahunt(F_hat, W, K = 3, dfspa_args = "nope"), "must be named lists")
})

test_that("print.metahunt runs", {
  set.seed(4)
  F_hat <- matrix(rnorm(30 * 20), 30, 20)
  W <- data.frame(w1 = rnorm(30))
  fit <- metahunt(F_hat, W, K = 3, dfspa_args = list(denoise = FALSE))
  expect_output(print(fit), "MetaHunt fit")
})
