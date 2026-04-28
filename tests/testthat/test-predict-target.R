test_that("predict_target reproduces the linear combination pi %*% bases", {
  set.seed(1)
  G <- 40; m <- 80; K <- 3
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  beta <- cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
  eta  <- as.matrix(W) %*% beta
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.01), m, G)

  fit     <- dfspa(F_hat, K = K)
  pi_hat  <- project_to_simplex(F_hat, fit$bases)
  wm      <- fit_weight_model(pi_hat, W)

  W_new   <- data.frame(w1 = c(0, 0.5), w2 = c(0.2, -0.3))
  f_pred  <- predict_target(fit, wm, W_new)
  pi_pred <- predict(wm, newdata = W_new)

  expect_equal(dim(f_pred), c(2L, G))
  expect_equal(f_pred, pi_pred %*% fit$bases)
})

test_that("predict_target errors on mismatched K between fit and weight model", {
  set.seed(2)
  G <- 20; m <- 40
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m))
  pi_true <- matrix(rgamma(m * 3, 1), m, 3); pi_true <- pi_true / rowSums(pi_true)
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.01), m, G)

  fit_K2 <- dfspa(F_hat, K = 2)
  pi_hat <- project_to_simplex(F_hat, dfspa(F_hat, K = 3)$bases)
  wm_K3  <- fit_weight_model(pi_hat, W)
  expect_error(predict_target(fit_K2, wm_K3, W), "mismatch")
})

test_that("apply_wrapper default equals weighted mean", {
  F_mat <- matrix(c(1, 2, 3, 4,
                    10, 20, 30, 40), nrow = 2, byrow = TRUE)
  expect_equal(apply_wrapper(F_mat), c(mean(1:4), mean(c(10, 20, 30, 40))))
  w <- c(0.1, 0.2, 0.3, 0.4)
  expect_equal(apply_wrapper(F_mat, grid_weights = w),
               as.numeric(F_mat %*% w) / sum(w))
})

test_that("apply_wrapper supports arbitrary wrappers", {
  F_mat <- matrix(c(1, 2, 3,
                    7, 5, 6), nrow = 2, byrow = TRUE)
  expect_equal(apply_wrapper(F_mat, wrapper = max), c(3, 7))
  expect_equal(apply_wrapper(F_mat, wrapper = function(f) f[2]), c(2, 5))
})

test_that("apply_wrapper errors on bad wrapper output", {
  F_mat <- matrix(1:6, 2, 3)
  bad_wrapper <- function(f) c(1, 2)
  expect_error(apply_wrapper(F_mat, wrapper = bad_wrapper), "single numeric")
})
