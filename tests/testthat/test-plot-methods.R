test_that("plot.metahunt runs without error", {
  set.seed(1)
  G <- 20; m <- 30
  F_hat <- matrix(rnorm(m * G), m, G)
  W <- data.frame(w1 = rnorm(m))
  fit <- metahunt(F_hat, W, K = 3, dfspa_args = list(denoise = FALSE))
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(fit))
  expect_invisible(plot(fit, x_axis = seq(0, 1, length.out = G)))
  expect_error(plot(fit, x_axis = 1:5), "length")
})

test_that("plot.metahunt_conformal runs for pointwise and scalar objects", {
  set.seed(2)
  G <- 25; m <- 100
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
  eta <- as.matrix(W) %*% cbind(c(1, -0.5), c(-0.4, 1), c(0, 0))
  pi_true <- exp(eta) / rowSums(exp(eta))
  F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

  res_pw <- split_conformal(F_hat, W, W[1:3, ], K = 3, seed = 1)
  res_sc <- split_conformal(F_hat, W, W[1:3, ], K = 3,
                            wrapper = mean, seed = 1)
  pdf(NULL)
  on.exit(dev.off())
  expect_invisible(plot(res_pw))
  expect_invisible(plot(res_pw, target_idx = 2))
  expect_invisible(plot(res_sc))
  expect_error(plot(res_pw, target_idx = 99), "out of range")
})
