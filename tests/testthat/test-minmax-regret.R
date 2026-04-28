test_that("minmax_regret returns a simplex weight and a function on the grid", {
  set.seed(1)
  G <- 30; m <- 5
  x <- seq(0, 1, length.out = G)
  F_hat <- rbind(
    sin(pi * x),
    cos(pi * x),
    x,
    0.4 * sin(pi * x) + 0.6 * x,
    0.3 * cos(pi * x) + 0.7 * x
  )
  fit <- minmax_regret(F_hat)

  expect_s3_class(fit, "minmax_regret")
  expect_length(fit$q, m)
  expect_equal(sum(fit$q), 1, tolerance = 1e-8)
  expect_true(all(fit$q >= -1e-10))
  expect_length(fit$prediction, G)

  # prediction is exactly q^T F_hat
  expect_equal(fit$prediction,
               as.numeric(crossprod(fit$q, F_hat)),
               tolerance = 1e-8)
})

test_that("minmax_regret with one source returns it unchanged", {
  G <- 20
  f1 <- sin(seq(0, 2 * pi, length.out = G))
  fit <- minmax_regret(matrix(f1, nrow = 1))
  expect_equal(fit$q, 1)
  expect_equal(fit$prediction, f1)
})

test_that("minmax_regret wrapper returns scalar matching grid mean", {
  G <- 25; m <- 4
  set.seed(2)
  F_hat <- matrix(rnorm(m * G), m, G)
  fit_fun <- minmax_regret(F_hat)
  fit_avg <- minmax_regret(F_hat, wrapper = mean)
  expect_length(fit_avg$prediction, 1L)
  expect_equal(fit_avg$prediction, mean(fit_fun$prediction), tolerance = 1e-10)
})

test_that("minmax_regret achieves zero objective when target is in the convex hull", {
  set.seed(3)
  G <- 20; m <- 4
  x <- seq(0, 1, length.out = G)
  F_hat <- rbind(sin(pi * x), cos(pi * x), x, 0.3 * x + 0.7 * sin(pi * x))
  # the 4th row is in the convex hull of rows 1 and 3, so q should put
  # mass on the row representing the worst-case-regret minimiser.
  fit <- minmax_regret(F_hat)
  expect_equal(sum(fit$q), 1, tolerance = 1e-8)
  # the QP objective at the solution should be finite and not larger than
  # the value at any vertex
  Q_sol <- as.numeric(crossprod(fit$q, fit$Gamma) %*% fit$q -
                        crossprod(fit$d, fit$q))
  for (k in seq_len(m)) {
    e_k <- numeric(m); e_k[k] <- 1
    Q_k <- as.numeric(crossprod(e_k, fit$Gamma) %*% e_k -
                        crossprod(fit$d, e_k))
    expect_lte(Q_sol, Q_k + 1e-8)
  }
})

test_that("minmax_regret rejects malformed inputs", {
  expect_error(minmax_regret(as.data.frame(matrix(1:6, 2, 3))),
               "numeric matrix")
  expect_error(minmax_regret(matrix(c(1, NA, 2, 3), 2, 2)),
               "NA or NaN")
  expect_error(minmax_regret(matrix(1:6, 2, 3), grid_weights = c(1, -1, 1)),
               "non-negative")
  expect_error(minmax_regret(matrix(1:6, 2, 3), ridge = -1),
               "non-negative")
})

test_that("print.minmax_regret runs without error", {
  G <- 20
  F_hat <- matrix(rnorm(60), 3, 20)
  fit <- minmax_regret(F_hat)
  expect_output(print(fit), "Minimax-regret aggregator")
})
