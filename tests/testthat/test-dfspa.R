test_that("dfspa recovers pure vertices in the noiseless case", {
  set.seed(1)
  G <- 50
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  pi_mat <- rbind(
    diag(3),
    c(0.5, 0.3, 0.2),
    c(0.2, 0.5, 0.3),
    c(0.3, 0.3, 0.4),
    c(0.1, 0.8, 0.1),
    c(0.1, 0.1, 0.8)
  )
  F_hat <- pi_mat %*% basis

  fit <- dfspa(F_hat, K = 3, denoise = FALSE)

  expect_s3_class(fit, "dfspa")
  expect_length(fit$selected, 3L)
  expect_setequal(fit$original_indices, c(1L, 2L, 3L))
  expect_equal(dim(fit$bases), c(3L, G))
})

test_that("dfspa with denoising still recovers bases under small noise", {
  set.seed(42)
  G <- 80
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x^2)
  K_true <- 3L
  m <- 40L
  # generate random simplex weights concentrated near vertices
  pi_mat <- matrix(0, m, K_true)
  pi_mat[1:3, ] <- diag(3)
  pi_mat[4:6, ] <- diag(3)
  for (i in 7:m) {
    w <- stats::rgamma(K_true, shape = 0.3)
    pi_mat[i, ] <- w / sum(w)
  }
  F_true <- pi_mat %*% basis
  F_hat  <- F_true + matrix(stats::rnorm(m * G, sd = 0.02), m, G)

  fit <- dfspa(F_hat, K = K_true)

  expect_equal(fit$K, K_true)
  expect_equal(length(fit$original_indices), K_true)
  # recovered bases should match true bases (up to permutation) within
  # tolerance set by noise and purity level
  dists <- outer(
    seq_len(K_true), seq_len(K_true),
    Vectorize(function(i, j) sqrt(mean((fit$bases[i, ] - basis[j, ])^2)))
  )
  best_match <- apply(dists, 1, min)
  expect_true(all(best_match < 0.05))
})

test_that("dfspa errors on invalid inputs", {
  F_hat <- matrix(1:20, nrow = 5, ncol = 4)
  expect_error(dfspa(as.data.frame(F_hat), K = 2), "numeric matrix")
  expect_error(dfspa(F_hat, K = 0), "positive integer")
  expect_error(dfspa(F_hat, K = 10), "must not exceed")
  expect_error(dfspa(F_hat, K = 2, grid_weights = c(1, 2, 3)),
               "length-`ncol\\(F_hat\\)`")
  expect_error(dfspa(F_hat, K = 2, grid_weights = c(-1, 1, 1, 1)),
               "non-negative")
})

test_that("grid_weights reshape the norm and drive first-iteration selection", {
  G <- 40
  x <- seq(0, 1, length.out = G)
  f1 <- ifelse(x < 0.5, 2, 0)        # mass on left half
  f2 <- ifelse(x < 0.5, 0, 2)        # mass on right half
  f3 <- rep(0.1, G)                   # small everywhere
  F_hat <- rbind(f1, f2, f3)

  # uniform: tied norms for f1, f2; which.max breaks the tie by index = 1
  fit_u <- dfspa(F_hat, K = 1, denoise = FALSE)
  expect_equal(fit_u$original_indices, 1L)

  # weights concentrated on the right half should make f2 win
  w_right <- ifelse(x < 0.5, 0.01, 1)
  fit_w <- dfspa(F_hat, K = 1, denoise = FALSE, grid_weights = w_right)
  expect_equal(fit_w$original_indices, 2L)
})
