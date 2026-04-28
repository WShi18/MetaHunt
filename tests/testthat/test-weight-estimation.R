test_that("project_to_simplex recovers exact weights in the noiseless case", {
  set.seed(1)
  G <- 40
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x), x)
  true_pi <- rbind(
    diag(3),
    c(0.5, 0.3, 0.2),
    c(0.1, 0.7, 0.2),
    c(0.25, 0.25, 0.5)
  )
  F_hat <- true_pi %*% basis

  pi_hat <- project_to_simplex(F_hat, basis)

  expect_equal(dim(pi_hat), c(6L, 3L))
  expect_true(all(pi_hat >= -1e-10))
  expect_equal(rowSums(pi_hat), rep(1, 6), tolerance = 1e-8)
  expect_equal(pi_hat, true_pi, tolerance = 1e-6)
})

test_that("project_to_simplex clips external points to the simplex boundary", {
  G <- 30
  x <- seq(0, 1, length.out = G)
  basis <- rbind(sin(pi * x), cos(pi * x))
  # a function living outside the convex hull: 1.5 * g1
  f_ext <- 1.5 * basis[1, , drop = FALSE]
  pi_hat <- project_to_simplex(f_ext, basis)
  # closest point on the simplex is (1, 0) (pure g1)
  expect_equal(as.numeric(pi_hat), c(1, 0), tolerance = 1e-6)
})

test_that("project_to_simplex handles K = 1", {
  G <- 20
  F_hat <- matrix(rnorm(4 * G), 4, G)
  basis <- matrix(rnorm(G), 1, G)
  pi_hat <- project_to_simplex(F_hat, basis)
  expect_equal(pi_hat, matrix(1, 4, 1))
})

test_that("project_to_simplex rejects mismatched grids", {
  expect_error(
    project_to_simplex(matrix(0, 3, 10), matrix(0, 2, 12)),
    "same number of columns"
  )
})

test_that("fit_weight_model + predict give simplex outputs close to truth", {
  set.seed(123)
  m <- 150; K <- 3; p <- 2
  W <- matrix(rnorm(m * p), m, p)
  colnames(W) <- c("w1", "w2")
  # weights driven by W via softmax of linear predictor
  beta <- cbind(c(1.0, -0.8), c(-0.5, 1.2), c(0, 0))
  eta <- W %*% beta
  pi_true <- exp(eta) / rowSums(exp(eta))
  # add small simplex noise via a Dirichlet draw centered at pi_true
  pi_hat <- t(apply(pi_true, 1, function(p) {
    a <- 80 * p + 1e-3
    g <- stats::rgamma(K, shape = a)
    g / sum(g)
  }))

  model <- fit_weight_model(pi_hat, W)
  expect_s3_class(model, "metahunt_weight_model")
  expect_equal(model$K, K)

  # predict on held-out points
  W_new <- matrix(c(0, 0, 1, -1, -1, 1), ncol = 2, byrow = TRUE)
  colnames(W_new) <- c("w1", "w2")
  pi_new <- predict(model, newdata = W_new)
  expect_equal(dim(pi_new), c(3L, K))
  expect_equal(rowSums(pi_new), rep(1, 3), tolerance = 1e-6)
  expect_true(all(pi_new >= 0))

  # compare to the true softmax weights at the same points
  eta_new <- W_new %*% beta
  pi_true_new <- exp(eta_new) / rowSums(exp(eta_new))
  expect_true(max(abs(pi_new - pi_true_new)) < 0.15)
})

test_that("fit_weight_model handles boundary weights via shrinkage", {
  set.seed(7)
  m <- 60
  W <- data.frame(w1 = rnorm(m))
  # include pure vertex rows (boundary!)
  pi_hat <- rbind(
    matrix(c(1, 0, 0), 5, 3, byrow = TRUE),
    matrix(c(0, 1, 0), 5, 3, byrow = TRUE),
    matrix(c(0, 0, 1), 5, 3, byrow = TRUE)
  )
  # fill the rest with interior Dirichlet draws
  rest <- t(replicate(m - 15, {
    g <- stats::rgamma(3, shape = 1); g / sum(g)
  }))
  pi_hat <- rbind(pi_hat, rest)
  # DirichletReg may emit harmless convergence warnings on some BLAS
  # implementations (e.g. Windows); we only require that the fit completes
  # and returns a model of the right class.
  suppressWarnings(model <- fit_weight_model(pi_hat, W))
  expect_s3_class(model, "metahunt_weight_model")
})

test_that("fit_weight_model rejects malformed pi_hat", {
  W <- data.frame(w1 = rnorm(5))
  bad_pi <- matrix(c(0.5, 0.2, 0.1), 5, 3, byrow = TRUE)  # rows do not sum to 1
  expect_error(fit_weight_model(bad_pi, W), "sum to 1")
  neg_pi <- matrix(c(-0.1, 0.5, 0.6), 5, 3, byrow = TRUE)
  expect_error(fit_weight_model(neg_pi, W), "non-negative")
})
