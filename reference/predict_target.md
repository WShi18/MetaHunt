# Predict the target function for new study-level covariates

Given a fitted d-fSPA basis decomposition and a fitted weight model,
compute the predicted target function on the shared grid as \$\$\tilde
f^{(0)}(\cdot) = \sum\_{k=1}^{\hat K} \tilde\pi\_{0k}\\ \hat g_k(\cdot),
\qquad \tilde\pi_0 = \widehat{\mathcal M}(\mathbf W_0).\$\$

## Usage

``` r
predict_target(dfspa_fit, weight_model, W_new)
```

## Arguments

- dfspa_fit:

  A `"dfspa"` object produced by
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md). Its
  `bases` slot (`K`-by-`G_grid` matrix) provides the basis functions on
  the grid.

- weight_model:

  A `"metahunt_weight_model"` object produced by
  [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md).
  Must have `weight_model$K == nrow(dfspa_fit$bases)`.

- W_new:

  A matrix or data frame of study-level covariates for the new target
  studies, with columns matching those used to fit `weight_model`.

## Value

An `nrow(W_new)`-by-`G_grid` numeric matrix; row `j` is the predicted
target function on the grid for the `j`-th new study.

## See also

[`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md)
to reduce predicted functions to scalars.

## Examples

``` r
set.seed(1)
G <- 40; m <- 60; K_true <- 3
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)

# generate study-level covariates and softmax weights
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
beta <- cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
eta  <- as.matrix(W) %*% beta
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.02), m, G)

fit    <- dfspa(F_hat, K = K_true)
pi_hat <- project_to_simplex(F_hat, fit$bases)
wm     <- fit_weight_model(pi_hat, W)

W_new   <- data.frame(w1 = c(0, 1), w2 = c(0, -1))
f_pred  <- predict_target(fit, wm, W_new)
dim(f_pred)   # 2 x G
#> [1]  2 40
```
