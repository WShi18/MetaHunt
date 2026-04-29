# Predict target functions (or scalar summaries) from a MetaHunt fit

Predict target functions (or scalar summaries) from a MetaHunt fit

## Usage

``` r
# S3 method for class 'metahunt'
predict(object, newdata, wrapper = NULL, grid_weights = NULL, ...)
```

## Arguments

- object:

  A `"metahunt"` object from
  [`metahunt()`](https://wshi18.github.io/MetaHunt/reference/metahunt.md).

- newdata:

  A matrix or data frame of new study-level covariates.

- wrapper:

  Optional reduction function. If `NULL`, returns the full predicted
  function matrix (`nrow(newdata)`-by-`G_grid`). If a function, applied
  to each predicted function to return a scalar per new target; see
  [`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md).

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector used only when
  `wrapper = NULL` with the default weighted-mean path. Defaults to the
  `grid_weights` stored in `object`.

- ...:

  Ignored.

## Value

Either an `nrow(newdata)`-by-`G_grid` matrix of predicted functions
(when `wrapper = NULL`) or a length-`nrow(newdata)` numeric vector of
scalar summaries.

## Examples

``` r
set.seed(1)
G <- 25; m <- 40
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

fit <- metahunt(F_hat, W, K = 3, dfspa_args = list(denoise = FALSE))
f_pred <- predict(fit, newdata = W[1:3, ])
dim(f_pred)                                       # 3 x G
#> [1]  3 25
predict(fit, newdata = W[1:3, ], wrapper = mean)  # scalar summaries
#> [1] 0.3084921 0.4289405 0.1178668
```
