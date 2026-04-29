# Predict simplex weights for new study-level covariates

Predict simplex weights for new study-level covariates

## Usage

``` r
# S3 method for class 'metahunt_weight_model'
predict(object, newdata, ...)
```

## Arguments

- object:

  A fitted `"metahunt_weight_model"` from
  [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md).

- newdata:

  A matrix or data frame of new study-level covariates with the same
  columns used at fitting.

- ...:

  Ignored.

## Value

An `nrow(newdata)`-by-`K` numeric matrix of predicted simplex weights
(component means); rows sum to 1.

## Examples

``` r
# \donttest{
set.seed(1)
m <- 40; K <- 3
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
eta <- cbind(0.5 * W$w1, -0.3 * W$w2, rep(0, m))
pi_true <- exp(eta) / rowSums(exp(eta))
pi_hat <- pi_true + matrix(rnorm(m * K, sd = 0.01), m, K)
pi_hat <- pmax(pi_hat, 0); pi_hat <- pi_hat / rowSums(pi_hat)
model <- fit_weight_model(pi_hat, W)
predict(model, newdata = data.frame(w1 = c(0, 1), w2 = c(0, -1)))
#>           [,1]      [,2]      [,3]
#> [1,] 0.3348964 0.3295546 0.3355490
#> [2,] 0.4044134 0.3385114 0.2570752
# }
```
