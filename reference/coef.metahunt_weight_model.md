# Extract coefficients from a MetaHunt weight model

Returns the regression coefficients from the underlying weight-model
fit. For the default `"dirichlet"` method this delegates to
[`DirichletReg::DirichReg()`](https://rdrr.io/pkg/DirichletReg/man/dirichreg.html)'s
[`coef()`](https://rdrr.io/r/stats/coef.html) method.

## Usage

``` r
# S3 method for class 'metahunt_weight_model'
coef(object, ...)
```

## Arguments

- object:

  A fitted `"metahunt_weight_model"` from
  [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md).

- ...:

  Passed through to
  [`stats::coef()`](https://rdrr.io/r/stats/coef.html).

## Value

The coefficient vector / matrix returned by
[`DirichletReg::DirichReg()`](https://rdrr.io/pkg/DirichletReg/man/dirichreg.html)'s
[`coef()`](https://rdrr.io/r/stats/coef.html) method (numeric vector or
matrix depending on the parametrisation used by the underlying fit).

## Examples

``` r
# \donttest{
set.seed(1)
m <- 60; K <- 3
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
eta <- cbind(0.5 * W$w1, -0.3 * W$w2, rep(0, m))
pi_true <- exp(eta) / rowSums(exp(eta))
pi_hat <- pi_true + matrix(rnorm(m * K, sd = 0.01), m, K)
pi_hat <- pmax(pi_hat, 0); pi_hat <- pi_hat / rowSums(pi_hat)
model <- fit_weight_model(pi_hat, W)
coef(model)
#> $v1
#> (Intercept)          w1          w2 
#>   6.7709026   0.7037691   0.1620994 
#> 
#> $v2
#> (Intercept)          w1          w2 
#>   6.7695064   0.1975436  -0.1387985 
#> 
#> $v3
#> (Intercept)          w1          w2 
#>   6.7736531   0.2093901   0.1684880 
#> 
# }
```
