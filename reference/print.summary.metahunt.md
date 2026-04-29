# Print a `summary.metahunt` object

Print a `summary.metahunt` object

## Usage

``` r
# S3 method for class 'summary.metahunt'
print(x, ...)
```

## Arguments

- x:

  A `"summary.metahunt"` object from
  [`summary.metahunt()`](https://wshi18.github.io/MetaHunt/reference/summary.metahunt.md).

- ...:

  Ignored.

## Value

Invisibly returns `x`.

## Examples

``` r
set.seed(1)
G <- 25; m <- 40
x_grid <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x_grid), cos(pi * x_grid), x_grid)
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

fit <- metahunt(F_hat, W, K = 3, dfspa_args = list(denoise = FALSE))
print(summary(fit))
#> MetaHunt fit summary
#>   m (studies):    40 
#>   G (grid size):  25 
#>   K (bases):      3 
#>   weight method:  dirichlet 
#>   predictors:     w1, w2 
#>   studies kept:   40 
#>   studies dropped: 0 
#>   denoising N:    NA 
#>   denoising Delta: NA 
#> 
#> Per-basis pi_hat summary:
#>         min   mean median max     sd
#> basis_1   0 0.3099 0.2376   1 0.2845
#> basis_2   0 0.3117 0.1956   1 0.3204
#> basis_3   0 0.3784 0.3983   1 0.2214
```
