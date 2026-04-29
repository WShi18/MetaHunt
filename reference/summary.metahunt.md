# Summarise a MetaHunt fit

Produces a compact summary of a `"metahunt"` object, including
study/grid sizes, the weight-model method, per-basis summary statistics
for the training simplex weights `pi_hat`, and denoising bookkeeping
from the underlying
[`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md) fit.

## Usage

``` r
# S3 method for class 'metahunt'
summary(object, ...)
```

## Arguments

- object:

  A `"metahunt"` object from
  [`metahunt()`](https://wshi18.github.io/MetaHunt/reference/metahunt.md).

- ...:

  Ignored.

## Value

An object of class `"summary.metahunt"`: a list with components

- `m`:

  Number of studies.

- `G_grid`:

  Grid size.

- `K`:

  Number of basis functions.

- `weight_method`:

  Method used by the weight model.

- `predictor_names`:

  Character vector of covariate names.

- `pi_summary`:

  A `K`-by-5 numeric matrix; each row gives `min`, `mean`, `median`,
  `max`, and `sd` of the corresponding column of `object$pi_hat`.

- `n_kept`:

  Number of studies retained after denoising.

- `n_dropped`:

  Number of studies dropped (`m - n_kept`).

- `denoising`:

  List with `N` and `Delta` from the `dfspa` fit (both `NA` when
  `denoise = FALSE`).

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
summary(fit)
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
