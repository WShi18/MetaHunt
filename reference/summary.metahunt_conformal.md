# Summarise a conformal prediction-interval object

Produces a small list of descriptive statistics about a
`"metahunt_conformal"` object: interval widths, quantile summaries, and
calibration diagnostics. Returns an object of class
`"summary.metahunt_conformal"` with a matching `print` method.

## Usage

``` r
# S3 method for class 'metahunt_conformal'
summary(object, ...)
```

## Arguments

- object:

  A `"metahunt_conformal"` object from
  [`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md),
  [`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md),
  or
  [`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md).

- ...:

  Unused; present for S3 generic consistency.

## Value

A list of class `"summary.metahunt_conformal"`. In pointwise mode (no
wrapper) the list contains `n_targets`, `G_grid`, `n_cal`, `alpha`,
`method`, `mean_interval_width`, `frac_finite_quantile`,
`quantile_summary`, and `wrapper`. In scalar mode (wrapper supplied) the
list contains `n_targets`, `n_cal`, `alpha`, `method`,
`mean_interval_width`, `quantile`, `quantile_finite`, and `wrapper`.

## Examples

``` r
set.seed(1)
G <- 25; m <- 60
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m))
eta <- cbind(0.6 * W$w1, -0.3 * W$w1, rep(0, m))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
W_new <- data.frame(w1 = c(0, 1))
res <- split_conformal(F_hat, W, W_new, K = 3,
                       dfspa_args = list(denoise = FALSE), seed = 1)
#> Warning: With n_cal = 18 and alpha = 0.05, the conformal quantile is infinite at 25 of 25 grid points; intervals are unbounded there. Increase calibration size (raise `cal_frac` or supply more studies) or use a larger `alpha`.
summary(res)
#> Summary of MetaHunt conformal prediction
#>   method:         split 
#>   alpha:          0.05 
#>   n calibration:  18 
#>   mode:           pointwise (per grid point)
#>   n targets:      2 
#>   grid size:      25 
#>   mean interval width: NA 
#>   fraction finite quantiles: 0 
#>   quantile:       all infinite (insufficient calibration)
```
