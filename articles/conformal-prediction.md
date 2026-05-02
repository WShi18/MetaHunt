# Conformal prediction with different choices

``` r

library(MetaHunt)
set.seed(1)
```

This vignette covers the three conformal-prediction interfaces exported
by MetaHunt. The validity of all three rests on the exchangeability
assumption A3 in
[`vignette("metahunt-intro", package = "MetaHunt")`](https://wshi18.github.io/MetaHunt/articles/metahunt-intro.md)
§“Key assumptions”; we do not re-derive it here.

## Why conformal prediction here

Conformal prediction wraps any black-box prediction rule and produces a
band around its forecast that, on average across new studies, will
contain the truth at least `(1 - alpha)` of the time. The key word is
*marginal*: the guarantee is over the random draw of the new study, not
conditional on a specific covariate value. All you need is for the
calibration data to be exchangeable with the new study (assumption A3) —
no distributional assumptions on the noise or on the weight model.

## A small standalone simulation

``` r

# m = 80 is large enough that with cal_frac = 0.5 and alpha = 0.05 the conformal quantile is finite.
m <- 80; G <- 20; K_true <- 3
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
beta <- cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
pi_true <- exp(as.matrix(W) %*% beta); pi_true <- pi_true / rowSums(pi_true)
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
```

``` r

W_new <- data.frame(w1 = c(0, 1, -1), w2 = c(0, -0.5, 1))
```

## The three flavours

All three return an object you can plot directly with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

| Function | When to use |
|----|----|
| [`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md) | Default. One train/calibration split. Fastest; some variance from the random split. |
| [`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md) | Many studies, want lower split-induced variance. Refits the pipeline `n_folds + 1` times. |
| [`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md) | You’ve already fit a pipeline (e.g. after tuning `K`) and want intervals without refitting. See [`?conformal_from_fit`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md). |

## Split conformal

[`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
does a single train / calibration split. With small `m`, set
`cal_frac = 0.5` so the calibration set is large enough for the chosen
`alpha`.

``` r

res_pw <- split_conformal(F_hat, W, W_new, K = K_true, alpha = 0.05,
                          cal_frac = 0.5, seed = 1,
                          dfspa_args = list(denoise = FALSE))
plot(res_pw, target_idx = 1, x_axis = x)
```

![](conformal-prediction_files/figure-html/cp-split-pointwise-1.png)

The shaded region is the pointwise 95% conformal band; it has finite
width because `n_cal = 40` is large enough for `α = 0.05`.

``` r

res_scalar <- split_conformal(F_hat, W, W_new, K = K_true,
                              wrapper = mean, alpha = 0.05,
                              cal_frac = 0.5, seed = 1,
                              dfspa_args = list(denoise = FALSE))
data.frame(prediction = res_scalar$prediction,
           lower      = res_scalar$lower,
           upper      = res_scalar$upper)
#>   prediction      lower     upper
#> 1 0.33312276 0.27629349 0.3899520
#> 2 0.50923298 0.45240371 0.5660622
#> 3 0.08545421 0.02862494 0.1422835
```

## Cross conformal

[`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md)
reduces the variance of the band that comes from the random split, at
the cost of refitting `n_folds + 1` times.

``` r

res_cross <- cross_conformal(F_hat, W, W_new, K = K_true, n_folds = 4,
                             wrapper = mean, alpha = 0.1, seed = 1,
                             dfspa_args = list(denoise = FALSE))
res_cross
#> MetaHunt conformal prediction
#>   method:        cross 
#>   alpha:         0.1 
#>   n calibration: 80 
#>   mode:          scalar (via wrapper)
#>   n targets:     3 
#>   quantile:      0.04072
```

## Pre-fit conformal

If you have already run
[`metahunt()`](https://wshi18.github.io/MetaHunt/reference/metahunt.md)
(for instance after tuning `K`) and do not want to refit,
[`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md)
recycles the existing fit to produce calibrated intervals. The example
below re-uses the training data as the calibration set *for
demonstration only*; in real use, hold out a separate calibration set so
the exchangeability argument applies to genuinely unseen studies.

``` r

fit <- metahunt(F_hat, W, K = K_true, dfspa_args = list(denoise = FALSE))
pi_hat <- project_to_simplex(F_hat, fit$dfspa_fit$bases)
res_pre <- conformal_from_fit(
  dfspa_fit = fit$dfspa_fit, weight_model = fit$weight_model,
  F_cal = F_hat, W_cal = W, W_new = W_new,
  wrapper = mean, alpha = 0.1
)
res_pre
#> MetaHunt conformal prediction
#>   method:        split 
#>   alpha:         0.1 
#>   n calibration: 80 
#>   mode:          scalar (via wrapper)
#>   n targets:     3 
#>   quantile:      0.033
```

## Pointwise vs scalar bands

A pointwise band returns a `(1 - alpha)` interval at each grid point but
does not give a joint guarantee across grid points: the probability that
the truth lies inside the entire band simultaneously is generally lower
than `1 - alpha`. A scalar wrapper (e.g. `wrapper = mean`) collapses the
function to a single number and gives one calibrated interval, which is
the right tool for joint inferential claims. If you need a joint
coverage statement across the grid, either apply a multiple-testing
correction (e.g. divide α by `G`) or replace the pointwise band with a
scalar wrapper — see
[`vignette('wrapper-scalar', package = 'MetaHunt')`](https://wshi18.github.io/MetaHunt/articles/wrapper-scalar.md).

## Small-`m` warning on `cal_frac`

> With too-few calibration studies for the chosen `alpha`, the conformal
> quantile is `Inf` and intervals are unbounded. The finite-sample
> formula needs `n_cal >= ceiling((1 - alpha)(n_cal + 1))` calibration
> studies; below that threshold the package warns and the bands
> degenerate. Either raise `cal_frac`, raise `alpha`, or switch to
> [`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md).

Below we deliberately reuse only the first 30 of our 80 studies so the
calibration set is too small for `α = 0.05`. The package issues a
warning and returns `Inf` quantiles; the corresponding intervals are
unbounded. The fix is to either supply more studies, raise `α`, or raise
`cal_frac`.

``` r

m_small <- 30  # too small for alpha = 0.05 with cal_frac = 0.5
F_small <- F_hat[1:m_small, , drop = FALSE]
W_small <- W[1:m_small, , drop = FALSE]
res_inf <- split_conformal(F_small, W_small, W_new, K = K_true,
                           alpha = 0.05, cal_frac = 0.5, seed = 1,
                           dfspa_args = list(denoise = FALSE))
#> Warning in .build_conformal_output(obs_cal = F_cal, pred_cal = pred_cal, : With
#> n_cal = 15 and alpha = 0.05, the conformal quantile is infinite at 20 of 20
#> grid points; intervals are unbounded there. Increase calibration size (raise
#> `cal_frac` or supply more studies) or use a larger `alpha`.
res_inf$quantile        # Inf — quantile is unbounded
#>  [1] Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf
#> [20] Inf
range(res_inf$lower)    # -Inf
#> [1] -Inf -Inf
range(res_inf$upper)    #  Inf
#> [1] Inf Inf
```

## See also

- [`vignette("metahunt-intro", package = "MetaHunt")`](https://wshi18.github.io/MetaHunt/articles/metahunt-intro.md)
  — pipeline context and the A3 exchangeability assumption.
- [`?split_conformal`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
  — single-split conformal calibration.
- [`?cross_conformal`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md)
  — cross-fitting conformal calibration.
- [`?conformal_from_fit`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md)
  — calibration using an existing fit.
- [`?coverage`](https://wshi18.github.io/MetaHunt/reference/coverage.md)
  — empirical coverage diagnostics for conformal bands.
- [`?plot.metahunt_conformal`](https://wshi18.github.io/MetaHunt/reference/plot.metahunt_conformal.md)
  — plotting method for the returned objects.
