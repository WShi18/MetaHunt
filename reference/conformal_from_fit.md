# Split conformal intervals from a pre-fit MetaHunt pipeline

Lower-level entry point that builds split conformal intervals for new
target covariates using an already-fitted d-fSPA basis decomposition, an
already-fitted weight model, and a user-supplied calibration set. Use
this when you have independently tuned `K` or want to reuse a pipeline
fit; otherwise the high-level
[`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
is usually more convenient.

## Usage

``` r
conformal_from_fit(
  dfspa_fit,
  weight_model,
  F_cal,
  W_cal,
  W_new,
  alpha = 0.05,
  wrapper = NULL,
  grid_weights = NULL
)
```

## Arguments

- dfspa_fit:

  A `"dfspa"` object from
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md).

- weight_model:

  A `"metahunt_weight_model"` object from
  [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md).
  Must satisfy `weight_model$K == nrow(dfspa_fit$bases)`.

- F_cal:

  An `n_cal`-by-`G_grid` numeric matrix of observed study-level function
  evaluations for the calibration set. Calibration studies must be
  independent of both `dfspa_fit` and `weight_model` for valid coverage.

- W_cal:

  A matrix or data frame of study-level covariates for the calibration
  set, with the same columns used to fit `weight_model`.

- W_new:

  A matrix or data frame of study-level covariates for new target
  studies.

- alpha:

  Miscoverage level; default `0.05`.

- wrapper:

  Optional reduction function (see
  [`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md)).
  If `NULL`, intervals are pointwise at every grid point.

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector used for the
  default wrapper and for weighted inner products.

## Value

An object of class `"metahunt_conformal"`; see
[`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
for a description of its fields.

## See also

[`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
for the high-level version that splits and fits internally,
[`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md)
for the K-fold variant.

## Examples

``` r
set.seed(1)
G <- 30; m <- 80
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m))
eta <- cbind(0.8 * W$w1, -0.4 * W$w1, rep(0, m))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

# user-controlled split and fit
tr  <- 1:50; cal <- 51:70; new <- 71:80
fit <- dfspa(F_hat[tr, ], K = 3)
pih <- project_to_simplex(F_hat[tr, ], fit$bases)
wm  <- fit_weight_model(pih, W[tr, , drop = FALSE])

res <- conformal_from_fit(
  fit, wm,
  F_cal = F_hat[cal, ], W_cal = W[cal, , drop = FALSE],
  W_new = W[new, , drop = FALSE],
  wrapper = mean
)
res
#> MetaHunt conformal prediction
#>   method:        split 
#>   alpha:         0.05 
#>   n calibration: 20 
#>   mode:          scalar (via wrapper)
#>   n targets:     10 
#>   quantile:      0.1492 
```
