# Fit the full MetaHunt pipeline

End-to-end convenience wrapper that runs the three training-time steps
of the MetaHunt pipeline in sequence: (1)
[`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md) for
basis hunting, (2)
[`project_to_simplex()`](https://wshi18.github.io/MetaHunt/reference/project_to_simplex.md)
for per-study weight recovery, and (3)
[`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md)
for modelling the weight-to-covariate map. The result supports
[`predict.metahunt()`](https://wshi18.github.io/MetaHunt/reference/predict.metahunt.md)
for generating target-function predictions on new study-level
covariates.

## Usage

``` r
metahunt(
  F_hat,
  W,
  K,
  grid_weights = NULL,
  dfspa_args = list(),
  weight_model_args = list()
)
```

## Arguments

- F_hat:

  An `m`-by-`G_grid` numeric matrix of study-level function evaluations
  on a shared grid; row `i` is \\\hat f^{(i)}\\.

- W:

  An `m`-by-`p` matrix or data frame of study-level covariates.

- K:

  Integer number of basis functions.

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector defining the
  \\L^2(\mu)\\ inner product; defaults to uniform.

- dfspa_args:

  Named list of extra arguments for
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md).

- weight_model_args:

  Named list of extra arguments for
  [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md).

## Value

An object of class `"metahunt"`: a list with the `dfspa_fit`,
`weight_model`, training `pi_hat`, `K`, and a stored copy of
`grid_weights`.

## Details

For uncertainty quantification, pair a `"metahunt"` fit with
[`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md)
(requires a separate calibration set). The high-level
[`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
and
[`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md)
functions perform their own fitting and do not consume a pre-fit
`"metahunt"` object.

## See also

[`predict.metahunt()`](https://wshi18.github.io/MetaHunt/reference/predict.metahunt.md),
[`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md),
[`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md),
[`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md),
[`reconstruction_error_curve()`](https://wshi18.github.io/MetaHunt/reference/reconstruction_error_curve.md),
[`cv_error_curve()`](https://wshi18.github.io/MetaHunt/reference/cv_error_curve.md).

## Examples

``` r
set.seed(1)
G <- 40; m <- 80
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

fit <- metahunt(F_hat, W, K = 3)
fit
#> MetaHunt fit
#>   m (studies):    80 
#>   G (grid size):  40 
#>   K (bases):      3 
#>   weight method:  dirichlet 
#>   predictors:     w1, w2 
f_pred <- predict(fit, newdata = W[1:3, ])
dim(f_pred)                                  # 3 x G: predicted functions
#> [1]  3 40
predict(fit, newdata = W[1:3, ], wrapper = mean)  # scalar summaries
#> [1] 0.39676714 0.38963182 0.08180082
```
