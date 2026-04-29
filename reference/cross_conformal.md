# Cross-conformal prediction intervals (pooled K-fold scores)

Computes K-fold split conformal intervals in which the calibration
scores are pooled across folds, while point predictions for new targets
are produced from a final pipeline fit on all studies. Equivalent to
running
[`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
`n_folds` times with different calibration sets and pooling all
conformity scores into a single empirical distribution.

## Usage

``` r
cross_conformal(
  F_hat,
  W,
  W_new,
  K,
  alpha = 0.05,
  n_folds = 5L,
  wrapper = NULL,
  grid_weights = NULL,
  dfspa_args = list(),
  weight_model_args = list(),
  seed = NULL
)
```

## Arguments

- F_hat:

  An `m`-by-`G_grid` numeric matrix of study-level function evaluations.

- W:

  An `m`-by-`p` matrix or data frame of study-level covariates.

- W_new:

  A matrix or data frame of new target covariates. Must contain columns
  matching `W`.

- K:

  Integer number of basis functions.

- alpha:

  Miscoverage level; interval has nominal coverage \\1-\alpha\\. Default
  `0.05`.

- n_folds:

  Integer number of folds (\>= 2). Default `5`.

- wrapper:

  Optional reduction function (see
  [`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md)).
  If `NULL`, intervals are constructed pointwise at every grid point.

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector used for the
  \\L^2(\mu)\\ norm and for the default wrapper.

- dfspa_args, weight_model_args:

  Named lists passed to
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md) and
  [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md)
  respectively.

- seed:

  Optional integer seed for reproducible train/calibration splits.

## Value

An object of class `"metahunt_conformal"` (see
[`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
for fields). `method` is `"cross"` and `n_cal` is the number of pooled
scores.

## Details

For each fold, the MetaHunt pipeline is fit on the out-of-fold studies,
and conformity scores are computed on the in-fold studies. After all
folds complete, the pooled scores yield a single \\(1-\alpha)\\-quantile
(or one per grid point). Point predictions for `W_new` use a pipeline
refit on the full dataset. The interval at grid point `g` for target `j`
is \\\[\tilde f^{(j)}(x_g) - q_g, \tilde f^{(j)}(x_g) + q_g\]\\ (or the
scalar analogue when `wrapper` is supplied).

This differs from Vovk's original cross-conformal predictor for
classification. For regression, pooling scores across folds is a common
practical extension of split conformal and reduces the variance due to
the single calibration split. Exact finite-sample coverage is not
guaranteed; see Barber et al. (2021, Jackknife+) for more conservative
alternatives.

## Examples

``` r
# \donttest{
set.seed(1)
G <- 30; m <- 60; K_true <- 3
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m))
eta <- cbind(0.8 * W$w1, -0.3 * W$w1, rep(0, m))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
W_new <- data.frame(w1 = c(0, 1))

res <- cross_conformal(F_hat, W, W_new, K = 3, n_folds = 4,
                       dfspa_args = list(denoise = FALSE), seed = 1)
res
#> MetaHunt conformal prediction
#>   method:        cross 
#>   alpha:         0.05 
#>   n calibration: 60 
#>   mode:          pointwise (per grid point)
#>   n targets:     2 
#>   grid size:     30 
#>   quantile range: 0.1011 to 0.2328 
# }
```
