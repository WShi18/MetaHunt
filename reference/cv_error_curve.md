# Cross-validated prediction-error curve for basis-rank selection

For each candidate `K`, perform k-fold cross-validation at the study
level. Within each fold, the full MetaHunt pipeline (d-fSPA +
constrained projection + weight model) is refit on the training studies,
and the held-out studies' functions are predicted from their study-level
covariates. The prediction error is \\\\\hat f^{(i)} - \tilde
f^{(i)}\\\_{L^2(\mu)}\\ averaged over held-out studies and then over
folds.

## Usage

``` r
cv_error_curve(
  F_hat,
  W,
  K_range = NULL,
  n_folds = 5L,
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

- K_range:

  Integer vector of candidate `K` values. Defaults to
  `2:min(nrow(F_hat) - 1, 10)`.

- n_folds:

  Integer number of CV folds (default `5`).

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector.

- dfspa_args:

  Named list of extra arguments for
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md).

- weight_model_args:

  Named list of extra arguments for
  [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md).

- seed:

  Optional integer seed for reproducible fold assignment; if `NULL` no
  seeding is performed.

## Value

A data frame with columns `K`, `cv_error` (mean over folds), and `cv_se`
(standard error across folds). The per-fold error matrix is attached as
the attribute `"fold_errors"` (`length(K_range)`-by-`n_folds`). Folds
where the pipeline fails contribute `NA` and are summarised in a single
warning.

## Details

This is the supervised rank-selection criterion of Section 3.2 of the
paper. Each held-out study is excluded from both basis hunting and
weight-model fitting.

## Examples

``` r
# \donttest{
set.seed(1)
G <- 40; m <- 80; K_true <- 3
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.02), m, G)
cv <- cv_error_curve(F_hat, W, K_range = 2:5, n_folds = 4, seed = 1)
cv
#>   K   cv_error       cv_se
#> 1 2 0.08336472 0.003590282
#> 2 3 0.06120156 0.007567407
#> 3 4 0.10205196 0.012299971
#> 4 5 0.11146587 0.012849123
# }
```
