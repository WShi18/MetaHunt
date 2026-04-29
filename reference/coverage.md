# Empirical coverage of a conformal prediction-interval object

Computes empirical coverage indicators of a fitted
`"metahunt_conformal"` object against held-out observed study-level
functions. The held-out studies must correspond positionally to the
targets used to build `object` (i.e. `F_obs[i, ]` is the observed
function for the same target whose prediction is in
`object$prediction[i, ]` or `object$prediction[i]`).

## Usage

``` r
coverage(object, F_obs, grid_weights = NULL)
```

## Arguments

- object:

  A `"metahunt_conformal"` object from
  [`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md),
  [`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md),
  or
  [`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md).

- F_obs:

  An `n_target`-by-`G_grid` numeric matrix of observed study-level
  functions for the target studies, in the same row order as
  `object$prediction`.

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector. Only used in
  scalar mode (passed to
  [`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md)).

## Value

A list. In **pointwise** mode the list contains

- `pointwise`:

  `n_target`-by-`G_grid` logical matrix of coverage indicators.

- `per_target`:

  Length-`n_target` numeric vector of mean coverage across the grid for
  each target.

- `per_grid_point`:

  Length-`G_grid` numeric vector of mean coverage across targets at each
  grid point.

- `overall`:

  Scalar mean coverage across all entries.

- `nominal`:

  Nominal coverage `1 - object$alpha`.

In **scalar** mode the list contains

- `pointwise`:

  Length-`n_target` logical vector of coverage indicators.

- `overall`:

  Scalar mean coverage.

- `nominal`:

  Nominal coverage `1 - object$alpha`.

## Details

In **pointwise** mode each entry `(i, g)` of `F_obs` is compared against
`[object$lower[i, g], object$upper[i, g]]`. In **scalar** mode `F_obs`
is first reduced to a length-`n_target` vector via
[`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md)
using `object$wrapper` and the supplied `grid_weights`, and then
compared against `[object$lower, object$upper]`.

Coverage is a finite-sample diagnostic. Nominal coverage is
`1 - object$alpha`; empirical coverage will fluctuate around this value
due to sampling.

## Examples

``` r
set.seed(1)
G <- 25; m <- 80
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m))
eta <- cbind(0.6 * W$w1, -0.3 * W$w1, rep(0, m))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

# held-out test set: same data-generating process, same W
test_idx <- 1:10
train_idx <- setdiff(seq_len(m), test_idx)
res <- split_conformal(F_hat[train_idx, ], W[train_idx, , drop = FALSE],
                       W[test_idx, , drop = FALSE], K = 3,
                       dfspa_args = list(denoise = FALSE), seed = 1)
cov <- coverage(res, F_obs = F_hat[test_idx, , drop = FALSE])
cov$overall
#> [1] 0.964
```
