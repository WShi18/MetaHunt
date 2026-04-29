# Reduce predicted functions to scalars via a user-supplied wrapper

Many downstream quantities of interest (average treatment effect,
pointwise predictions, or other functionals of \\f^{(0)}\\) are scalar
summaries of the predicted function. `apply_wrapper()` applies any
user-supplied reduction to each row of a function matrix, with a default
of the weighted mean with respect to `grid_weights` (which coincides
with \\\int f\\d\mu\\ when `grid_weights` represents \\\mu\\).

## Usage

``` r
apply_wrapper(F_mat, wrapper = NULL, grid_weights = NULL)
```

## Arguments

- F_mat:

  An `n`-by-`G_grid` numeric matrix; each row a function on a shared
  grid (e.g. the output of
  [`predict_target()`](https://wshi18.github.io/MetaHunt/reference/predict_target.md)).

- wrapper:

  Either `NULL` (default, weighted mean) or a function that takes a
  single numeric vector of length `G_grid` and returns a scalar.
  Examples: `mean`, `median`, `max`, `function(f) f[17]` (point
  evaluation), `function(f) sum(f^2)`.

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector. Used only when
  `wrapper = NULL`. Defaults to uniform `1 / G_grid`.

## Value

A length-`n` numeric vector of scalar summaries.

## Examples

``` r
F <- matrix(1:12, nrow = 3, byrow = TRUE)  # 3 "functions" on a 4-point grid
apply_wrapper(F)                            # row means (uniform grid weights)
#> [1]  2.5  6.5 10.5
apply_wrapper(F, wrapper = max)             # row maxes
#> [1]  4  8 12
apply_wrapper(F, wrapper = function(f) f[2]) # point evaluation at grid idx 2
#> [1]  2  6 10
```
