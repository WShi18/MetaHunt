# Build a shared evaluation grid from a reference dataset

Constructs a data frame of grid points suitable for
[`f_hat_from_models()`](https://wshi18.github.io/MetaHunt/reference/f_hat_from_models.md)
from any reference patient-level dataset. This is convenient when the
patient-level covariate space is multidimensional and there is no
obvious one-dimensional grid.

## Usage

``` r
build_grid(reference_data, n_grid = NULL, seed = NULL)
```

## Arguments

- reference_data:

  A data frame (or matrix) of reference patient-level covariates. May be
  a held-out target-population sample, the pooled source covariates, or
  any plausible reference distribution.

- n_grid:

  Optional integer giving the desired grid size. If `NULL` or
  `>= nrow(reference_data)`, the full reference data is returned.

- seed:

  Optional integer seed for reproducibility (used only when
  sub-sampling).

## Value

A data frame of grid points.

## Details

If `n_grid` is `NULL` or at least `nrow(reference_data)`, the reference
data is returned unchanged. Otherwise `n_grid` rows are sampled
uniformly at random (without replacement). The reference data should be
on the same scale and have the same columns as the data each centre's
model was fitted on.

The empirical distribution of the returned grid implicitly defines the
\\\mu\\ measure used downstream. Pass uniform `grid_weights` (the
default) to weight each grid point equally; pass non-uniform
`grid_weights` to weight by an external reference distribution.

## Examples

``` r
set.seed(1)
ref <- data.frame(age = rnorm(500, 60, 10),
                  bp  = rnorm(500, 130, 15),
                  sex = sample(c("F", "M"), 500, replace = TRUE))
grid <- build_grid(ref, n_grid = 50, seed = 1)
nrow(grid)
#> [1] 50
head(grid)
#>          age        bp sex
#> 324 41.30211 131.21499   M
#> 167 57.44973 109.91799   M
#> 129 53.18340 125.20321   M
#> 418 57.48835 140.33866   F
#> 471 51.86756 136.83815   M
#> 299 59.49434  91.05833   F
```
