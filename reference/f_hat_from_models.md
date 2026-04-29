# Build the `F_hat` matrix from a list of fitted study-level models

Many users arrive at MetaHunt with one fitted model per study (e.g. a
`ranger::ranger()` random forest or a
[`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html))
and a chosen evaluation grid. `f_hat_from_models()` evaluates each model
on the shared grid and stacks the predictions into the `m`-by-`G_grid`
matrix the rest of the package expects.

## Usage

``` r
f_hat_from_models(models, grid, predict_fn = NULL)
```

## Arguments

- models:

  A non-empty list of fitted model objects, one per study.

- grid:

  A data frame (or matrix) of grid points at which to evaluate each
  model. Same columns as the data each model was fitted on.

- predict_fn:

  Optional `function(model, grid)` returning a numeric vector; defaults
  to the class-aware dispatcher above.

## Value

An `length(models)`-by-`nrow(grid)` numeric matrix; row `i` is model `i`
evaluated at every row of `grid`.

## Details

By default the function dispatches on the model class:

- `ranger` objects are evaluated as
  `predict(model, data = grid)$predictions`.

- Objects inheriting from `causal_forest` or `grf` are evaluated as
  `predict(model, newdata = grid)$predictions`.

- All other classes fall through to
  `as.numeric(predict(model, newdata = grid))`, which works for `lm`,
  `glm`, `randomForest`, and most other R model objects.

Override the dispatch with `predict_fn = function(model, grid) ...` if
your models need bespoke handling. The function must return a
length-`G_grid` numeric vector for each model.

All rows of the returned matrix must have the same length (`G_grid`) and
contain no `NA` values; the function errors otherwise.

## See also

[`build_grid()`](https://wshi18.github.io/MetaHunt/reference/build_grid.md)
to construct `grid` from a reference dataset.

## Examples

``` r
# Toy example: each "centre" fits a polynomial regression
set.seed(1)
make_centre_data <- function(slope) {
  x <- runif(60)
  data.frame(x = x, y = slope * x + rnorm(60, sd = 0.1))
}
models <- lapply(c(-1, 0, 1, 0.5, -0.5), function(s)
  stats::lm(y ~ poly(x, 2), data = make_centre_data(s)))

grid  <- data.frame(x = seq(0, 1, length.out = 30))
F_hat <- f_hat_from_models(models, grid)
dim(F_hat)   # 5 x 30
#> [1]  5 30
```
