# Plot a conformal prediction-interval object

For pointwise objects (no `wrapper` was used), draws the predicted
function for one target study together with its pointwise band. For
scalar objects, draws point predictions with whisker error bars over all
targets.

## Usage

``` r
# S3 method for class 'metahunt_conformal'
plot(
  x,
  target_idx = 1L,
  x_axis = NULL,
  fill = grDevices::adjustcolor("steelblue", 0.2),
  line_col = "steelblue",
  ...
)
```

## Arguments

- x:

  A `"metahunt_conformal"` object.

- target_idx:

  For pointwise objects, the integer index of the target study to plot
  (default `1`). Ignored for scalar objects.

- x_axis:

  Optional numeric vector of length `G_grid` giving the x-axis values
  for pointwise plotting. Defaults to `seq_len(G_grid)`.

- fill:

  Polygon fill colour for the band. Default semi-transparent steel blue.

- line_col:

  Line colour for the predicted function (or points in scalar mode).
  Default steel blue.

- ...:

  Additional graphical parameters passed to the underlying plotting
  calls.

## Value

Invisibly returns `x`.
