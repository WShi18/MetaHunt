# Minimax-regret aggregator for multisite function-valued estimands

Implements the minimax-regret estimator of Zhang, Huang, and Imai
([arXiv:2412.11136](https://arxiv.org/abs/2412.11136)) for aggregating
site-level function estimates. Given study-level functions \\\hat
f^{(i)}\\, the estimator is \$\$\hat q = \arg\min\_{q \in
\Delta\_{m-1}}\\ q^\top \hat\Gamma\\ q - \hat d^\top q, \qquad
\hat\Gamma\_{ij} = \sum_g w_g\\ \hat f^{(i)}(x_g)\hat f^{(j)}(x_g),\\ \\
\hat d_i = \sum_g w_g\\ (\hat f^{(i)}(x_g))^2,\$\$ yielding the
predicted target function \\\tilde f(x) = \sum\_{i=1}^m \hat q_i\\ \hat
f^{(i)}(x)\\. Unlike
[`metahunt()`](https://wshi18.github.io/MetaHunt/reference/metahunt.md),
this method does not use study-level covariates; the target is the
worst-case-regret aggregator over the convex hull of source functions.

## Usage

``` r
minmax_regret(F_hat, grid_weights = NULL, ridge = 1e-10, wrapper = NULL)
```

## Arguments

- F_hat:

  An `m`-by-`G_grid` numeric matrix of source-site function evaluations
  on a shared grid; row `i` is \\\hat f^{(i)}\\.

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector defining the
  target measure used to compute \\\hat\Gamma\\ and \\\hat d\\. Defaults
  to uniform weights `1 / G_grid`.

- ridge:

  Non-negative scalar; replaces \\\hat\Gamma\\ with \\\hat\Gamma +
  \text{ridge}\cdot I\\ for numerical stability. Defaults to `1e-10`.

- wrapper:

  Optional reduction function (see
  [`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md)).
  If `NULL`, `prediction` is the length-`G_grid` predicted target
  function; if a function, `prediction` is the scalar
  `wrapper(prediction)` (e.g. an ATE when `wrapper = mean` and rows of
  `F_hat` are CATE evaluations).

## Value

An object of class `"minmax_regret"`: a list with

- `prediction`:

  Predicted target. Length-`G_grid` vector when `wrapper = NULL`; scalar
  otherwise.

- `q`:

  Length-`m` simplex weights from the minimax-regret QP.

- `Gamma`:

  The `m`-by-`m` Gram matrix used (post-ridge).

- `d`:

  Length-`m` linear coefficient vector.

- `grid_weights`:

  Grid weights used.

- `ridge`:

  Ridge value used.

- `wrapper`:

  Wrapper function or `NULL`.

## Details

The simplex-constrained QP is solved with
[`quadprog::solve.QP()`](https://rdrr.io/pkg/quadprog/man/solve.QP.html).
A small ridge is added to \\\hat\Gamma\\ for numerical stability when
source functions are highly collinear. The resulting `q` is clipped to
be non-negative and renormalised to sum to 1 to absorb floating-point
drift.

## References

Zhang, Y., Huang, M., and Imai, K. (2024). Minimax regret estimation for
generalizing heterogeneous treatment effects with multisite data.
[arXiv:2412.11136](https://arxiv.org/abs/2412.11136).

## Examples

``` r
set.seed(1)
G <- 30; m <- 6
x <- seq(0, 1, length.out = G)
F_hat <- rbind(
  sin(pi * x),
  cos(pi * x),
  x,
  0.5 * sin(pi * x) + 0.5 * x,
  0.3 * cos(pi * x) + 0.7 * x,
  0.4 * sin(pi * x) + 0.4 * cos(pi * x) + 0.2 * x
)

fit <- minmax_regret(F_hat)
fit$q                                        # simplex weights over sources
#> [1] 0.0 0.5 0.5 0.0 0.0 0.0
length(fit$prediction)                       # G: predicted target function
#> [1] 30

# ATE-style scalar via wrapper
minmax_regret(F_hat, wrapper = mean)$prediction
#> [1] 0.25
```
