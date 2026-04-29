# Project study-level functions onto the simplex spanned by basis functions

For each study `i`, solves the constrained projection \$\$\hat\pi_i =
\arg\min\_{\pi \in \Delta\_{K-1}} \left\\ \hat f^{(i)} - \sum\_{k=1}^K
\pi_k \hat g_k \right\\\_{L^2(\mu)}\$\$ where the norm is the weighted
\\L^2\\ norm defined by `grid_weights`. This is Equation (3) of the
paper and yields the study-specific weights `hat pi_i` used downstream
for weight-model fitting and prediction.

## Usage

``` r
project_to_simplex(F_hat, bases, grid_weights = NULL, ridge = 1e-10)
```

## Arguments

- F_hat:

  An `m`-by-`G_grid` numeric matrix; row `i` is the study function
  \\\hat f^{(i)}\\ evaluated on the shared grid.

- bases:

  A `K`-by-`G_grid` numeric matrix of basis functions on the same grid,
  typically the `bases` slot of a
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md)
  result.

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector of grid weights
  defining the \\L^2(\mu)\\ inner product. Defaults to uniform weights
  `1 / G_grid`.

- ridge:

  Small non-negative scalar added to the diagonal of the QP Hessian for
  numerical stability. Defaults to `1e-10`.

## Value

An `m`-by-`K` numeric matrix of simplex weights; rows sum to 1 and
entries are non-negative.

## Details

The projection reduces to the quadratic program \$\$\min\_{\pi \in
\mathbb R^K}\\ \pi^\top D\\\pi - 2\\d^\top \pi \quad \text{s.t.} \quad
\mathbf 1^\top \pi = 1,\\ \pi \ge 0\$\$ with \\D = G W G^\top\\ and \\d
= G W f^{(i)}\\, where \\G\\ is the `K`-by-`G_grid` basis matrix, \\W =
\mathrm{diag}(\\`grid_weights`\\)\\, and \\f^{(i)}\\ is the `i`-th row
of `F_hat`. Solved via
[`quadprog::solve.QP()`](https://rdrr.io/pkg/quadprog/man/solve.QP.html).
A tiny ridge is added to `D` for numerical stability.

## Examples

``` r
set.seed(1)
G <- 40
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
true_pi <- rbind(diag(3), c(0.4, 0.3, 0.3), c(0.1, 0.7, 0.2))
F_hat <- true_pi %*% basis
fit <- dfspa(F_hat, K = 3, denoise = FALSE)
pi_hat <- project_to_simplex(F_hat, fit$bases)
round(pi_hat, 3)
#>      [,1] [,2] [,3]
#> [1,]  0.0  1.0  0.0
#> [2,]  1.0  0.0  0.0
#> [3,]  0.0  0.0  1.0
#> [4,]  0.3  0.4  0.3
#> [5,]  0.7  0.1  0.2
```
