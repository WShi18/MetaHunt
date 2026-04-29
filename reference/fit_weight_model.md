# Fit a weight model mapping study-level covariates to simplex weights

Given a matrix of simplex-valued weights \\\hat\pi_1,\ldots,\hat\pi_m\\
(e.g. from
[`project_to_simplex()`](https://wshi18.github.io/MetaHunt/reference/project_to_simplex.md))
and associated study-level covariates \\\mathbf W_1,\ldots,\mathbf
W_m\\, fit a model \\\widehat{\mathcal M}:\mathbf W \mapsto
\boldsymbol\pi\\. The default method is Dirichlet regression via the
`DirichletReg` package.

## Usage

``` r
fit_weight_model(
  pi_hat,
  W,
  method = c("dirichlet"),
  boundary_eps = 1e-04,
  formula = NULL,
  ...
)
```

## Arguments

- pi_hat:

  An `m`-by-`K` numeric matrix of simplex weights; rows must be
  non-negative and sum to 1 (up to tolerance `1e-6`).

- W:

  An `m`-by-`p` matrix or data frame of study-level covariates.

- method:

  Weight-model method. Currently only `"dirichlet"` is supported.

- boundary_eps:

  Small positive scalar used to shrink weights away from the simplex
  boundary before Dirichlet fitting. Defaults to `1e-4`.

- formula:

  Optional RHS-only formula (e.g. `~ x1 + I(x2^2)`) describing the
  covariate part of the Dirichlet regression. Defaults to `~ .` (all
  columns of `W`).

- ...:

  Passed through to
  [`DirichletReg::DirichReg()`](https://rdrr.io/pkg/DirichletReg/man/dirichreg.html).

## Value

An object of class `"metahunt_weight_model"`: a list with the fitted
model, formula, method, `K`, and training covariate names.

## Details

Dirichlet regression cannot handle weights exactly at the simplex
boundary (`0` or `1`), which frequently arise after constrained
projection. Before fitting, rows of `pi_hat` are shrunk toward the
barycenter via \\\tilde\pi = (\pi + \varepsilon) / (1 + K\varepsilon)\\,
with \\\varepsilon\\ set by `boundary_eps`.

## Examples

``` r
# \donttest{
set.seed(1)
m <- 80; K <- 3; p <- 2
W <- matrix(rnorm(m * p), m, p); colnames(W) <- c("w1", "w2")
# generate simplex weights driven by W
eta <- cbind(0.5 * W[, 1], -0.3 * W[, 2], rep(0, m))
pi_true <- exp(eta) / rowSums(exp(eta))
pi_hat <- pi_true + matrix(rnorm(m * K, sd = 0.01), m, K)
pi_hat <- pmax(pi_hat, 0); pi_hat <- pi_hat / rowSums(pi_hat)
model <- fit_weight_model(pi_hat, W)
predict(model, newdata = matrix(c(0, 0), 1, 2, dimnames = list(NULL, c("w1","w2"))))
#>           [,1]      [,2]      [,3]
#> [1,] 0.3328722 0.3340615 0.3330662
# }
```
