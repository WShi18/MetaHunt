# Reconstruction-error curve for basis-rank selection

For each candidate number of bases `K`, run
[`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md)
followed by
[`project_to_simplex()`](https://wshi18.github.io/MetaHunt/reference/project_to_simplex.md)
and report the average projection residual \$\$\mathcal E(K) =
\frac{1}{m}\sum\_{i=1}^m \left\\\hat f^{(i)} - \sum\_{k=1}^K
\hat\pi\_{ik}\hat g_k\right\\\_{L^2(\mu)}.\$\$ Plotting `error` against
`K` typically shows an elbow.

## Usage

``` r
reconstruction_error_curve(
  F_hat,
  K_range = NULL,
  grid_weights = NULL,
  dfspa_args = list()
)
```

## Arguments

- F_hat:

  An `m`-by-`G_grid` numeric matrix of study-level function evaluations
  on the shared grid.

- K_range:

  Integer vector of candidate `K` values. Defaults to
  `2:min(nrow(F_hat) - 1, 10)`.

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector used for the
  \\L^2(\mu)\\ norm; defaults to uniform.

- dfspa_args:

  Named list of extra arguments passed to
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md),
  e.g. `list(denoise = FALSE, N = 2)`.

## Value

A data frame with columns `K` (integer) and `error` (numeric). Rows
where [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md)
or the projection fails are reported with `error = NA` and a single
warning summarising the failures.

## Details

This is the unsupervised rank-selection criterion of Section 3.2 of the
paper (Equation for \\\mathcal E(K)\\). It does not require study-level
covariates.

## Examples

``` r
set.seed(1)
G <- 40
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
m <- 50
pi_mat <- matrix(stats::rgamma(m * 3, shape = 0.5), m, 3)
pi_mat <- pi_mat / rowSums(pi_mat)
F_hat  <- pi_mat %*% basis + matrix(stats::rnorm(m * G, sd = 0.02), m, G)

elbow <- reconstruction_error_curve(F_hat, K_range = 2:6)
elbow
#>   K      error
#> 1 2 0.10309809
#> 2 3 0.05037484
#> 3 4 0.04141392
#> 4 5 0.03386482
#> 5 6 0.02471533
```
