# Choose d-fSPA denoising parameters by cross-validation

At a fixed `K`, performs k-fold cross-validation over a grid of
denoising parameter pairs `(N, Delta)` for
[`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md). For
each candidate pair and each fold, the full MetaHunt pipeline is fit on
the out-of-fold studies and predicts the held-out studies' functions.
The pair with the lowest average prediction error is selected.

## Usage

``` r
select_denoising_params(
  F_hat,
  W,
  K,
  N_grid = NULL,
  Delta_grid = NULL,
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

- K:

  Integer number of basis functions (fixed during this search).

- N_grid:

  Optional numeric vector of candidate `N` values. Defaults to
  `c(0.2, 0.5, 1.0) * log(nrow(F_hat))`.

- Delta_grid:

  Optional numeric vector of candidate `Delta` values. Defaults to
  `c(0.05, 0.10, 0.20, 0.30)` times the maximum pairwise \\L^2(\mu)\\
  distance among studies.

- n_folds:

  Integer number of folds (default `5`).

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector.

- dfspa_args:

  Named list of additional arguments for
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md)
  (e.g. [`list()`](https://rdrr.io/r/base/list.html); do not include
  `N`, `Delta`, or `denoise` here, they are set by the search).

- weight_model_args:

  Named list of additional arguments for
  [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md).

- seed:

  Optional integer seed for reproducible fold assignment.

## Value

An object of class `metahunt_denoising_search`: a list with

- `grid`:

  A data frame with one row per `(N, Delta)` pair, columns `N`, `Delta`,
  `cv_error`, `cv_se`, `n_folds_ok`.

- `best`:

  A list with the `(N, Delta)` minimising `cv_error`.

- `K`, `n_folds`, `grid_weights`:

  Inputs echoed back for traceability.

## Details

This is the cross-validated tuning of the denoising parameters discussed
in Section 3.1 of the paper. Joint tuning over `(K, N, Delta)` is not
supported because it scales poorly; if you also want to choose `K`, do
it first via
[`cv_error_curve()`](https://wshi18.github.io/MetaHunt/reference/cv_error_curve.md)
and then call this function at the selected `K`.

Default candidate grids:

- `N_grid = m * c(NA, NA, NA)` resolved at runtime to
  `c(0.2, 0.5, 1.0) * log(m)`.

- `Delta_grid = max_pairwise_dist * c(0.05, 0.10, 0.20, 0.30)`.

Pass your own `N_grid` / `Delta_grid` (in original units) to override.

## See also

[`cv_error_curve()`](https://wshi18.github.io/MetaHunt/reference/cv_error_curve.md)
for selecting `K`,
[`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md) for
the underlying basis-hunting algorithm.

## Examples

``` r
# \donttest{
set.seed(1)
G <- 30; m <- 80
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
eta <- as.matrix(W) %*% cbind(c(1, -0.5), c(-0.4, 1), c(0, 0))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

tune <- select_denoising_params(F_hat, W, K = 3, n_folds = 4, seed = 1)
#> Warning: 8 (N, Delta, fold) combinations failed:
#>   (N=2.19101, Delta=0.0479349, fold 1): Only 0 studies survive denoising but K = 3; consider increasing `Delta`, decreasing `N`, or reducing `K`.
#>   (N=2.19101, Delta=0.0479349, fold 2): Only 0 studies survive denoising but K = 3; consider increasing `Delta`, decreasing `N`, or reducing `K`.
#>   (N=2.19101, Delta=0.0479349, fold 3): Only 0 studies survive denoising but K = 3; consider increasing `Delta`, decreasing `N`, or reducing `K`.
#>   (N=2.19101, Delta=0.0479349, fold 4): Only 0 studies survive denoising but K = 3; consider increasing `Delta`, decreasing `N`, or reducing `K`.
#>   (N=4.38203, Delta=0.0479349, fold 1): Only 0 studies survive denoising but K = 3; consider increasing `Delta`, decreasing `N`, or reducing `K`.
#>   (N=4.38203, Delta=0.0479349, fold 2): Only 0 studies survive denoising but K = 3; consider increasing `Delta`, decreasing `N`, or reducing `K`.
#>   (N=4.38203, Delta=0.0479349, fold 3): Only 0 studies survive denoising but K = 3; consider increasing `Delta`, decreasing `N`, or reducing `K`.
#>   (N=4.38203, Delta=0.0479349, fold 4): Only 0 studies survive denoising but K = 3; consider increasing `Delta`, decreasing `N`, or reducing `K`.
tune$grid
#>            N      Delta   cv_error       cv_se n_folds_ok
#> 1  0.8764053 0.04793491 0.07080549 0.002541663          4
#> 2  2.1910133 0.04793491        NaN          NA          0
#> 3  4.3820266 0.04793491        NaN          NA          0
#> 4  0.8764053 0.09586983 0.07454212 0.003150261          4
#> 5  2.1910133 0.09586983 0.07890784 0.003902666          4
#> 6  4.3820266 0.09586983 0.08239793 0.003161490          4
#> 7  0.8764053 0.19173966 0.08429097 0.003326955          4
#> 8  2.1910133 0.19173966 0.08429097 0.003326955          4
#> 9  4.3820266 0.19173966 0.08525947 0.003824418          4
#> 10 0.8764053 0.28760949 0.09502216 0.006381783          4
#> 11 2.1910133 0.28760949 0.09502216 0.006381783          4
#> 12 4.3820266 0.28760949 0.09502216 0.006381783          4
tune$best
#> $N
#> [1] 0.8764053
#> 
#> $Delta
#> [1] 0.04793491
#> 
#> $cv_error
#> [1] 0.07080549
#> 
# }
```
