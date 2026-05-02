# Denoised functional Successive Projection Algorithm (d-fSPA)

Recovers a set of `K` latent basis functions from a collection of
study-level function estimates under the low-rank cross-study
heterogeneity assumption of Shi, Imai, and Zhang. Implements Algorithm 1
of the paper ("The d-fSPA Algorithm for basis hunting").

## Usage

``` r
dfspa(F_hat, K, grid_weights = NULL, N = NULL, Delta = NULL, denoise = TRUE)
```

## Arguments

- F_hat:

  An `m`-by-`G` numeric matrix where row `i` is the evaluation of the
  estimated function \\\hat f^{(i)}\\ at `G` grid points.

- K:

  Integer number of basis functions to recover. Must satisfy
  `1 <= K <= m` after denoising.

- grid_weights:

  Optional length-`G` non-negative numeric vector of grid weights
  defining the \\L^2(\mu)\\ inner product. Defaults to uniform weights
  `1 / G`.

- N, Delta:

  Optional numeric tuning parameters controlling denoising. See Details.

- denoise:

  Logical; if `FALSE`, the denoising step is skipped and plain fSPA is
  run. Defaults to `TRUE`.

## Value

An object of class `"dfspa"`: a list containing

- `bases`:

  A `K`-by-`G` matrix whose rows are the recovered basis functions
  evaluated on the grid (denoised, if applicable).

- `selected`:

  Length-`K` integer vector of the selected row indices into the
  post-denoising function matrix.

- `original_indices`:

  Length-`K` integer vector of the selected study indices in the
  original input `F_hat` (before any rows were dropped by denoising).

- `kept`:

  Integer vector of row indices of `F_hat` that survived denoising.

- `F_denoised`:

  The post-denoising function matrix (`length(kept)`-by-`G`).

- `grid_weights`:

  Grid weights used.

- `N`, `Delta`:

  Tuning parameters actually used (or `NA` when `denoise = FALSE`).

- `K`:

  Number of bases requested.

- `call`:

  The matched call.

## Details

Each study-level function is represented by its evaluations on a shared
grid of `G` points. The (weighted) \\L^2(\mu)\\ inner product is
\\\langle f,g\rangle = \sum\_{j=1}^G w_j f(x_j) g(x_j)\\, where the
`grid_weights` `w_j` are proportional to the measure \\\mu\\. If not
supplied, uniform weights `1 / G` are used.

Denoising follows Jin (2024): for each study `i`, let \\B\_\Delta(\hat
f^{(i)}) = \\j : \\\hat f^{(j)} - \hat f^{(i)}\\ \le \Delta\\\\. If
\\\|B\_\Delta(\hat f^{(i)})\| \< N\\, study `i` is discarded; otherwise
\\\hat f^{(i)}\\ is replaced by the average of the functions in
\\B\_\Delta\\. After denoising, the functional SPA step iteratively
selects, at each of the `K` iterations, the remaining function with the
largest norm after projecting out the span of previously selected bases.

Default tuning parameters follow the heuristics of the paper:
`N = 0.5 * log(m)` and \\\Delta = \max\_{ij} \\\hat f^{(i)} - \hat
f^{(j)}\\ / 10\\.

## Examples

``` r
set.seed(1)
G <- 50
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)          # 3 true bases
pi_mat <- rbind(diag(3),                              # 3 pure studies
                c(0.5, 0.3, 0.2),
                c(0.2, 0.5, 0.3),
                c(0.3, 0.3, 0.4))
F_hat <- pi_mat %*% basis                             # m = 6, G = 50
fit <- dfspa(F_hat, K = 3, denoise = FALSE)
fit$original_indices    # should be a permutation of 1, 2, 3
#> [1] 2 1 3
```
