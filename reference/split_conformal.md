# Split conformal prediction intervals for target-function predictions

Implements Algorithm 2 of the paper (split conformal prediction) over
the MetaHunt pipeline. Studies are partitioned into training and
calibration sets. The training set is used to fit d-fSPA, the
constrained projection, and the weight model; the calibration set
supplies conformity scores, which determine the width of the intervals.

## Usage

``` r
split_conformal(
  F_hat,
  W,
  W_new,
  K,
  alpha = 0.05,
  cal_frac = 0.3,
  wrapper = NULL,
  grid_weights = NULL,
  calibration_idx = NULL,
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

- W_new:

  A matrix or data frame of new target covariates. Must contain columns
  matching `W`.

- K:

  Integer number of basis functions.

- alpha:

  Miscoverage level; interval has nominal coverage \\1-\alpha\\. Default
  `0.05`.

- cal_frac:

  Numeric in `(0, 1)` giving the fraction of studies in the calibration
  set. Default `0.3`. Ignored if `calibration_idx` is supplied.

- wrapper:

  Optional reduction function (see
  [`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md)).
  If `NULL`, intervals are constructed pointwise at every grid point.

- grid_weights:

  Optional length-`G_grid` non-negative numeric vector used for the
  \\L^2(\mu)\\ norm and for the default wrapper.

- calibration_idx:

  Optional integer vector of row indices in `F_hat` to use as the
  calibration set. If supplied, `cal_frac` is ignored.

- dfspa_args, weight_model_args:

  Named lists passed to
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md) and
  [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md)
  respectively.

- seed:

  Optional integer seed for reproducible train/calibration splits.

## Value

An object of class `"metahunt_conformal"`: a list with

- `prediction`:

  Point predictions for `W_new`. A numeric vector of length
  `nrow(W_new)` in the scalar case, or an `nrow(W_new)`-by-`G_grid`
  matrix in the pointwise case.

- `lower`, `upper`:

  Interval endpoints, same shape as `prediction`.

- `alpha`:

  Miscoverage level used.

- `method`:

  `"split"`.

- `n_cal`:

  Calibration sample size.

- `quantile`:

  The conformal quantile: a scalar (scalar case) or a length-`G_grid`
  vector (pointwise case).

- `wrapper`:

  The wrapper used, or `NULL`.

## Details

Given a target function, one can either construct intervals **pointwise
at every grid point** (when `wrapper = NULL`) or for a **scalar
summary** of the target function (when `wrapper` is a function).

- **Pointwise** (`wrapper = NULL`): for each grid point `g` the
  conformity score is \\R\_{i,g} = \|\hat f^{(i)}(x_g) - \tilde
  f^{(i)}(x_g)\|\\ across calibration studies `i`. A separate
  \\(1-\alpha)\\-quantile \\q_g\\ is computed per grid point, and the
  interval at grid point `g` for target `j` is \\\[\tilde
  f^{(j)}(x_g)-q_g,\\ \tilde f^{(j)}(x_g)+q_g\]\\.

- **Scalar** (`wrapper` supplied): conformity scores are \\R_i =
  \|s(\hat f^{(i)}) - s(\tilde f^{(i)})\|\\ with `s = wrapper`, and the
  interval for each target is \\\[s(\tilde f^{(j)})-q,\\ s(\tilde
  f^{(j)})+q\]\\ with a single shared quantile \\q\\.

The finite-sample quantile is \\q = R\_{(k)}\\ with \\k =
\lceil(1-\alpha)(n\_\mathrm{cal}+1)\rceil\\; if \\k \>
n\_\mathrm{cal}\\, `q = Inf` and intervals are \\(-\infty, \infty)\\.

## Examples

``` r
set.seed(1)
G <- 40; m <- 80; K_true <- 3
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

W_new <- data.frame(w1 = c(0, 1), w2 = c(0, -1))
# pointwise intervals at every grid point
pi_grid <- split_conformal(F_hat, W, W_new, K = 3, seed = 1)
dim(pi_grid$lower)  # 2 x 40
#> [1]  2 40
# scalar intervals for the grid-weighted mean (ATE-style)
pi_ate  <- split_conformal(F_hat, W, W_new, K = 3, wrapper = mean, seed = 1)
pi_ate$prediction
#> [1] 0.3668234 0.5851352
```
