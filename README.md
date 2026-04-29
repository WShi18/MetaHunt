# MetaHunt

<!-- badges: start -->
[![R-CMD-check](https://github.com/WShi18/MetaHunt/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/WShi18/MetaHunt/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/WShi18/MetaHunt/actions/workflows/pkgdown.yaml/badge.svg)](https://wshi18.github.io/MetaHunt/)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**MetaHunt** is an R package for privacy-preserving meta-analysis of
function-valued quantities (e.g. regression functions, conditional
average treatment effect functions) across heterogeneous studies
(Shi, Imai, & Zhang 2024, [arXiv:2604.23847](https://arxiv.org/abs/2604.23847)).

Documentation: <https://wshi18.github.io/MetaHunt/>

The package implements the methodology of

> Shi, W., Imai, K., and Zhang, Y. *Privacy-preserving Meta-analysis
> through Low-Rank Basis Hunting.*

Given aggregate-level estimates of a function of interest from several
source studies, together with study-level covariates, MetaHunt recovers
a small set of latent basis functions that span all studies under a
low-rank assumption, models how the mixing weights depend on study-level
covariates via Dirichlet regression, and predicts the corresponding
function for a new target population with valid conformal prediction
intervals. Individual-level data from source studies are not required,
which makes the procedure suitable for federated and multi-site
settings.

The package also exports `minmax_regret()`, a covariate-free baseline
based on the worst-case-regret aggregator of Zhang, Huang, and Imai
(2024, [arXiv:2412.11136](https://arxiv.org/abs/2412.11136)). It is
useful when the number of source studies is small or the study-level
covariates are uninformative.

## Installation

```r
# install.packages("remotes")
remotes::install_github("WShi18/MetaHunt")
```

## Quick example

```r
library(MetaHunt)
set.seed(1)

# Simulate m studies, each represented by its function on a shared grid.
G <- 30; m <- 60
x <- seq(0, 1, length.out = G)
basis <- rbind(sin(pi * x), cos(pi * x), x)
W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
eta <- as.matrix(W) %*% cbind(c(1, -0.5), c(-0.4, 1), c(0, 0))
pi_true <- exp(eta) / rowSums(exp(eta))
F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)

# Fit and predict for a new study with W_new = (0, 0).
fit   <- metahunt(F_hat, W, K = 3, dfspa_args = list(denoise = FALSE))
W_new <- data.frame(w1 = 0, w2 = 0)
pred  <- predict(fit, newdata = W_new, wrapper = mean)

# 95% conformal interval for the same target.
ci <- split_conformal(F_hat, W, W_new, K = 3,
                      wrapper = mean, alpha = 0.05, seed = 1,
                      dfspa_args = list(denoise = FALSE))
ci$prediction
c(ci$lower, ci$upper)
```

See `vignette("metahunt-intro")` for a full walkthrough including
data preparation from fitted model objects, rank selection, conformal
intervals (split / cross / pre-fit), and the minimax-regret baseline.

## Methods provided

1. **`metahunt()`** + `predict()`, `plot()`, `summary()` — end-to-end fitter.
2. **`split_conformal()`, `cross_conformal()`, `conformal_from_fit()`** — prediction intervals.
3. **`coverage()`, `summary.metahunt_conformal()`** — conformal diagnostics.
4. **`minmax_regret()`** — covariate-free worst-case-regret baseline (Zhang, Huang, & Imai 2024).
5. **`reconstruction_error_curve()`, `cv_error_curve()`** — rank-selection diagnostics.
6. **`select_denoising_params()`** — cross-validated tuning of `(N, Delta)`.
7. **`f_hat_from_models()`, `build_grid()`** — onramp from fitted-model lists to the package's matrix inputs.

The individual pipeline building blocks (`dfspa()`,
`project_to_simplex()`, `fit_weight_model()`, `predict_target()`,
`apply_wrapper()`) are also exported and can be composed independently.

## Citation

If you use MetaHunt in academic work, please cite:

> Shi, W., Imai, K., and Zhang, Y. (2024). *Privacy-preserving
> meta-analysis through low-rank basis hunting.* arXiv:2604.23847.
> <https://arxiv.org/abs/2604.23847>

## Status

`v0.1.0` is the first GitHub release. The API is stabilising; user
feedback is welcome via [issues](https://github.com/WShi18/MetaHunt/issues).

## License

MIT © Wenqi Shi, Kosuke Imai, Yi Zhang
