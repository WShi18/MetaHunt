# MetaHunt

**MetaHunt** is an R package for privacy-preserving meta-analysis of
function-valued quantities (e.g. regression functions, conditional
average treatment effect functions) across heterogeneous studies (Shi,
Imai, & Zhang 2024,
[arXiv:2604.23847](https://arxiv.org/abs/2604.23847)).

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

The package also exports
[`minmax_regret()`](https://wshi18.github.io/MetaHunt/reference/minmax_regret.md),
a covariate-free baseline based on the worst-case-regret aggregator of
Zhang, Huang, and Imai (2024,
[arXiv:2412.11136](https://arxiv.org/abs/2412.11136)). It is useful when
the number of source studies is small or the study-level covariates are
uninformative.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("WShi18/MetaHunt")
```

## Quick example

``` r
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

See
[`vignette("metahunt-intro")`](https://wshi18.github.io/MetaHunt/articles/metahunt-intro.md)
for a full walkthrough including data preparation from fitted model
objects, rank selection, conformal intervals (split / cross / pre-fit),
and the minimax-regret baseline.

## Methods provided

1.  **[`metahunt()`](https://wshi18.github.io/MetaHunt/reference/metahunt.md)** +
    [`predict()`](https://rdrr.io/r/stats/predict.html),
    [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
    [`summary()`](https://rdrr.io/r/base/summary.html) — end-to-end
    fitter.
2.  **[`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md),
    [`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md),
    [`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md)**
    — prediction intervals.
3.  **[`coverage()`](https://wshi18.github.io/MetaHunt/reference/coverage.md),
    [`summary.metahunt_conformal()`](https://wshi18.github.io/MetaHunt/reference/summary.metahunt_conformal.md)**
    — conformal diagnostics.
4.  **[`minmax_regret()`](https://wshi18.github.io/MetaHunt/reference/minmax_regret.md)**
    — covariate-free worst-case-regret baseline (Zhang, Huang, & Imai
    2024).
5.  **[`reconstruction_error_curve()`](https://wshi18.github.io/MetaHunt/reference/reconstruction_error_curve.md),
    [`cv_error_curve()`](https://wshi18.github.io/MetaHunt/reference/cv_error_curve.md)**
    — rank-selection diagnostics.
6.  **[`select_denoising_params()`](https://wshi18.github.io/MetaHunt/reference/select_denoising_params.md)**
    — cross-validated tuning of `(N, Delta)`.
7.  **[`f_hat_from_models()`](https://wshi18.github.io/MetaHunt/reference/f_hat_from_models.md),
    [`build_grid()`](https://wshi18.github.io/MetaHunt/reference/build_grid.md)**
    — onramp from fitted-model lists to the package’s matrix inputs.

The individual pipeline building blocks
([`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md),
[`project_to_simplex()`](https://wshi18.github.io/MetaHunt/reference/project_to_simplex.md),
[`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md),
[`predict_target()`](https://wshi18.github.io/MetaHunt/reference/predict_target.md),
[`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md))
are also exported and can be composed independently.

## Citation

If you use MetaHunt in academic work, please cite:

> Shi, W., Imai, K., and Zhang, Y. (2024). *Privacy-preserving
> meta-analysis through low-rank basis hunting.* arXiv:2604.23847.
> <https://arxiv.org/abs/2604.23847>

## Status

`v0.1.0` is the first GitHub release. The API is stabilising; user
feedback is welcome via
[issues](https://github.com/WShi18/MetaHunt/issues).

## License

MIT © Wenqi Shi, Kosuke Imai, Yi Zhang
