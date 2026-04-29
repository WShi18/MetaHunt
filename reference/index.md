# Package index

## End-to-end fitter

One-call interface; covers the full pipeline.

- [`metahunt()`](https://wshi18.github.io/MetaHunt/reference/metahunt.md)
  : Fit the full MetaHunt pipeline

- [`predict(`*`<metahunt>`*`)`](https://wshi18.github.io/MetaHunt/reference/predict.metahunt.md)
  : Predict target functions (or scalar summaries) from a MetaHunt fit

- [`summary(`*`<metahunt>`*`)`](https://wshi18.github.io/MetaHunt/reference/summary.metahunt.md)
  : Summarise a MetaHunt fit

- [`print(`*`<summary.metahunt>`*`)`](https://wshi18.github.io/MetaHunt/reference/print.summary.metahunt.md)
  :

  Print a `summary.metahunt` object

- [`plot(`*`<metahunt>`*`)`](https://wshi18.github.io/MetaHunt/reference/plot.metahunt.md)
  : Plot recovered basis functions from a MetaHunt fit

## Conformal prediction

Distribution-free prediction intervals around the target function.

- [`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
  : Split conformal prediction intervals for target-function predictions
- [`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md)
  : Cross-conformal prediction intervals (pooled K-fold scores)
- [`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md)
  : Split conformal intervals from a pre-fit MetaHunt pipeline
- [`coverage()`](https://wshi18.github.io/MetaHunt/reference/coverage.md)
  : Empirical coverage of a conformal prediction-interval object
- [`summary(`*`<metahunt_conformal>`*`)`](https://wshi18.github.io/MetaHunt/reference/summary.metahunt_conformal.md)
  : Summarise a conformal prediction-interval object
- [`plot(`*`<metahunt_conformal>`*`)`](https://wshi18.github.io/MetaHunt/reference/plot.metahunt_conformal.md)
  : Plot a conformal prediction-interval object

## Rank and tuning selection

Diagnostics for choosing the number of bases K and the d-fSPA denoising
knobs.

- [`reconstruction_error_curve()`](https://wshi18.github.io/MetaHunt/reference/reconstruction_error_curve.md)
  : Reconstruction-error curve for basis-rank selection
- [`cv_error_curve()`](https://wshi18.github.io/MetaHunt/reference/cv_error_curve.md)
  : Cross-validated prediction-error curve for basis-rank selection
- [`select_denoising_params()`](https://wshi18.github.io/MetaHunt/reference/select_denoising_params.md)
  : Choose d-fSPA denoising parameters by cross-validation
- [`print(`*`<metahunt_denoising_search>`*`)`](https://wshi18.github.io/MetaHunt/reference/print.metahunt_denoising_search.md)
  : Print method for d-fSPA denoising parameter search results

## Pipeline building blocks

Lower-level primitives that can be composed independently.

- [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md) :
  Denoised functional Successive Projection Algorithm (d-fSPA)
- [`project_to_simplex()`](https://wshi18.github.io/MetaHunt/reference/project_to_simplex.md)
  : Project study-level functions onto the simplex spanned by basis
  functions
- [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md)
  : Fit a weight model mapping study-level covariates to simplex weights
- [`predict(`*`<metahunt_weight_model>`*`)`](https://wshi18.github.io/MetaHunt/reference/predict.metahunt_weight_model.md)
  : Predict simplex weights for new study-level covariates
- [`coef(`*`<metahunt_weight_model>`*`)`](https://wshi18.github.io/MetaHunt/reference/coef.metahunt_weight_model.md)
  : Extract coefficients from a MetaHunt weight model
- [`predict_target()`](https://wshi18.github.io/MetaHunt/reference/predict_target.md)
  : Predict the target function for new study-level covariates
- [`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md)
  : Reduce predicted functions to scalars via a user-supplied wrapper

## Data preparation

Onramp from fitted-model lists to MetaHunt’s matrix inputs.

- [`build_grid()`](https://wshi18.github.io/MetaHunt/reference/build_grid.md)
  : Build a shared evaluation grid from a reference dataset

- [`f_hat_from_models()`](https://wshi18.github.io/MetaHunt/reference/f_hat_from_models.md)
  :

  Build the `F_hat` matrix from a list of fitted study-level models

## Baselines

Covariate-free worst-case-regret aggregator (Zhang, Huang & Imai 2024).

- [`minmax_regret()`](https://wshi18.github.io/MetaHunt/reference/minmax_regret.md)
  : Minimax-regret aggregator for multisite function-valued estimands

## Package

- [`MetaHunt`](https://wshi18.github.io/MetaHunt/reference/MetaHunt-package.md)
  [`MetaHunt-package`](https://wshi18.github.io/MetaHunt/reference/MetaHunt-package.md)
  : MetaHunt: Privacy-Preserving Meta-Analysis via Low-Rank Basis
  Hunting
