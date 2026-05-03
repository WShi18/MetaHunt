# Changelog

## MetaHunt 0.1.0

Initial public release.

### End-to-end fitter

- [`metahunt()`](https://wshi18.github.io/MetaHunt/reference/metahunt.md)
  chains denoised functional SPA basis hunting, constrained simplex
  projection, and Dirichlet weight modelling in a single call. Method
  dispatch for [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the
  returned `"metahunt"` object.

### Conformal prediction

- [`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md)
  and
  [`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md)
  return distribution-free prediction intervals around the target
  function (pointwise on the grid or, with a `wrapper`, around a scalar
  summary).
- [`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md)
  adds intervals to an already-fit pipeline using a held-out calibration
  set.
- [`coverage()`](https://wshi18.github.io/MetaHunt/reference/coverage.md),
  [`summary()`](https://rdrr.io/r/base/summary.html), and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods for
  the `"metahunt_conformal"` class.

### Rank and tuning selection

- [`reconstruction_error_curve()`](https://wshi18.github.io/MetaHunt/reference/reconstruction_error_curve.md)
  (unsupervised elbow) and
  [`cv_error_curve()`](https://wshi18.github.io/MetaHunt/reference/cv_error_curve.md)
  (supervised CV) for picking `K`.
- [`select_denoising_params()`](https://wshi18.github.io/MetaHunt/reference/select_denoising_params.md)
  cross-validates the `(N, Delta)` knobs of
  [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md).

### Pipeline building blocks

- [`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md)
  denoised functional Successive Projection Algorithm (Algorithm 1 of
  the paper).
- [`project_to_simplex()`](https://wshi18.github.io/MetaHunt/reference/project_to_simplex.md)
  constrained simplex projection of each study’s function onto the
  recovered bases (quadratic program via `quadprog`).
- [`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md)
  and
  [`predict.metahunt_weight_model()`](https://wshi18.github.io/MetaHunt/reference/predict.metahunt_weight_model.md)
  for Dirichlet regression of simplex weights on study-level covariates,
  with
  [`coef.metahunt_weight_model()`](https://wshi18.github.io/MetaHunt/reference/coef.metahunt_weight_model.md)
  for inspecting coefficients.
- [`predict_target()`](https://wshi18.github.io/MetaHunt/reference/predict_target.md)
  and
  [`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md)
  for composing predictions and scalar summaries by hand.

### Data preparation

- [`build_grid()`](https://wshi18.github.io/MetaHunt/reference/build_grid.md)
  constructs a shared evaluation grid from any reference patient-level
  dataset.
- [`f_hat_from_models()`](https://wshi18.github.io/MetaHunt/reference/f_hat_from_models.md)
  evaluates a list of fitted models on the shared grid with class-aware
  dispatch for `ranger`, `grf` (`causal_forest`, `regression_forest`),
  and a default branch that covers `lm`/`glm`/`randomForest`. Custom S4
  classes can supply their own `predict_fn`.

### Baselines

- [`minmax_regret()`](https://wshi18.github.io/MetaHunt/reference/minmax_regret.md)
  implements the covariate-free worst-case-regret aggregator of Zhang,
  Huang, and Imai (2024,
  [arXiv:2412.11136](https://arxiv.org/abs/2412.11136)).

### Documentation

- Tutorials: `metahunt-intro`, `data-prep`, `grid-weights`,
  `wrapper-scalar`, plus `get-started`.
- Companion paper: Shi, Imai, and Zhang (2024,
  [arXiv:2604.23847](https://arxiv.org/abs/2604.23847)).
