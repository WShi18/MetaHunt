# MetaHunt 0.1.0

Initial public release.

## End-to-end fitter

- `metahunt()` chains denoised functional SPA basis hunting,
  constrained simplex projection, and Dirichlet weight modelling in a
  single call. Method dispatch for `predict()`, `summary()`, and
  `plot()` on the returned `"metahunt"` object.

## Conformal prediction

- `split_conformal()` and `cross_conformal()` return distribution-free
  prediction intervals around the target function (pointwise on the
  grid or, with a `wrapper`, around a scalar summary).
- `conformal_from_fit()` adds intervals to an already-fit pipeline
  using a held-out calibration set.
- `coverage()`, `summary()`, and `plot()` methods for the
  `"metahunt_conformal"` class.

## Rank and tuning selection

- `reconstruction_error_curve()` (unsupervised elbow) and
  `cv_error_curve()` (supervised CV) for picking `K`.
- `select_denoising_params()` cross-validates the `(N, Delta)` knobs
  of `dfspa()`.

## Pipeline building blocks

- `dfspa()` denoised functional Successive Projection Algorithm
  (Algorithm 1 of the paper).
- `project_to_simplex()` constrained simplex projection of each
  study's function onto the recovered bases (quadratic program via
  `quadprog`).
- `fit_weight_model()` and `predict.metahunt_weight_model()` for
  Dirichlet regression of simplex weights on study-level covariates,
  with `coef.metahunt_weight_model()` for inspecting coefficients.
- `predict_target()` and `apply_wrapper()` for composing predictions
  and scalar summaries by hand.

## Data preparation

- `build_grid()` constructs a shared evaluation grid from any
  reference patient-level dataset.
- `f_hat_from_models()` evaluates a list of fitted models on the
  shared grid with class-aware dispatch for `ranger`, `grf`
  (`causal_forest`, `regression_forest`), and a default branch that
  covers `lm`/`glm`/`randomForest`. Custom S4 classes can supply
  their own `predict_fn`.

## Baselines

- `minmax_regret()` implements the covariate-free worst-case-regret
  aggregator of Zhang, Huang, and Imai (2024,
  [arXiv:2412.11136](https://arxiv.org/abs/2412.11136)).

## Documentation

- Tutorials: `metahunt-intro`, `data-prep`, `grid-weights`,
  `wrapper-scalar`, plus `get-started`.
- Companion paper: Shi, Imai, and Zhang (2024,
  [arXiv:2604.23847](https://arxiv.org/abs/2604.23847)).
