# MetaHunt: Privacy-Preserving Meta-Analysis via Low-Rank Basis Hunting

Tools for privacy-preserving meta-analysis of function-valued quantities
(e.g. regression, conditional average treatment effect functions) across
heterogeneous studies. Implements the MetaHunt pipeline of Shi, Imai,
and Zhang: denoised functional Successive Projection Algorithm (d-fSPA)
for basis hunting, constrained projection of study functions onto the
recovered simplex, Dirichlet regression of mixing weights on study-level
covariates, target prediction, and split or cross conformal prediction
intervals.

## Main entry points

- [`metahunt()`](https://wshi18.github.io/MetaHunt/reference/metahunt.md)
  — fit the full pipeline end-to-end.

- [`predict.metahunt()`](https://wshi18.github.io/MetaHunt/reference/predict.metahunt.md)
  — predict target functions for new study-level covariates.

- [`split_conformal()`](https://wshi18.github.io/MetaHunt/reference/split_conformal.md),
  [`cross_conformal()`](https://wshi18.github.io/MetaHunt/reference/cross_conformal.md),
  [`conformal_from_fit()`](https://wshi18.github.io/MetaHunt/reference/conformal_from_fit.md)
  — prediction intervals.

- [`minmax_regret()`](https://wshi18.github.io/MetaHunt/reference/minmax_regret.md)
  — covariate-free worst-case-regret aggregator (Zhang, Huang, and
  Imai 2024) as a baseline.

- [`f_hat_from_models()`](https://wshi18.github.io/MetaHunt/reference/f_hat_from_models.md),
  [`build_grid()`](https://wshi18.github.io/MetaHunt/reference/build_grid.md)
  — onramp from fitted-model lists and reference data to the package's
  matrix inputs.

## Pipeline building blocks

[`dfspa()`](https://wshi18.github.io/MetaHunt/reference/dfspa.md),
[`project_to_simplex()`](https://wshi18.github.io/MetaHunt/reference/project_to_simplex.md),
[`fit_weight_model()`](https://wshi18.github.io/MetaHunt/reference/fit_weight_model.md),
[`predict_target()`](https://wshi18.github.io/MetaHunt/reference/predict_target.md),
[`apply_wrapper()`](https://wshi18.github.io/MetaHunt/reference/apply_wrapper.md),
[`reconstruction_error_curve()`](https://wshi18.github.io/MetaHunt/reference/reconstruction_error_curve.md),
[`cv_error_curve()`](https://wshi18.github.io/MetaHunt/reference/cv_error_curve.md),
[`select_denoising_params()`](https://wshi18.github.io/MetaHunt/reference/select_denoising_params.md).

## See also

Useful links:

- <https://github.com/WShi18/MetaHunt>

- <https://wshi18.github.io/MetaHunt/>

- <https://arxiv.org/abs/2604.23847>

- Report bugs at <https://github.com/WShi18/MetaHunt/issues>

## Author

**Maintainer**: Wenqi Shi <wenqishi18@gmail.com>

Authors:

- Kosuke Imai <imai@harvard.edu>

- Yi Zhang <yizhang0017@gmail.com>
