#' MetaHunt: Privacy-Preserving Meta-Analysis via Low-Rank Basis Hunting
#'
#' Tools for privacy-preserving meta-analysis of function-valued quantities
#' (e.g. regression, conditional average treatment effect functions) across
#' heterogeneous studies. Implements the MetaHunt pipeline of Shi, Imai,
#' and Zhang: denoised functional Successive Projection Algorithm (d-fSPA)
#' for basis hunting, constrained projection of study functions onto the
#' recovered simplex, Dirichlet regression of mixing weights on study-level
#' covariates, target prediction, and split or cross conformal prediction
#' intervals.
#'
#' @section Main entry points:
#' * [metahunt()] — fit the full pipeline end-to-end.
#' * [predict.metahunt()] — predict target functions for new study-level covariates.
#' * [split_conformal()], [cross_conformal()], [conformal_from_fit()] —
#'   prediction intervals.
#' * [minmax_regret()] — covariate-free worst-case-regret aggregator
#'   (Zhang, Huang, and Imai 2024) as a baseline.
#' * [f_hat_from_models()], [build_grid()] — onramp from fitted-model
#'   lists and reference data to the package's matrix inputs.
#'
#' @section Pipeline building blocks:
#' [dfspa()], [project_to_simplex()], [fit_weight_model()],
#' [predict_target()], [apply_wrapper()],
#' [reconstruction_error_curve()], [cv_error_curve()],
#' [select_denoising_params()].
#'
#' @keywords internal
"_PACKAGE"
