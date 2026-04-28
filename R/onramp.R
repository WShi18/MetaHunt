#' Build the `F_hat` matrix from a list of fitted study-level models
#'
#' Many users arrive at MetaHunt with one fitted model per study (e.g. a
#' `ranger::ranger()` random forest or a `grf::causal_forest()`) and a
#' chosen evaluation grid. `f_hat_from_models()` evaluates each model on
#' the shared grid and stacks the predictions into the `m`-by-`G_grid`
#' matrix the rest of the package expects.
#'
#' @details
#' By default the function dispatches on the model class:
#' * `ranger` objects are evaluated as `predict(model, data = grid)$predictions`.
#' * Objects inheriting from `causal_forest` or `grf` are evaluated as
#'   `predict(model, newdata = grid)$predictions`.
#' * All other classes fall through to
#'   `as.numeric(predict(model, newdata = grid))`, which works for `lm`,
#'   `glm`, `randomForest`, and most other R model objects.
#'
#' Override the dispatch with `predict_fn = function(model, grid) ...` if
#' your models need bespoke handling. The function must return a length-`G_grid`
#' numeric vector for each model.
#'
#' All rows of the returned matrix must have the same length (`G_grid`)
#' and contain no `NA` values; the function errors otherwise.
#'
#' @param models A non-empty list of fitted model objects, one per study.
#' @param grid A data frame (or matrix) of grid points at which to evaluate
#'   each model. Same columns as the data each model was fitted on.
#' @param predict_fn Optional `function(model, grid)` returning a numeric
#'   vector; defaults to the class-aware dispatcher above.
#'
#' @return An `length(models)`-by-`nrow(grid)` numeric matrix; row `i` is
#'   model `i` evaluated at every row of `grid`.
#'
#' @seealso [build_grid()] to construct `grid` from a reference dataset.
#'
#' @examples
#' # Toy example: each "centre" fits a polynomial regression
#' set.seed(1)
#' make_centre_data <- function(slope) {
#'   x <- runif(60)
#'   data.frame(x = x, y = slope * x + rnorm(60, sd = 0.1))
#' }
#' models <- lapply(c(-1, 0, 1, 0.5, -0.5), function(s)
#'   stats::lm(y ~ poly(x, 2), data = make_centre_data(s)))
#'
#' grid  <- data.frame(x = seq(0, 1, length.out = 30))
#' F_hat <- f_hat_from_models(models, grid)
#' dim(F_hat)   # 5 x 30
#'
#' @export
f_hat_from_models <- function(models, grid, predict_fn = NULL) {
  if (!is.list(models) || length(models) < 1L) {
    stop("`models` must be a non-empty list.")
  }
  if (is.matrix(grid)) grid <- as.data.frame(grid)
  if (!is.data.frame(grid)) {
    stop("`grid` must be a matrix or data frame.")
  }
  if (nrow(grid) < 1L) stop("`grid` must have at least one row.")

  if (is.null(predict_fn)) {
    predict_fn <- .predict_dispatch
  } else if (!is.function(predict_fn)) {
    stop("`predict_fn` must be a function or NULL.")
  }

  m <- length(models)
  G_grid <- nrow(grid)
  out <- matrix(NA_real_, nrow = m, ncol = G_grid)

  for (i in seq_len(m)) {
    pred <- tryCatch(predict_fn(models[[i]], grid),
                     error = function(e)
                       stop(sprintf("predict_fn failed on model %d: %s",
                                    i, conditionMessage(e))))
    pred <- as.numeric(pred)
    if (length(pred) != G_grid) {
      stop(sprintf("Model %d returned %d predictions but grid has %d rows.",
                   i, length(pred), G_grid))
    }
    if (anyNA(pred)) {
      stop(sprintf("Model %d returned NA predictions.", i))
    }
    out[i, ] <- pred
  }
  out
}

#' @keywords internal
#' @noRd
.predict_dispatch <- function(model, grid) {
  if (inherits(model, "ranger")) {
    pred <- stats::predict(model, data = grid)$predictions
  } else if (inherits(model, "causal_forest") ||
             inherits(model, "regression_forest") ||
             inherits(model, "grf")) {
    pred <- stats::predict(model, newdata = grid)$predictions
  } else {
    pred <- stats::predict(model, newdata = grid)
  }
  as.numeric(pred)
}

#' Build a shared evaluation grid from a reference dataset
#'
#' Constructs a data frame of grid points suitable for [f_hat_from_models()]
#' from any reference patient-level dataset. This is convenient when the
#' patient-level covariate space is multidimensional and there is no
#' obvious one-dimensional grid.
#'
#' @details
#' If `n_grid` is `NULL` or at least `nrow(reference_data)`, the reference
#' data is returned unchanged. Otherwise `n_grid` rows are sampled uniformly
#' at random (without replacement). The reference data should be on the
#' same scale and have the same columns as the data each centre's model
#' was fitted on.
#'
#' The empirical distribution of the returned grid implicitly defines the
#' \eqn{\mu} measure used downstream. Pass uniform `grid_weights` (the
#' default) to weight each grid point equally; pass non-uniform
#' `grid_weights` to weight by an external reference distribution.
#'
#' @param reference_data A data frame (or matrix) of reference patient-level
#'   covariates. May be a held-out target-population sample, the pooled
#'   source covariates, or any plausible reference distribution.
#' @param n_grid Optional integer giving the desired grid size. If `NULL`
#'   or `>= nrow(reference_data)`, the full reference data is returned.
#' @param seed Optional integer seed for reproducibility (used only when
#'   sub-sampling).
#'
#' @return A data frame of grid points.
#'
#' @examples
#' set.seed(1)
#' ref <- data.frame(age = rnorm(500, 60, 10),
#'                   bp  = rnorm(500, 130, 15),
#'                   sex = sample(c("F", "M"), 500, replace = TRUE))
#' grid <- build_grid(ref, n_grid = 50, seed = 1)
#' nrow(grid)
#' head(grid)
#'
#' @export
build_grid <- function(reference_data, n_grid = NULL, seed = NULL) {
  if (is.matrix(reference_data)) {
    reference_data <- as.data.frame(reference_data)
  }
  if (!is.data.frame(reference_data)) {
    stop("`reference_data` must be a matrix or data frame.")
  }
  N <- nrow(reference_data)
  if (N < 1L) stop("`reference_data` is empty.")

  if (is.null(n_grid) || n_grid >= N) {
    return(reference_data)
  }
  if (!is.numeric(n_grid) || length(n_grid) != 1L ||
      n_grid < 1L || n_grid != as.integer(n_grid)) {
    stop("`n_grid` must be a positive integer or NULL.")
  }
  if (!is.null(seed)) withr::local_seed(seed)
  reference_data[sample(N, n_grid), , drop = FALSE]
}
