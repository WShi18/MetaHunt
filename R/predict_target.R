#' Predict the target function for new study-level covariates
#'
#' Given a fitted d-fSPA basis decomposition and a fitted weight model,
#' compute the predicted target function on the shared grid as
#' \deqn{\tilde f^{(0)}(\cdot) = \sum_{k=1}^{\hat K} \tilde\pi_{0k}\, \hat g_k(\cdot),
#'   \qquad \tilde\pi_0 = \widehat{\mathcal M}(\mathbf W_0).}
#'
#' @param dfspa_fit A `"dfspa"` object produced by [dfspa()]. Its `bases`
#'   slot (`K`-by-`G_grid` matrix) provides the basis functions on the grid.
#' @param weight_model A `"metahunt_weight_model"` object produced by
#'   [fit_weight_model()]. Must have `weight_model$K == nrow(dfspa_fit$bases)`.
#' @param W_new A matrix or data frame of study-level covariates for the new
#'   target studies, with columns matching those used to fit `weight_model`.
#'
#' @return An `nrow(W_new)`-by-`G_grid` numeric matrix; row `j` is the
#'   predicted target function on the grid for the `j`-th new study.
#'
#' @seealso [apply_wrapper()] to reduce predicted functions to scalars.
#'
#' @examples
#' set.seed(1)
#' G <- 40; m <- 60; K_true <- 3
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#'
#' # generate study-level covariates and softmax weights
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' beta <- cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
#' eta  <- as.matrix(W) %*% beta
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.02), m, G)
#'
#' fit    <- dfspa(F_hat, K = K_true)
#' pi_hat <- project_to_simplex(F_hat, fit$bases)
#' wm     <- fit_weight_model(pi_hat, W)
#'
#' W_new   <- data.frame(w1 = c(0, 1), w2 = c(0, -1))
#' f_pred  <- predict_target(fit, wm, W_new)
#' dim(f_pred)   # 2 x G
#'
#' @export
predict_target <- function(dfspa_fit, weight_model, W_new) {
  if (!inherits(dfspa_fit, "dfspa")) {
    stop("`dfspa_fit` must be a `dfspa` object (from `dfspa()`).")
  }
  if (!inherits(weight_model, "metahunt_weight_model")) {
    stop("`weight_model` must be a `metahunt_weight_model` object (from `fit_weight_model()`).")
  }
  if (weight_model$K != nrow(dfspa_fit$bases)) {
    stop(sprintf(
      "Basis/weight-model mismatch: dfspa_fit has K = %d bases but weight_model has K = %d.",
      nrow(dfspa_fit$bases), weight_model$K))
  }

  pi_new <- stats::predict(weight_model, newdata = W_new)
  pi_new %*% dfspa_fit$bases
}

#' Reduce predicted functions to scalars via a user-supplied wrapper
#'
#' Many downstream quantities of interest (average treatment effect,
#' pointwise predictions, or other functionals of \eqn{f^{(0)}}) are scalar
#' summaries of the predicted function. `apply_wrapper()` applies any
#' user-supplied reduction to each row of a function matrix, with a default
#' of the weighted mean with respect to `grid_weights` (which coincides with
#' \eqn{\int f\,d\mu} when `grid_weights` represents \eqn{\mu}).
#'
#' @param F_mat An `n`-by-`G_grid` numeric matrix; each row a function
#'   on a shared grid (e.g. the output of [predict_target()]).
#' @param wrapper Either `NULL` (default, weighted mean) or a function that
#'   takes a single numeric vector of length `G_grid` and returns a scalar.
#'   Examples: `mean`, `median`, `max`, `function(f) f[17]` (point
#'   evaluation), `function(f) sum(f^2)`.
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector.
#'   Used only when `wrapper = NULL`. Defaults to uniform `1 / G_grid`.
#'
#' @return A length-`n` numeric vector of scalar summaries.
#'
#' @examples
#' F <- matrix(1:12, nrow = 3, byrow = TRUE)  # 3 "functions" on a 4-point grid
#' apply_wrapper(F)                            # row means (uniform grid weights)
#' apply_wrapper(F, wrapper = max)             # row maxes
#' apply_wrapper(F, wrapper = function(f) f[2]) # point evaluation at grid idx 2
#'
#' @export
apply_wrapper <- function(F_mat, wrapper = NULL, grid_weights = NULL) {
  if (!is.matrix(F_mat) || !is.numeric(F_mat)) {
    stop("`F_mat` must be a numeric matrix.")
  }
  G_grid <- ncol(F_mat)

  if (is.null(wrapper)) {
    if (is.null(grid_weights)) {
      grid_weights <- rep(1 / G_grid, G_grid)
    } else {
      if (!is.numeric(grid_weights) || length(grid_weights) != G_grid ||
          any(!is.finite(grid_weights)) || any(grid_weights < 0)) {
        stop("`grid_weights` must be a length-`ncol(F_mat)` non-negative numeric vector.")
      }
      if (sum(grid_weights) <= 0) {
        stop("`grid_weights` must sum to a positive value.")
      }
    }
    return(as.numeric(F_mat %*% grid_weights) / sum(grid_weights))
  }

  if (!is.function(wrapper)) {
    stop("`wrapper` must be a function or NULL.")
  }
  vals <- apply(F_mat, 1, wrapper)
  if (!is.numeric(vals) || length(vals) != nrow(F_mat)) {
    stop("`wrapper` must return a single numeric value for each row.")
  }
  as.numeric(vals)
}
