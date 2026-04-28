#' Fit the full MetaHunt pipeline
#'
#' End-to-end convenience wrapper that runs the three training-time steps
#' of the MetaHunt pipeline in sequence:
#' (1) [dfspa()] for basis hunting,
#' (2) [project_to_simplex()] for per-study weight recovery, and
#' (3) [fit_weight_model()] for modelling the weight-to-covariate map.
#' The result supports [predict.metahunt()] for generating target-function
#' predictions on new study-level covariates.
#'
#' For uncertainty quantification, pair a `"metahunt"` fit with
#' [conformal_from_fit()] (requires a separate calibration set). The
#' high-level [split_conformal()] and [cross_conformal()] functions
#' perform their own fitting and do not consume a pre-fit `"metahunt"`
#' object.
#'
#' @param F_hat An `m`-by-`G_grid` numeric matrix of study-level function
#'   evaluations on a shared grid; row `i` is \eqn{\hat f^{(i)}}.
#' @param W An `m`-by-`p` matrix or data frame of study-level covariates.
#' @param K Integer number of basis functions.
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector
#'   defining the \eqn{L^2(\mu)} inner product; defaults to uniform.
#' @param dfspa_args Named list of extra arguments for [dfspa()].
#' @param weight_model_args Named list of extra arguments for
#'   [fit_weight_model()].
#'
#' @return An object of class `"metahunt"`: a list with the `dfspa_fit`,
#'   `weight_model`, training `pi_hat`, `K`, and a stored copy of
#'   `grid_weights`.
#'
#' @seealso [predict.metahunt()], [split_conformal()],
#'   [cross_conformal()], [conformal_from_fit()],
#'   [reconstruction_error_curve()], [cv_error_curve()].
#'
#' @examples
#' set.seed(1)
#' G <- 40; m <- 80
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#'
#' fit <- metahunt(F_hat, W, K = 3)
#' fit
#' f_pred <- predict(fit, newdata = W[1:3, ])
#' dim(f_pred)                                  # 3 x G: predicted functions
#' predict(fit, newdata = W[1:3, ], wrapper = mean)  # scalar summaries
#'
#' @export
metahunt <- function(F_hat, W, K,
                     grid_weights = NULL,
                     dfspa_args = list(),
                     weight_model_args = list()) {
  call <- match.call()

  if (!is.matrix(F_hat) || !is.numeric(F_hat)) {
    stop("`F_hat` must be a numeric matrix with studies as rows.")
  }
  if (is.matrix(W)) W <- as.data.frame(W)
  if (!is.data.frame(W)) stop("`W` must be a matrix or data frame.")
  if (nrow(W) != nrow(F_hat)) {
    stop("`W` must have the same number of rows as `F_hat`.")
  }
  if (!is.numeric(K) || length(K) != 1L || K < 1L || K != as.integer(K)) {
    stop("`K` must be a positive integer.")
  }
  K <- as.integer(K)
  if (!is.list(dfspa_args) || !is.list(weight_model_args)) {
    stop("`dfspa_args` and `weight_model_args` must be named lists.")
  }

  fit <- do.call(dfspa,
                 c(list(F_hat = F_hat, K = K,
                        grid_weights = grid_weights), dfspa_args))
  pi_hat <- project_to_simplex(F_hat, fit$bases,
                               grid_weights = grid_weights)
  wm <- do.call(fit_weight_model,
                c(list(pi_hat = pi_hat, W = W), weight_model_args))

  structure(
    list(
      dfspa_fit    = fit,
      weight_model = wm,
      pi_hat       = pi_hat,
      K            = K,
      grid_weights = grid_weights,
      m            = nrow(F_hat),
      G_grid       = ncol(F_hat),
      call         = call
    ),
    class = "metahunt"
  )
}

#' Predict target functions (or scalar summaries) from a MetaHunt fit
#'
#' @param object A `"metahunt"` object from [metahunt()].
#' @param newdata A matrix or data frame of new study-level covariates.
#' @param wrapper Optional reduction function. If `NULL`, returns the full
#'   predicted function matrix (`nrow(newdata)`-by-`G_grid`). If a function,
#'   applied to each predicted function to return a scalar per new target;
#'   see [apply_wrapper()].
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector
#'   used only when `wrapper = NULL` with the default weighted-mean path.
#'   Defaults to the `grid_weights` stored in `object`.
#' @param ... Ignored.
#'
#' @return Either an `nrow(newdata)`-by-`G_grid` matrix of predicted
#'   functions (when `wrapper = NULL`) or a length-`nrow(newdata)` numeric
#'   vector of scalar summaries.
#'
#' @examples
#' set.seed(1)
#' G <- 25; m <- 40
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#'
#' fit <- metahunt(F_hat, W, K = 3, dfspa_args = list(denoise = FALSE))
#' f_pred <- predict(fit, newdata = W[1:3, ])
#' dim(f_pred)                                       # 3 x G
#' predict(fit, newdata = W[1:3, ], wrapper = mean)  # scalar summaries
#'
#' @export
predict.metahunt <- function(object, newdata, wrapper = NULL,
                             grid_weights = NULL, ...) {
  f_pred <- predict_target(object$dfspa_fit, object$weight_model, newdata)
  if (is.null(wrapper)) return(f_pred)
  gw <- if (is.null(grid_weights)) object$grid_weights else grid_weights
  apply_wrapper(f_pred, wrapper = wrapper, grid_weights = gw)
}

#' Plot recovered basis functions from a MetaHunt fit
#'
#' @param x A `"metahunt"` object from [metahunt()].
#' @param x_axis Optional numeric vector of length `G_grid` giving the
#'   x-axis values. Defaults to `seq_len(G_grid)`.
#' @param ... Passed to [graphics::matplot()].
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' set.seed(1)
#' G <- 25; m <- 40
#' x_grid <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x_grid), cos(pi * x_grid), x_grid)
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#'
#' fit <- metahunt(F_hat, W, K = 3, dfspa_args = list(denoise = FALSE))
#' plot(fit)
#' plot(fit, x_axis = x_grid)
#'
#' @export
plot.metahunt <- function(x, x_axis = NULL, ...) {
  bases <- x$dfspa_fit$bases
  K <- nrow(bases); G <- ncol(bases)
  if (is.null(x_axis)) x_axis <- seq_len(G)
  if (length(x_axis) != G) {
    stop("`x_axis` must have length equal to ncol(bases).")
  }
  args <- list(x = x_axis, y = t(bases), type = "l", lty = 1,
               xlab = "grid", ylab = expression(g[k](x)),
               main = sprintf("Recovered basis functions (K = %d)", K))
  user <- list(...)
  args[names(user)] <- user
  do.call(graphics::matplot, args)
  graphics::legend("topright",
                   legend = paste("basis", seq_len(K)),
                   col = seq_len(K), lty = 1, bty = "n")
  invisible(x)
}

#' @export
print.metahunt <- function(x, ...) {
  cat("MetaHunt fit\n")
  cat("  m (studies):   ", x$m, "\n")
  cat("  G (grid size): ", x$G_grid, "\n")
  cat("  K (bases):     ", x$K, "\n")
  cat("  weight method: ", x$weight_model$method, "\n")
  cat("  predictors:    ",
      paste(x$weight_model$predictor_names, collapse = ", "), "\n")
  invisible(x)
}

#' Summarise a MetaHunt fit
#'
#' Produces a compact summary of a `"metahunt"` object, including study/grid
#' sizes, the weight-model method, per-basis summary statistics for the
#' training simplex weights `pi_hat`, and denoising bookkeeping from the
#' underlying [dfspa()] fit.
#'
#' @param object A `"metahunt"` object from [metahunt()].
#' @param ... Ignored.
#'
#' @return An object of class `"summary.metahunt"`: a list with components
#'   \describe{
#'     \item{`m`}{Number of studies.}
#'     \item{`G_grid`}{Grid size.}
#'     \item{`K`}{Number of basis functions.}
#'     \item{`weight_method`}{Method used by the weight model.}
#'     \item{`predictor_names`}{Character vector of covariate names.}
#'     \item{`pi_summary`}{A `K`-by-5 numeric matrix; each row gives `min`,
#'       `mean`, `median`, `max`, and `sd` of the corresponding column of
#'       `object$pi_hat`.}
#'     \item{`n_kept`}{Number of studies retained after denoising.}
#'     \item{`n_dropped`}{Number of studies dropped (`m - n_kept`).}
#'     \item{`denoising`}{List with `N` and `Delta` from the `dfspa` fit
#'       (both `NA` when `denoise = FALSE`).}
#'   }
#'
#' @examples
#' set.seed(1)
#' G <- 25; m <- 40
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#'
#' fit <- metahunt(F_hat, W, K = 3, dfspa_args = list(denoise = FALSE))
#' summary(fit)
#'
#' @export
summary.metahunt <- function(object, ...) {
  if (!inherits(object, "metahunt")) {
    stop("`object` must be a `metahunt` object.")
  }

  pi_hat <- object$pi_hat
  K <- object$K
  pi_summary <- matrix(NA_real_, nrow = K, ncol = 5L,
                       dimnames = list(paste0("basis_", seq_len(K)),
                                       c("min", "mean", "median",
                                         "max", "sd")))
  for (k in seq_len(K)) {
    col_k <- pi_hat[, k]
    pi_summary[k, ] <- c(min(col_k), mean(col_k), stats::median(col_k),
                         max(col_k), stats::sd(col_k))
  }

  n_kept <- length(object$dfspa_fit$kept)
  n_dropped <- object$m - n_kept

  denoising <- list(
    N     = object$dfspa_fit$N,
    Delta = object$dfspa_fit$Delta
  )

  out <- list(
    m               = object$m,
    G_grid          = object$G_grid,
    K               = K,
    weight_method   = object$weight_model$method,
    predictor_names = object$weight_model$predictor_names,
    pi_summary      = pi_summary,
    n_kept          = n_kept,
    n_dropped       = n_dropped,
    denoising       = denoising
  )
  class(out) <- "summary.metahunt"
  out
}

#' Print a `summary.metahunt` object
#'
#' @param x A `"summary.metahunt"` object from [summary.metahunt()].
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @examples
#' set.seed(1)
#' G <- 25; m <- 40
#' x_grid <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x_grid), cos(pi * x_grid), x_grid)
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#'
#' fit <- metahunt(F_hat, W, K = 3, dfspa_args = list(denoise = FALSE))
#' print(summary(fit))
#'
#' @export
print.summary.metahunt <- function(x, ...) {
  cat("MetaHunt fit summary\n")
  cat("  m (studies):   ", x$m, "\n")
  cat("  G (grid size): ", x$G_grid, "\n")
  cat("  K (bases):     ", x$K, "\n")
  cat("  weight method: ", x$weight_method, "\n")
  cat("  predictors:    ",
      paste(x$predictor_names, collapse = ", "), "\n")
  cat("  studies kept:  ", x$n_kept, "\n")
  cat("  studies dropped:", x$n_dropped, "\n")
  N_val     <- x$denoising$N
  Delta_val <- x$denoising$Delta
  cat("  denoising N:   ",
      if (is.finite(N_val)) signif(N_val, 4) else "NA", "\n")
  cat("  denoising Delta:",
      if (is.finite(Delta_val)) signif(Delta_val, 4) else "NA", "\n")
  cat("\nPer-basis pi_hat summary:\n")
  print(round(x$pi_summary, 4))
  invisible(x)
}
