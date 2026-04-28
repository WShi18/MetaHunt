#' @keywords internal
#' @noRd
.conformal_quantile <- function(scores, alpha) {
  scores <- scores[is.finite(scores)]
  n <- length(scores)
  if (n == 0L) return(Inf)
  k <- ceiling((1 - alpha) * (n + 1))
  if (k > n) return(Inf)
  sort(scores)[k]
}

#' @keywords internal
#' @noRd
.run_pipeline <- function(F_hat_tr, W_tr, K, grid_weights,
                          dfspa_args, weight_model_args) {
  fit <- do.call(dfspa,
                 c(list(F_hat = F_hat_tr, K = K,
                        grid_weights = grid_weights), dfspa_args))
  pi_tr <- project_to_simplex(F_hat_tr, fit$bases,
                              grid_weights = grid_weights)
  wm <- do.call(fit_weight_model,
                c(list(pi_hat = pi_tr, W = W_tr), weight_model_args))
  list(fit = fit, wm = wm)
}

#' Split conformal intervals from a pre-fit MetaHunt pipeline
#'
#' Lower-level entry point that builds split conformal intervals for new
#' target covariates using an already-fitted d-fSPA basis decomposition,
#' an already-fitted weight model, and a user-supplied calibration set.
#' Use this when you have independently tuned `K` or want to reuse a
#' pipeline fit; otherwise the high-level [split_conformal()] is usually
#' more convenient.
#'
#' @param dfspa_fit A `"dfspa"` object from [dfspa()].
#' @param weight_model A `"metahunt_weight_model"` object from
#'   [fit_weight_model()]. Must satisfy
#'   `weight_model$K == nrow(dfspa_fit$bases)`.
#' @param F_cal An `n_cal`-by-`G_grid` numeric matrix of observed study-level
#'   function evaluations for the calibration set. Calibration studies must
#'   be independent of both `dfspa_fit` and `weight_model` for valid
#'   coverage.
#' @param W_cal A matrix or data frame of study-level covariates for the
#'   calibration set, with the same columns used to fit `weight_model`.
#' @param W_new A matrix or data frame of study-level covariates for new
#'   target studies.
#' @param alpha Miscoverage level; default `0.05`.
#' @param wrapper Optional reduction function (see [apply_wrapper()]). If
#'   `NULL`, intervals are pointwise at every grid point.
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector
#'   used for the default wrapper and for weighted inner products.
#'
#' @return An object of class `"metahunt_conformal"`; see [split_conformal()]
#'   for a description of its fields.
#'
#' @seealso [split_conformal()] for the high-level version that splits and
#'   fits internally, [cross_conformal()] for the K-fold variant.
#'
#' @examples
#' set.seed(1)
#' G <- 30; m <- 80
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m))
#' eta <- cbind(0.8 * W$w1, -0.4 * W$w1, rep(0, m))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#'
#' # user-controlled split and fit
#' tr  <- 1:50; cal <- 51:70; new <- 71:80
#' fit <- dfspa(F_hat[tr, ], K = 3)
#' pih <- project_to_simplex(F_hat[tr, ], fit$bases)
#' wm  <- fit_weight_model(pih, W[tr, , drop = FALSE])
#'
#' res <- conformal_from_fit(
#'   fit, wm,
#'   F_cal = F_hat[cal, ], W_cal = W[cal, , drop = FALSE],
#'   W_new = W[new, , drop = FALSE],
#'   wrapper = mean
#' )
#' res
#'
#' @export
conformal_from_fit <- function(dfspa_fit, weight_model,
                               F_cal, W_cal, W_new,
                               alpha = 0.05,
                               wrapper = NULL,
                               grid_weights = NULL) {
  if (!inherits(dfspa_fit, "dfspa")) {
    stop("`dfspa_fit` must be a `dfspa` object (from `dfspa()`).")
  }
  if (!inherits(weight_model, "metahunt_weight_model")) {
    stop("`weight_model` must be a `metahunt_weight_model` object.")
  }
  if (weight_model$K != nrow(dfspa_fit$bases)) {
    stop(sprintf(
      "Basis/weight-model mismatch: dfspa_fit has K = %d bases but weight_model has K = %d.",
      nrow(dfspa_fit$bases), weight_model$K))
  }
  if (!is.matrix(F_cal) || !is.numeric(F_cal)) {
    stop("`F_cal` must be a numeric matrix with calibration studies as rows.")
  }
  if (ncol(F_cal) != ncol(dfspa_fit$bases)) {
    stop("`F_cal` must have the same number of columns as `dfspa_fit$bases`.")
  }
  if (is.matrix(W_cal)) W_cal <- as.data.frame(W_cal)
  if (!is.data.frame(W_cal)) stop("`W_cal` must be a matrix or data frame.")
  if (nrow(W_cal) != nrow(F_cal)) {
    stop("`W_cal` and `F_cal` must have the same number of rows.")
  }
  if (nrow(F_cal) < 2L) {
    stop("Calibration set must have at least 2 studies.")
  }
  if (is.matrix(W_new)) W_new <- as.data.frame(W_new)
  if (!is.data.frame(W_new)) stop("`W_new` must be a matrix or data frame.")
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).")
  }

  pred_cal <- predict_target(dfspa_fit, weight_model, W_cal)

  .build_conformal_output(
    obs_cal = F_cal, pred_cal = pred_cal,
    pipe = list(fit = dfspa_fit, wm = weight_model),
    W_new = W_new,
    alpha = alpha, wrapper = wrapper, grid_weights = grid_weights,
    method = "split",
    n_cal_label = nrow(F_cal)
  )
}

#' Split conformal prediction intervals for target-function predictions
#'
#' Implements Algorithm 2 of the paper (split conformal prediction) over
#' the MetaHunt pipeline. Studies are partitioned into training and
#' calibration sets. The training set is used to fit d-fSPA, the constrained
#' projection, and the weight model; the calibration set supplies
#' conformity scores, which determine the width of the intervals.
#'
#' @details
#' Given a target function, one can either construct intervals **pointwise
#' at every grid point** (when `wrapper = NULL`) or for a **scalar summary**
#' of the target function (when `wrapper` is a function).
#'
#' * **Pointwise** (`wrapper = NULL`): for each grid point `g` the
#'   conformity score is \eqn{R_{i,g} = |\hat f^{(i)}(x_g) - \tilde
#'   f^{(i)}(x_g)|} across calibration studies `i`. A separate
#'   \eqn{(1-\alpha)}-quantile \eqn{q_g} is computed per grid point, and
#'   the interval at grid point `g` for target `j` is
#'   \eqn{[\tilde f^{(j)}(x_g)-q_g,\ \tilde f^{(j)}(x_g)+q_g]}.
#' * **Scalar** (`wrapper` supplied): conformity scores are
#'   \eqn{R_i = |s(\hat f^{(i)}) - s(\tilde f^{(i)})|} with
#'   `s = wrapper`, and the interval for each target is
#'   \eqn{[s(\tilde f^{(j)})-q,\ s(\tilde f^{(j)})+q]} with a single
#'   shared quantile \eqn{q}.
#'
#' The finite-sample quantile is
#' \eqn{q = R_{(k)}} with \eqn{k = \lceil(1-\alpha)(n_\mathrm{cal}+1)\rceil};
#' if \eqn{k > n_\mathrm{cal}}, `q = Inf` and intervals are
#' \eqn{(-\infty, \infty)}.
#'
#' @param F_hat An `m`-by-`G_grid` numeric matrix of study-level function
#'   evaluations.
#' @param W An `m`-by-`p` matrix or data frame of study-level covariates.
#' @param W_new A matrix or data frame of new target covariates. Must
#'   contain columns matching `W`.
#' @param K Integer number of basis functions.
#' @param alpha Miscoverage level; interval has nominal coverage
#'   \eqn{1-\alpha}. Default `0.05`.
#' @param cal_frac Numeric in `(0, 1)` giving the fraction of studies in
#'   the calibration set. Default `0.3`. Ignored if `calibration_idx` is
#'   supplied.
#' @param wrapper Optional reduction function (see [apply_wrapper()]). If
#'   `NULL`, intervals are constructed pointwise at every grid point.
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector
#'   used for the \eqn{L^2(\mu)} norm and for the default wrapper.
#' @param calibration_idx Optional integer vector of row indices in `F_hat`
#'   to use as the calibration set. If supplied, `cal_frac` is ignored.
#' @param dfspa_args,weight_model_args Named lists passed to [dfspa()] and
#'   [fit_weight_model()] respectively.
#' @param seed Optional integer seed for reproducible train/calibration
#'   splits.
#'
#' @return An object of class `"metahunt_conformal"`: a list with
#'   \describe{
#'     \item{`prediction`}{Point predictions for `W_new`. A numeric vector
#'       of length `nrow(W_new)` in the scalar case, or an
#'       `nrow(W_new)`-by-`G_grid` matrix in the pointwise case.}
#'     \item{`lower`, `upper`}{Interval endpoints, same shape as
#'       `prediction`.}
#'     \item{`alpha`}{Miscoverage level used.}
#'     \item{`method`}{`"split"`.}
#'     \item{`n_cal`}{Calibration sample size.}
#'     \item{`quantile`}{The conformal quantile: a scalar (scalar case)
#'       or a length-`G_grid` vector (pointwise case).}
#'     \item{`wrapper`}{The wrapper used, or `NULL`.}
#'   }
#'
#' @examples
#' set.seed(1)
#' G <- 40; m <- 80; K_true <- 3
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#'
#' W_new <- data.frame(w1 = c(0, 1), w2 = c(0, -1))
#' # pointwise intervals at every grid point
#' pi_grid <- split_conformal(F_hat, W, W_new, K = 3, seed = 1)
#' dim(pi_grid$lower)  # 2 x 40
#' # scalar intervals for the grid-weighted mean (ATE-style)
#' pi_ate  <- split_conformal(F_hat, W, W_new, K = 3, wrapper = mean, seed = 1)
#' pi_ate$prediction
#'
#' @export
split_conformal <- function(F_hat, W, W_new, K,
                            alpha = 0.05, cal_frac = 0.3,
                            wrapper = NULL, grid_weights = NULL,
                            calibration_idx = NULL,
                            dfspa_args = list(),
                            weight_model_args = list(),
                            seed = NULL) {
  .validate_conformal_inputs(F_hat, W, W_new, K, alpha, grid_weights,
                             dfspa_args, weight_model_args)
  m <- nrow(F_hat)
  if (is.matrix(W))     W     <- as.data.frame(W)
  if (is.matrix(W_new)) W_new <- as.data.frame(W_new)

  if (is.null(calibration_idx)) {
    if (!is.numeric(cal_frac) || length(cal_frac) != 1L ||
        cal_frac <= 0 || cal_frac >= 1) {
      stop("`cal_frac` must be a single number in (0, 1).")
    }
    if (!is.null(seed)) withr::local_seed(seed)
    n_cal <- max(1L, floor(cal_frac * m))
    calibration_idx <- sort(sample(seq_len(m), n_cal))
  } else {
    calibration_idx <- as.integer(calibration_idx)
    if (any(calibration_idx < 1L) || any(calibration_idx > m) ||
        anyDuplicated(calibration_idx)) {
      stop("`calibration_idx` must be unique integers in 1:nrow(F_hat).")
    }
  }
  training_idx <- setdiff(seq_len(m), calibration_idx)
  if (length(training_idx) < K) {
    stop(sprintf("Training set has only %d studies but K = %d.",
                 length(training_idx), K))
  }
  if (length(calibration_idx) < 2L) {
    stop("Calibration set must have at least 2 studies.")
  }

  pipe <- .run_pipeline(F_hat[training_idx, , drop = FALSE],
                        W[training_idx, , drop = FALSE],
                        K, grid_weights, dfspa_args, weight_model_args)

  conformal_from_fit(
    dfspa_fit    = pipe$fit,
    weight_model = pipe$wm,
    F_cal        = F_hat[calibration_idx, , drop = FALSE],
    W_cal        = W[calibration_idx, , drop = FALSE],
    W_new        = W_new,
    alpha        = alpha,
    wrapper      = wrapper,
    grid_weights = grid_weights
  )
}

#' Cross-conformal prediction intervals (pooled K-fold scores)
#'
#' Computes K-fold split conformal intervals in which the calibration scores
#' are pooled across folds, while point predictions for new targets are
#' produced from a final pipeline fit on all studies. Equivalent to running
#' [split_conformal()] `n_folds` times with different calibration sets and
#' pooling all conformity scores into a single empirical distribution.
#'
#' @details
#' For each fold, the MetaHunt pipeline is fit on the out-of-fold studies,
#' and conformity scores are computed on the in-fold studies. After all
#' folds complete, the pooled scores yield a single
#' \eqn{(1-\alpha)}-quantile (or one per grid point). Point predictions for
#' `W_new` use a pipeline refit on the full dataset. The interval at grid
#' point `g` for target `j` is
#' \eqn{[\tilde f^{(j)}(x_g) - q_g, \tilde f^{(j)}(x_g) + q_g]} (or the
#' scalar analogue when `wrapper` is supplied).
#'
#' This differs from Vovk's original cross-conformal predictor for
#' classification. For regression, pooling scores across folds is a
#' common practical extension of split conformal and reduces the variance
#' due to the single calibration split. Exact finite-sample coverage is
#' not guaranteed; see Barber et al. (2021, Jackknife+) for more
#' conservative alternatives.
#'
#' @inheritParams split_conformal
#' @param n_folds Integer number of folds (>= 2). Default `5`.
#'
#' @return An object of class `"metahunt_conformal"` (see
#'   [split_conformal()] for fields). `method` is `"cross"` and `n_cal` is
#'   the number of pooled scores.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' G <- 30; m <- 60; K_true <- 3
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m))
#' eta <- cbind(0.8 * W$w1, -0.3 * W$w1, rep(0, m))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#' W_new <- data.frame(w1 = c(0, 1))
#'
#' res <- cross_conformal(F_hat, W, W_new, K = 3, n_folds = 4,
#'                        dfspa_args = list(denoise = FALSE), seed = 1)
#' res
#' }
#'
#' @export
cross_conformal <- function(F_hat, W, W_new, K,
                            alpha = 0.05, n_folds = 5L,
                            wrapper = NULL, grid_weights = NULL,
                            dfspa_args = list(),
                            weight_model_args = list(),
                            seed = NULL) {
  .validate_conformal_inputs(F_hat, W, W_new, K, alpha, grid_weights,
                             dfspa_args, weight_model_args)
  m <- nrow(F_hat); G_grid <- ncol(F_hat)
  if (is.matrix(W))     W     <- as.data.frame(W)
  if (is.matrix(W_new)) W_new <- as.data.frame(W_new)

  if (!is.numeric(n_folds) || length(n_folds) != 1L ||
      n_folds < 2L || n_folds != as.integer(n_folds)) {
    stop("`n_folds` must be an integer >= 2.")
  }
  n_folds <- as.integer(n_folds)
  if (n_folds > m) stop("`n_folds` cannot exceed the number of studies.")

  if (!is.null(seed)) withr::local_seed(seed)
  fold_id <- sample(rep(seq_len(n_folds), length.out = m))

  if (is.null(wrapper)) {
    pooled_scores <- matrix(NA_real_, m, G_grid)
  } else {
    pooled_scores <- rep(NA_real_, m)
  }
  failures <- character(0)
  for (f in seq_len(n_folds)) {
    train <- which(fold_id != f)
    test  <- which(fold_id == f)
    result <- tryCatch({
      pipe_f <- .run_pipeline(F_hat[train, , drop = FALSE],
                              W[train, , drop = FALSE],
                              K, grid_weights, dfspa_args, weight_model_args)
      pred_f <- predict_target(pipe_f$fit, pipe_f$wm,
                               W[test, , drop = FALSE])
      if (is.null(wrapper)) {
        abs(F_hat[test, , drop = FALSE] - pred_f)
      } else {
        abs(apply_wrapper(F_hat[test, , drop = FALSE], wrapper, grid_weights) -
              apply_wrapper(pred_f, wrapper, grid_weights))
      }
    }, error = function(e) e)
    if (inherits(result, "error")) {
      failures <- c(failures, sprintf("fold %d: %s", f,
                                      conditionMessage(result)))
      next
    }
    if (is.null(wrapper)) {
      pooled_scores[test, ] <- result
    } else {
      pooled_scores[test]   <- result
    }
  }
  if (length(failures)) {
    warning(sprintf("%d of %d folds failed:\n%s",
                    length(failures), n_folds,
                    paste0("  ", failures, collapse = "\n")))
  }

  pipe_all <- .run_pipeline(F_hat, W, K, grid_weights,
                            dfspa_args, weight_model_args)

  # Fake observed/predicted matrices to reuse output builder; we already
  # have the scores, so we'll inject them directly.
  .build_conformal_output(
    obs_cal = NULL, pred_cal = NULL,
    pipe = pipe_all, W_new = W_new,
    alpha = alpha, wrapper = wrapper, grid_weights = grid_weights,
    method = "cross",
    n_cal_label = NA_integer_,
    precomputed_scores = pooled_scores
  )
}

#' @keywords internal
#' @noRd
.validate_conformal_inputs <- function(F_hat, W, W_new, K, alpha,
                                       grid_weights,
                                       dfspa_args, weight_model_args) {
  if (!is.matrix(F_hat) || !is.numeric(F_hat)) {
    stop("`F_hat` must be a numeric matrix with studies as rows.")
  }
  if (!is.data.frame(W) && !is.matrix(W)) {
    stop("`W` must be a matrix or data frame.")
  }
  if (!is.data.frame(W_new) && !is.matrix(W_new)) {
    stop("`W_new` must be a matrix or data frame.")
  }
  if (nrow(W) != nrow(F_hat)) {
    stop("`W` must have the same number of rows as `F_hat`.")
  }
  if (!is.numeric(K) || length(K) != 1L || K < 1L || K != as.integer(K)) {
    stop("`K` must be a positive integer.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number in (0, 1).")
  }
  if (!is.list(dfspa_args) || !is.list(weight_model_args)) {
    stop("`dfspa_args` and `weight_model_args` must be named lists.")
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.build_conformal_output <- function(obs_cal, pred_cal, pipe, W_new, alpha,
                                    wrapper, grid_weights, method,
                                    n_cal_label,
                                    precomputed_scores = NULL) {
  if (is.null(precomputed_scores)) {
    if (is.null(wrapper)) {
      scores <- abs(obs_cal - pred_cal)
    } else {
      scores <- abs(apply_wrapper(obs_cal, wrapper, grid_weights) -
                      apply_wrapper(pred_cal, wrapper, grid_weights))
    }
  } else {
    scores <- precomputed_scores
  }

  if (is.null(wrapper)) {
    q <- apply(scores, 2, .conformal_quantile, alpha = alpha)
    n_cal_used <- sum(stats::complete.cases(scores))
    if (any(is.infinite(q))) {
      n_inf <- sum(is.infinite(q))
      warning(sprintf(
        paste0("With n_cal = %d and alpha = %g, the conformal quantile is infinite ",
               "at %d of %d grid points; intervals are unbounded there. Increase ",
               "calibration size (raise `cal_frac` or supply more studies) or use ",
               "a larger `alpha`."),
        n_cal_used, alpha, n_inf, length(q)))
    }
    pred_new <- predict_target(pipe$fit, pipe$wm, W_new)
    lower <- sweep(pred_new, 2, q, `-`)
    upper <- sweep(pred_new, 2, q, `+`)
    result <- list(
      prediction = pred_new, lower = lower, upper = upper,
      alpha = alpha, method = method,
      n_cal = if (is.na(n_cal_label)) n_cal_used else n_cal_label,
      quantile = q, wrapper = NULL
    )
  } else {
    scores_vec <- as.numeric(scores)
    q <- .conformal_quantile(scores_vec, alpha)
    n_cal_used <- sum(is.finite(scores_vec))
    if (is.infinite(q)) {
      warning(sprintf(
        paste0("With n_cal = %d and alpha = %g, the conformal quantile is infinite; ",
               "intervals are unbounded. Increase calibration size or use a larger `alpha`."),
        n_cal_used, alpha))
    }
    pred_new <- predict_target(pipe$fit, pipe$wm, W_new)
    s_new <- apply_wrapper(pred_new, wrapper, grid_weights)
    result <- list(
      prediction = s_new, lower = s_new - q, upper = s_new + q,
      alpha = alpha, method = method,
      n_cal = if (is.na(n_cal_label)) n_cal_used else n_cal_label,
      quantile = q, wrapper = wrapper
    )
  }
  class(result) <- "metahunt_conformal"
  result
}

#' Plot a conformal prediction-interval object
#'
#' For pointwise objects (no `wrapper` was used), draws the predicted
#' function for one target study together with its pointwise band.
#' For scalar objects, draws point predictions with whisker error bars
#' over all targets.
#'
#' @param x A `"metahunt_conformal"` object.
#' @param target_idx For pointwise objects, the integer index of the
#'   target study to plot (default `1`). Ignored for scalar objects.
#' @param x_axis Optional numeric vector of length `G_grid` giving the
#'   x-axis values for pointwise plotting. Defaults to `seq_len(G_grid)`.
#' @param fill Polygon fill colour for the band. Default semi-transparent
#'   steel blue.
#' @param line_col Line colour for the predicted function (or points in
#'   scalar mode). Default steel blue.
#' @param ... Additional graphical parameters passed to the underlying
#'   plotting calls.
#'
#' @return Invisibly returns `x`.
#'
#' @export
plot.metahunt_conformal <- function(x, target_idx = 1L, x_axis = NULL,
                                    fill = grDevices::adjustcolor("steelblue", 0.2),
                                    line_col = "steelblue", ...) {
  alpha <- x$alpha
  if (is.matrix(x$prediction)) {
    if (target_idx < 1L || target_idx > nrow(x$prediction)) {
      stop("`target_idx` is out of range for this object.")
    }
    G <- ncol(x$prediction)
    if (is.null(x_axis)) x_axis <- seq_len(G)
    if (length(x_axis) != G) {
      stop("`x_axis` must have length equal to ncol(x$prediction).")
    }
    pred <- x$prediction[target_idx, ]
    lo   <- x$lower[target_idx, ]
    up   <- x$upper[target_idx, ]
    finite_y <- c(pred, lo[is.finite(lo)], up[is.finite(up)])
    plot(x_axis, pred, type = "n", ylim = range(finite_y),
         xlab = "grid", ylab = expression(tilde(f)(x)),
         main = sprintf("%.0f%% pointwise band (target %d)",
                        100 * (1 - alpha), target_idx), ...)
    if (all(is.finite(c(lo, up)))) {
      graphics::polygon(c(x_axis, rev(x_axis)), c(lo, rev(up)),
                        col = fill, border = NA)
    }
    graphics::lines(x_axis, pred, lwd = 2, col = line_col)
  } else {
    n <- length(x$prediction)
    yrange <- range(c(x$prediction, x$lower, x$upper),
                    finite = TRUE)
    plot(seq_len(n), x$prediction, pch = 16, col = line_col,
         ylim = yrange,
         xlab = "target index", ylab = "scalar summary",
         main = sprintf("%.0f%% intervals (n = %d)",
                        100 * (1 - alpha), n), ...)
    graphics::arrows(seq_len(n), x$lower, seq_len(n), x$upper,
                     length = 0.05, angle = 90, code = 3,
                     col = line_col)
  }
  invisible(x)
}

#' Summarise a conformal prediction-interval object
#'
#' Produces a small list of descriptive statistics about a
#' `"metahunt_conformal"` object: interval widths, quantile summaries, and
#' calibration diagnostics. Returns an object of class
#' `"summary.metahunt_conformal"` with a matching `print` method.
#'
#' @param object A `"metahunt_conformal"` object from [split_conformal()],
#'   [cross_conformal()], or [conformal_from_fit()].
#' @param ... Unused; present for S3 generic consistency.
#'
#' @return A list of class `"summary.metahunt_conformal"`. In pointwise mode
#'   (no wrapper) the list contains `n_targets`, `G_grid`, `n_cal`, `alpha`,
#'   `method`, `mean_interval_width`, `frac_finite_quantile`,
#'   `quantile_summary`, and `wrapper`. In scalar mode (wrapper supplied)
#'   the list contains `n_targets`, `n_cal`, `alpha`, `method`,
#'   `mean_interval_width`, `quantile`, `quantile_finite`, and `wrapper`.
#'
#' @examples
#' set.seed(1)
#' G <- 25; m <- 60
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m))
#' eta <- cbind(0.6 * W$w1, -0.3 * W$w1, rep(0, m))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#' W_new <- data.frame(w1 = c(0, 1))
#' res <- split_conformal(F_hat, W, W_new, K = 3,
#'                        dfspa_args = list(denoise = FALSE), seed = 1)
#' summary(res)
#'
#' @export
summary.metahunt_conformal <- function(object, ...) {
  if (!inherits(object, "metahunt_conformal")) {
    stop("`object` must be a `metahunt_conformal` object.")
  }
  if (is.null(object$wrapper)) {
    widths <- object$upper - object$lower
    finite_widths <- widths[is.finite(widths)]
    qfinite <- object$quantile[is.finite(object$quantile)]
    if (length(qfinite)) {
      qsum <- c(min = min(qfinite), median = stats::median(qfinite),
                mean = mean(qfinite), max = max(qfinite))
    } else {
      qsum <- c(min = NA_real_, median = NA_real_,
                mean = NA_real_, max = NA_real_)
    }
    out <- list(
      mode = "pointwise",
      n_targets = nrow(object$prediction),
      G_grid = ncol(object$prediction),
      n_cal = object$n_cal,
      alpha = object$alpha,
      method = object$method,
      mean_interval_width = if (length(finite_widths))
        mean(finite_widths) else NA_real_,
      frac_finite_quantile = mean(is.finite(object$quantile)),
      quantile_summary = qsum,
      wrapper = NULL
    )
  } else {
    widths <- object$upper - object$lower
    finite_widths <- widths[is.finite(widths)]
    out <- list(
      mode = "scalar",
      n_targets = length(object$prediction),
      n_cal = object$n_cal,
      alpha = object$alpha,
      method = object$method,
      mean_interval_width = if (length(finite_widths))
        mean(finite_widths) else NA_real_,
      quantile = object$quantile,
      quantile_finite = is.finite(object$quantile),
      wrapper = object$wrapper
    )
  }
  class(out) <- "summary.metahunt_conformal"
  out
}

#' @export
print.summary.metahunt_conformal <- function(x, ...) {
  cat("Summary of MetaHunt conformal prediction\n")
  cat("  method:        ", x$method, "\n")
  cat("  alpha:         ", x$alpha, "\n")
  cat("  n calibration: ", x$n_cal, "\n")
  if (identical(x$mode, "pointwise")) {
    cat("  mode:           pointwise (per grid point)\n")
    cat("  n targets:     ", x$n_targets, "\n")
    cat("  grid size:     ", x$G_grid, "\n")
    cat("  mean interval width:", signif(x$mean_interval_width, 4), "\n")
    cat("  fraction finite quantiles:",
        signif(x$frac_finite_quantile, 4), "\n")
    qs <- x$quantile_summary
    if (all(is.na(qs))) {
      cat("  quantile:       all infinite (insufficient calibration)\n")
    } else {
      cat(sprintf(
        "  quantile (finite):  min=%.4g, median=%.4g, mean=%.4g, max=%.4g\n",
        qs["min"], qs["median"], qs["mean"], qs["max"]))
    }
  } else {
    cat("  mode:           scalar (via wrapper)\n")
    cat("  n targets:     ", x$n_targets, "\n")
    if (isTRUE(x$quantile_finite)) {
      cat("  quantile:      ", signif(x$quantile, 4), "\n")
      cat("  mean interval width:", signif(x$mean_interval_width, 4), "\n")
    } else {
      cat("  quantile:       Inf (insufficient calibration)\n")
    }
  }
  invisible(x)
}

#' Empirical coverage of a conformal prediction-interval object
#'
#' Computes empirical coverage indicators of a fitted
#' `"metahunt_conformal"` object against held-out observed study-level
#' functions. The held-out studies must correspond positionally to the
#' targets used to build `object` (i.e. `F_obs[i, ]` is the observed
#' function for the same target whose prediction is in
#' `object$prediction[i, ]` or `object$prediction[i]`).
#'
#' @details
#' In **pointwise** mode each entry `(i, g)` of `F_obs` is compared against
#' `[object$lower[i, g], object$upper[i, g]]`. In **scalar** mode `F_obs`
#' is first reduced to a length-`n_target` vector via
#' [apply_wrapper()] using `object$wrapper` and the supplied
#' `grid_weights`, and then compared against
#' `[object$lower, object$upper]`.
#'
#' Coverage is a finite-sample diagnostic. Nominal coverage is
#' `1 - object$alpha`; empirical coverage will fluctuate around this
#' value due to sampling.
#'
#' @param object A `"metahunt_conformal"` object from [split_conformal()],
#'   [cross_conformal()], or [conformal_from_fit()].
#' @param F_obs An `n_target`-by-`G_grid` numeric matrix of observed
#'   study-level functions for the target studies, in the same row order
#'   as `object$prediction`.
#' @param grid_weights Optional length-`G_grid` non-negative numeric
#'   vector. Only used in scalar mode (passed to [apply_wrapper()]).
#'
#' @return A list. In **pointwise** mode the list contains
#'   \describe{
#'     \item{`pointwise`}{`n_target`-by-`G_grid` logical matrix of
#'       coverage indicators.}
#'     \item{`per_target`}{Length-`n_target` numeric vector of mean
#'       coverage across the grid for each target.}
#'     \item{`per_grid_point`}{Length-`G_grid` numeric vector of mean
#'       coverage across targets at each grid point.}
#'     \item{`overall`}{Scalar mean coverage across all entries.}
#'     \item{`nominal`}{Nominal coverage `1 - object$alpha`.}
#'   }
#'   In **scalar** mode the list contains
#'   \describe{
#'     \item{`pointwise`}{Length-`n_target` logical vector of coverage
#'       indicators.}
#'     \item{`overall`}{Scalar mean coverage.}
#'     \item{`nominal`}{Nominal coverage `1 - object$alpha`.}
#'   }
#'
#' @examples
#' set.seed(1)
#' G <- 25; m <- 80
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m))
#' eta <- cbind(0.6 * W$w1, -0.3 * W$w1, rep(0, m))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#'
#' # held-out test set: same data-generating process, same W
#' test_idx <- 1:10
#' train_idx <- setdiff(seq_len(m), test_idx)
#' res <- split_conformal(F_hat[train_idx, ], W[train_idx, , drop = FALSE],
#'                        W[test_idx, , drop = FALSE], K = 3,
#'                        dfspa_args = list(denoise = FALSE), seed = 1)
#' cov <- coverage(res, F_obs = F_hat[test_idx, , drop = FALSE])
#' cov$overall
#'
#' @export
coverage <- function(object, F_obs, grid_weights = NULL) {
  if (!inherits(object, "metahunt_conformal")) {
    stop("`object` must be a `metahunt_conformal` object (from `split_conformal()`, `cross_conformal()`, or `conformal_from_fit()`).")
  }
  if (!is.matrix(F_obs) || !is.numeric(F_obs)) {
    stop("`F_obs` must be a numeric matrix with target studies as rows.")
  }

  if (is.null(object$wrapper)) {
    n_target <- nrow(object$prediction)
    G_grid   <- ncol(object$prediction)
    if (nrow(F_obs) != n_target) {
      stop(sprintf(
        "`F_obs` must have %d rows (one per target in `object$prediction`); received %d.",
        n_target, nrow(F_obs)))
    }
    if (ncol(F_obs) != G_grid) {
      stop(sprintf(
        "`F_obs` must have %d columns (matching `ncol(object$prediction)`); received %d.",
        G_grid, ncol(F_obs)))
    }
    pw <- (F_obs >= object$lower) & (F_obs <= object$upper)
    list(
      pointwise = pw,
      per_target = rowMeans(pw),
      per_grid_point = colMeans(pw),
      overall = mean(pw),
      nominal = 1 - object$alpha
    )
  } else {
    n_target <- length(object$prediction)
    if (nrow(F_obs) != n_target) {
      stop(sprintf(
        "`F_obs` must have %d rows (one per target in `object$prediction`); received %d.",
        n_target, nrow(F_obs)))
    }
    s_obs <- apply_wrapper(F_obs, object$wrapper, grid_weights)
    pw <- (s_obs >= object$lower) & (s_obs <= object$upper)
    list(
      pointwise = pw,
      overall = mean(pw),
      nominal = 1 - object$alpha
    )
  }
}

#' @export
print.metahunt_conformal <- function(x, ...) {
  cat("MetaHunt conformal prediction\n")
  cat("  method:       ", x$method, "\n")
  cat("  alpha:        ", x$alpha, "\n")
  cat("  n calibration:", x$n_cal, "\n")
  if (is.null(x$wrapper)) {
    cat("  mode:          pointwise (per grid point)\n")
    cat("  n targets:    ", nrow(x$prediction), "\n")
    cat("  grid size:    ", ncol(x$prediction), "\n")
    qfinite <- x$quantile[is.finite(x$quantile)]
    if (length(qfinite)) {
      cat("  quantile range:", signif(min(qfinite), 4),
          "to", signif(max(qfinite), 4), "\n")
    } else {
      cat("  quantile:      Inf (insufficient calibration)\n")
    }
  } else {
    cat("  mode:          scalar (via wrapper)\n")
    cat("  n targets:    ", length(x$prediction), "\n")
    cat("  quantile:     ", signif(x$quantile, 4), "\n")
  }
  invisible(x)
}
