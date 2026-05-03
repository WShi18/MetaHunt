#' @keywords internal
#' @noRd
.row_l2 <- function(F_mat, grid_weights = NULL) {
  if (is.null(grid_weights)) {
    grid_weights <- rep(1 / ncol(F_mat), ncol(F_mat))
  }
  sqrt(as.numeric((F_mat^2) %*% grid_weights))
}

#' Reconstruction-error curve for basis-rank selection
#'
#' For each candidate number of bases `K`, run [dfspa()] followed by
#' [project_to_simplex()] and report the average projection residual
#' \deqn{\mathcal E(K) = \frac{1}{m}\sum_{i=1}^m
#'        \left\|\hat f^{(i)} - \sum_{k=1}^K \hat\pi_{ik}\hat g_k\right\|_{L^2(\mu)}.}
#' Plotting `error` against `K` typically shows an elbow.
#'
#' @details
#' This is the unsupervised rank-selection criterion of Section 3.2 of the
#' paper (Equation for \eqn{\mathcal E(K)}). It does not require
#' study-level covariates.
#'
#' @param F_hat An `m`-by-`G_grid` numeric matrix of study-level function
#'   evaluations on the shared grid.
#' @param K_range Integer vector of candidate `K` values. Defaults to
#'   `2:min(nrow(F_hat) - 1, 10)`.
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector
#'   used for the \eqn{L^2(\mu)} norm; defaults to uniform.
#' @param dfspa_args Named list of extra arguments passed to [dfspa()],
#'   e.g. `list(denoise = FALSE, N = 2)`.
#'
#' @return A data frame with columns `K` (integer) and `error` (numeric).
#'   Rows where `dfspa()` or the projection fails are reported with
#'   `error = NA` and a single warning summarising the failures.
#'
#' @examples
#' set.seed(1)
#' G <- 40
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' m <- 50
#' pi_mat <- matrix(stats::rgamma(m * 3, shape = 0.5), m, 3)
#' pi_mat <- pi_mat / rowSums(pi_mat)
#' F_hat  <- pi_mat %*% basis + matrix(stats::rnorm(m * G, sd = 0.02), m, G)
#'
#' elbow <- reconstruction_error_curve(F_hat, K_range = 2:6)
#' elbow
#'
#' @export
reconstruction_error_curve <- function(F_hat, K_range = NULL,
                                       grid_weights = NULL,
                                       dfspa_args = list()) {
  if (!is.matrix(F_hat) || !is.numeric(F_hat)) {
    stop("`F_hat` must be a numeric matrix with studies as rows.")
  }
  m <- nrow(F_hat)
  if (is.null(K_range)) {
    K_range <- 2L:min(m - 1L, 10L)
  }
  K_range <- as.integer(K_range)
  if (any(K_range < 1L) || any(K_range > m)) {
    stop("`K_range` values must lie in 1:nrow(F_hat).")
  }
  if (!is.list(dfspa_args)) stop("`dfspa_args` must be a named list.")

  errs <- rep(NA_real_, length(K_range))
  failures <- character(0)
  for (j in seq_along(K_range)) {
    K <- K_range[j]
    fit <- tryCatch(
      do.call(dfspa, c(list(F_hat = F_hat, K = K,
                            grid_weights = grid_weights), dfspa_args)),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      failures <- c(failures, sprintf("K = %d: dfspa: %s", K, conditionMessage(fit)))
      next
    }
    pi_hat <- tryCatch(
      project_to_simplex(F_hat, fit$bases, grid_weights = grid_weights),
      error = function(e) e
    )
    if (inherits(pi_hat, "error")) {
      failures <- c(failures, sprintf("K = %d: projection: %s", K, conditionMessage(pi_hat)))
      next
    }
    resid <- F_hat - pi_hat %*% fit$bases
    errs[j] <- mean(.row_l2(resid, grid_weights))
  }

  if (length(failures)) {
    warning(sprintf("%d of %d K values failed:\n%s",
                    length(failures), length(K_range),
                    paste0("  ", failures, collapse = "\n")))
  }

  data.frame(K = K_range, error = errs)
}

#' Cross-validated prediction-error curve for basis-rank selection
#'
#' For each candidate `K`, perform k-fold cross-validation at the study
#' level. Within each fold, the full MetaHunt pipeline (d-fSPA + constrained
#' projection + weight model) is refit on the training studies, and the
#' held-out studies' functions are predicted from their study-level
#' covariates. The prediction error is
#' \eqn{\|\hat f^{(i)} - \tilde f^{(i)}\|_{L^2(\mu)}} averaged over
#' held-out studies and then over folds.
#'
#' @details
#' This is the supervised rank-selection criterion of Section 3.2 of the
#' paper. Each held-out study is excluded from both basis hunting and
#' weight-model fitting.
#'
#' @param F_hat An `m`-by-`G_grid` numeric matrix of study-level function
#'   evaluations.
#' @param W An `m`-by-`p` matrix or data frame of study-level covariates.
#' @param K_range Integer vector of candidate `K` values. Defaults to
#'   `2:min(nrow(F_hat) - 1, 10)`.
#' @param n_folds Integer number of CV folds (default `5`).
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector.
#' @param dfspa_args Named list of extra arguments for [dfspa()].
#' @param weight_model_args Named list of extra arguments for
#'   [fit_weight_model()].
#' @param seed Optional integer seed for reproducible fold assignment; if
#'   `NULL` no seeding is performed.
#'
#' @return A data frame with columns `K`, `cv_error` (mean over folds),
#'   and `cv_se` (standard error across folds). The per-fold error matrix
#'   is attached as the attribute `"fold_errors"`
#'   (`length(K_range)`-by-`n_folds`). Folds where the pipeline fails
#'   contribute `NA` and are summarised in a single warning.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' G <- 40; m <- 80; K_true <- 3
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- as.matrix(W) %*% cbind(c(1, -0.8), c(-0.5, 1.2), c(0, 0))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.02), m, G)
#' cv <- cv_error_curve(F_hat, W, K_range = 2:5, n_folds = 4, seed = 1)
#' cv
#' }
#'
#' @export
cv_error_curve <- function(F_hat, W, K_range = NULL, n_folds = 5L,
                           grid_weights = NULL,
                           dfspa_args = list(),
                           weight_model_args = list(),
                           seed = NULL) {
  if (!is.matrix(F_hat) || !is.numeric(F_hat)) {
    stop("`F_hat` must be a numeric matrix with studies as rows.")
  }
  m <- nrow(F_hat)
  if (is.matrix(W)) W <- as.data.frame(W)
  if (!is.data.frame(W)) stop("`W` must be a matrix or data frame.")
  if (nrow(W) != m) stop("`W` must have the same number of rows as `F_hat`.")
  if (!is.list(dfspa_args) || !is.list(weight_model_args)) {
    stop("`dfspa_args` and `weight_model_args` must be named lists.")
  }
  if (!is.numeric(n_folds) || length(n_folds) != 1L ||
      n_folds < 2L || n_folds != as.integer(n_folds)) {
    stop("`n_folds` must be an integer >= 2.")
  }
  n_folds <- as.integer(n_folds)
  if (n_folds > m) stop("`n_folds` cannot exceed the number of studies.")

  if (is.null(K_range)) {
    K_range <- 2L:min(m - 1L, 10L)
  }
  K_range <- as.integer(K_range)

  if (!is.null(seed)) withr::local_seed(seed)
  fold_id <- sample(rep(seq_len(n_folds), length.out = m))

  fold_errors <- matrix(NA_real_, nrow = length(K_range), ncol = n_folds,
                        dimnames = list(NULL, paste0("fold", seq_len(n_folds))))
  failures <- character(0)

  for (j in seq_along(K_range)) {
    K <- K_range[j]
    for (f in seq_len(n_folds)) {
      train <- which(fold_id != f)
      test  <- which(fold_id == f)
      result <- tryCatch({
        fit_f <- do.call(dfspa,
                         c(list(F_hat = F_hat[train, , drop = FALSE], K = K,
                                grid_weights = grid_weights), dfspa_args))
        pi_tr <- project_to_simplex(F_hat[train, , drop = FALSE],
                                    fit_f$bases, grid_weights = grid_weights)
        wm <- do.call(fit_weight_model,
                      c(list(pi_hat = pi_tr, W = W[train, , drop = FALSE]),
                        weight_model_args))
        f_pred <- predict_target(fit_f, wm, W[test, , drop = FALSE])
        mean(.row_l2(F_hat[test, , drop = FALSE] - f_pred, grid_weights))
      }, error = function(e) e)
      if (inherits(result, "error")) {
        failures <- c(failures, sprintf("K = %d, fold %d: %s",
                                        K, f, conditionMessage(result)))
      } else {
        fold_errors[j, f] <- result
      }
    }
  }

  if (length(failures)) {
    warning(sprintf("%d (K, fold) combinations failed:\n%s",
                    length(failures),
                    paste0("  ", failures, collapse = "\n")))
  }

  n_ok     <- rowSums(!is.na(fold_errors))
  cv_error <- ifelse(n_ok > 0L, rowMeans(fold_errors, na.rm = TRUE), NA_real_)
  cv_sd    <- ifelse(n_ok > 1L, apply(fold_errors, 1, stats::sd, na.rm = TRUE),
                     NA_real_)
  cv_se    <- ifelse(n_ok > 1L, cv_sd / sqrt(n_ok), NA_real_)

  if (any(n_ok == 0L)) {
    warning(sprintf(
      "K = %s had no successful folds; cv_error and cv_se reported as NA.",
      paste(K_range[n_ok == 0L], collapse = ", ")
    ), call. = FALSE)
  }

  out <- data.frame(K = K_range, cv_error = cv_error, cv_se = cv_se)
  attr(out, "fold_errors") <- fold_errors
  out
}
