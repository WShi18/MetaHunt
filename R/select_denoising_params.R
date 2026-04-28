#' Choose d-fSPA denoising parameters by cross-validation
#'
#' At a fixed `K`, performs k-fold cross-validation over a grid of denoising
#' parameter pairs `(N, Delta)` for [dfspa()]. For each candidate pair and
#' each fold, the full MetaHunt pipeline is fit on the out-of-fold studies
#' and predicts the held-out studies' functions. The pair with the lowest
#' average prediction error is selected.
#'
#' @details
#' This is the cross-validated tuning of the denoising parameters discussed
#' in Section 3.1 of the paper. Joint tuning over `(K, N, Delta)` is not
#' supported because it scales poorly; if you also want to choose `K`, do
#' it first via [cv_error_curve()] and then call this function at the
#' selected `K`.
#'
#' Default candidate grids:
#' * `N_grid    = m * c(NA, NA, NA)` resolved at runtime to
#'   `c(0.2, 0.5, 1.0) * log(m)`.
#' * `Delta_grid = max_pairwise_dist * c(0.05, 0.10, 0.20, 0.30)`.
#'
#' Pass your own `N_grid` / `Delta_grid` (in original units) to override.
#'
#' @param F_hat An `m`-by-`G_grid` numeric matrix of study-level function
#'   evaluations.
#' @param W An `m`-by-`p` matrix or data frame of study-level covariates.
#' @param K Integer number of basis functions (fixed during this search).
#' @param N_grid Optional numeric vector of candidate `N` values. Defaults
#'   to `c(0.2, 0.5, 1.0) * log(nrow(F_hat))`.
#' @param Delta_grid Optional numeric vector of candidate `Delta` values.
#'   Defaults to `c(0.05, 0.10, 0.20, 0.30)` times the maximum pairwise
#'   \eqn{L^2(\mu)} distance among studies.
#' @param n_folds Integer number of folds (default `5`).
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector.
#' @param dfspa_args Named list of additional arguments for [dfspa()]
#'   (e.g. `list()`; do not include `N`, `Delta`, or `denoise` here, they
#'   are set by the search).
#' @param weight_model_args Named list of additional arguments for
#'   [fit_weight_model()].
#' @param seed Optional integer seed for reproducible fold assignment.
#'
#' @return An object of class `metahunt_denoising_search`: a list with
#'   \describe{
#'     \item{`grid`}{A data frame with one row per `(N, Delta)` pair,
#'       columns `N`, `Delta`, `cv_error`, `cv_se`, `n_folds_ok`.}
#'     \item{`best`}{A list with the `(N, Delta)` minimising `cv_error`.}
#'     \item{`K`, `n_folds`, `grid_weights`}{Inputs echoed back for
#'       traceability.}
#'   }
#'
#' @seealso [cv_error_curve()] for selecting `K`, [dfspa()] for the
#'   underlying basis-hunting algorithm.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' G <- 30; m <- 80
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- as.matrix(W) %*% cbind(c(1, -0.5), c(-0.4, 1), c(0, 0))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' F_hat <- pi_true %*% basis + matrix(rnorm(m * G, sd = 0.05), m, G)
#'
#' tune <- select_denoising_params(F_hat, W, K = 3, n_folds = 4, seed = 1)
#' tune$grid
#' tune$best
#' }
#'
#' @export
select_denoising_params <- function(F_hat, W, K,
                                    N_grid = NULL,
                                    Delta_grid = NULL,
                                    n_folds = 5L,
                                    grid_weights = NULL,
                                    dfspa_args = list(),
                                    weight_model_args = list(),
                                    seed = NULL) {
  if (!is.matrix(F_hat) || !is.numeric(F_hat)) {
    stop("`F_hat` must be a numeric matrix with studies as rows.")
  }
  m <- nrow(F_hat); G_grid <- ncol(F_hat)
  if (is.matrix(W)) W <- as.data.frame(W)
  if (!is.data.frame(W)) stop("`W` must be a matrix or data frame.")
  if (nrow(W) != m) stop("`W` must have the same number of rows as `F_hat`.")
  if (!is.numeric(K) || length(K) != 1L || K < 1L || K != as.integer(K)) {
    stop("`K` must be a positive integer.")
  }
  K <- as.integer(K)
  if (!is.list(dfspa_args) || !is.list(weight_model_args)) {
    stop("`dfspa_args` and `weight_model_args` must be named lists.")
  }
  reserved <- c("N", "Delta", "denoise")
  bad <- intersect(names(dfspa_args), reserved)
  if (length(bad)) {
    stop("`dfspa_args` must not contain ", paste(bad, collapse = ", "),
         "; these are set by the search.")
  }
  if (!is.numeric(n_folds) || length(n_folds) != 1L ||
      n_folds < 2L || n_folds != as.integer(n_folds)) {
    stop("`n_folds` must be an integer >= 2.")
  }
  n_folds <- as.integer(n_folds)
  if (n_folds > m) stop("`n_folds` cannot exceed the number of studies.")

  if (is.null(N_grid)) {
    N_grid <- c(0.2, 0.5, 1.0) * log(m)
  }
  if (!is.numeric(N_grid) || any(!is.finite(N_grid)) || any(N_grid <= 0)) {
    stop("`N_grid` must be a vector of positive finite numbers.")
  }

  if (is.null(Delta_grid)) {
    if (is.null(grid_weights)) {
      gw <- rep(1 / G_grid, G_grid)
    } else {
      gw <- grid_weights
    }
    sqw <- sqrt(gw)
    F_tilde <- sweep(F_hat, 2, sqw, `*`)
    max_dist <- max(stats::dist(F_tilde))
    Delta_grid <- max_dist * c(0.05, 0.10, 0.20, 0.30)
  }
  if (!is.numeric(Delta_grid) || any(!is.finite(Delta_grid)) ||
      any(Delta_grid < 0)) {
    stop("`Delta_grid` must be a vector of non-negative finite numbers.")
  }

  combos <- expand.grid(N = N_grid, Delta = Delta_grid,
                        KEEP.OUT.ATTRS = FALSE)
  if (!is.null(seed)) withr::local_seed(seed)
  fold_id <- sample(rep(seq_len(n_folds), length.out = m))

  fold_errors <- matrix(NA_real_, nrow = nrow(combos), ncol = n_folds)
  failures <- character(0)

  for (j in seq_len(nrow(combos))) {
    Nj <- combos$N[j]; Dj <- combos$Delta[j]
    for (f in seq_len(n_folds)) {
      train <- which(fold_id != f)
      test  <- which(fold_id == f)
      result <- tryCatch({
        fit_f <- do.call(dfspa,
          c(list(F_hat = F_hat[train, , drop = FALSE], K = K,
                 grid_weights = grid_weights,
                 N = Nj, Delta = Dj, denoise = TRUE),
            dfspa_args))
        pi_tr <- project_to_simplex(F_hat[train, , drop = FALSE],
                                    fit_f$bases,
                                    grid_weights = grid_weights)
        wm <- do.call(fit_weight_model,
                      c(list(pi_hat = pi_tr,
                             W = W[train, , drop = FALSE]),
                        weight_model_args))
        f_pred <- predict_target(fit_f, wm, W[test, , drop = FALSE])
        mean(.row_l2(F_hat[test, , drop = FALSE] - f_pred, grid_weights))
      }, error = function(e) e)
      if (inherits(result, "error")) {
        failures <- c(failures,
                      sprintf("(N=%g, Delta=%g, fold %d): %s",
                              Nj, Dj, f, conditionMessage(result)))
      } else {
        fold_errors[j, f] <- result
      }
    }
  }
  if (length(failures)) {
    warning(sprintf("%d (N, Delta, fold) combinations failed:\n%s",
                    length(failures),
                    paste0("  ", failures, collapse = "\n")))
  }

  cv_error <- rowMeans(fold_errors, na.rm = TRUE)
  cv_sd    <- apply(fold_errors, 1, stats::sd, na.rm = TRUE)
  n_ok     <- rowSums(!is.na(fold_errors))
  cv_se    <- ifelse(n_ok > 1L, cv_sd / sqrt(n_ok), NA_real_)

  grid <- data.frame(N = combos$N, Delta = combos$Delta,
                     cv_error = cv_error, cv_se = cv_se,
                     n_folds_ok = n_ok)

  if (all(is.na(grid$cv_error))) {
    stop("All (N, Delta) combinations failed for every fold.")
  }
  best_idx <- which.min(grid$cv_error)
  best <- list(N = grid$N[best_idx], Delta = grid$Delta[best_idx],
               cv_error = grid$cv_error[best_idx])

  structure(
    list(grid = grid, best = best,
         K = K, n_folds = n_folds, grid_weights = grid_weights),
    class = "metahunt_denoising_search"
  )
}

#' Print method for d-fSPA denoising parameter search results
#'
#' @param x An object of class `metahunt_denoising_search` returned by
#'   [select_denoising_params()].
#' @param ... Unused; present for S3 generic compatibility.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.metahunt_denoising_search <- function(x, ...) {
  cat("Denoising parameter search (CV at fixed K)\n")
  cat("  K:           ", x$K, "\n")
  cat("  n_folds:     ", x$n_folds, "\n")
  n_N     <- length(unique(x$grid$N))
  n_Delta <- length(unique(x$grid$Delta))
  cat(sprintf("  grid size:    %d (%d N values x %d Delta values)\n",
              nrow(x$grid), n_N, n_Delta))
  cat(sprintf("  best:         N = %s, Delta = %s, cv_error = %s\n",
              signif(x$best$N, 4),
              signif(x$best$Delta, 4),
              signif(x$best$cv_error, 4)))
  cat("\n")
  top_n <- min(5L, nrow(x$grid))
  cat(sprintf("  Top %d (lowest cv_error):\n", top_n))
  ord <- order(x$grid$cv_error, na.last = TRUE)
  top_grid <- x$grid[ord[seq_len(top_n)], , drop = FALSE]
  for (col in c("N", "Delta", "cv_error", "cv_se")) {
    top_grid[[col]] <- signif(top_grid[[col]], 4)
  }
  print(top_grid, row.names = FALSE)
  invisible(x)
}
