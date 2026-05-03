#' Project study-level functions onto the simplex spanned by basis functions
#'
#' For each study `i`, solves the constrained projection
#' \deqn{\hat\pi_i = \arg\min_{\pi \in \Delta_{K-1}}
#'   \left\| \hat f^{(i)} - \sum_{k=1}^K \pi_k \hat g_k \right\|_{L^2(\mu)}}
#' where the norm is the weighted \eqn{L^2} norm defined by `grid_weights`.
#' This is Equation (3) of the paper and yields the study-specific weights
#' `hat pi_i` used downstream for weight-model fitting and prediction.
#'
#' @details
#' The projection reduces to the quadratic program
#' \deqn{\min_{\pi \in \mathbb R^K}\ \pi^\top D\,\pi - 2\,d^\top \pi
#'   \quad \text{s.t.} \quad \mathbf 1^\top \pi = 1,\ \pi \ge 0}
#' with \eqn{D = G W G^\top} and \eqn{d = G W f^{(i)}}, where
#' \eqn{G} is the `K`-by-`G_grid` basis matrix, \eqn{W = \mathrm{diag}(}`grid_weights`\eqn{)},
#' and \eqn{f^{(i)}} is the `i`-th row of `F_hat`. Solved via
#' [quadprog::solve.QP()]. A tiny ridge is added to `D` for numerical stability.
#'
#' @param F_hat An `m`-by-`G_grid` numeric matrix; row `i` is the study
#'   function \eqn{\hat f^{(i)}} evaluated on the shared grid.
#' @param bases A `K`-by-`G_grid` numeric matrix of basis functions on the
#'   same grid, typically the `bases` slot of a [dfspa()] result.
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector of
#'   grid weights defining the \eqn{L^2(\mu)} inner product. Defaults to
#'   uniform weights `1 / G_grid`.
#' @param ridge Small non-negative scalar added to the diagonal of the QP
#'   Hessian for numerical stability. Defaults to `1e-10`.
#'
#' @return An `m`-by-`K` numeric matrix of simplex weights; rows sum to 1
#'   and entries are non-negative.
#'
#' @examples
#' set.seed(1)
#' G <- 40
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)
#' true_pi <- rbind(diag(3), c(0.4, 0.3, 0.3), c(0.1, 0.7, 0.2))
#' F_hat <- true_pi %*% basis
#' fit <- dfspa(F_hat, K = 3, denoise = FALSE)
#' pi_hat <- project_to_simplex(F_hat, fit$bases)
#' round(pi_hat, 3)
#'
#' @export
project_to_simplex <- function(F_hat, bases, grid_weights = NULL,
                               ridge = 1e-10) {
  if (!is.matrix(F_hat) || !is.numeric(F_hat)) {
    stop("`F_hat` must be a numeric matrix with studies as rows.")
  }
  if (!is.matrix(bases) || !is.numeric(bases)) {
    stop("`bases` must be a numeric matrix with bases as rows.")
  }
  if (anyNA(F_hat) || anyNA(bases)) {
    stop("`F_hat` and `bases` must not contain NA values.")
  }
  if (ncol(F_hat) != ncol(bases)) {
    stop("`F_hat` and `bases` must have the same number of columns (grid size).")
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || ridge < 0 || !is.finite(ridge)) {
    stop("`ridge` must be a non-negative finite scalar.")
  }

  m <- nrow(F_hat)
  G_grid <- ncol(F_hat)
  K <- nrow(bases)

  if (is.null(grid_weights)) {
    grid_weights <- rep(1 / G_grid, G_grid)
  } else {
    if (!is.numeric(grid_weights) || length(grid_weights) != G_grid ||
        any(!is.finite(grid_weights)) || any(grid_weights < 0)) {
      stop("`grid_weights` must be a length-`ncol(F_hat)` non-negative numeric vector.")
    }
  }

  if (K == 1L) {
    return(matrix(1, nrow = m, ncol = 1L))
  }

  Gw   <- sweep(bases, 2, grid_weights, `*`)   # K x G_grid
  Dmat <- tcrossprod(Gw, bases)                # K x K, = G W G^T
  Dmat <- Dmat + ridge * diag(K)

  Amat <- cbind(rep(1, K), diag(K))            # K x (K + 1)
  bvec <- c(1, rep(0, K))

  pi_hat <- matrix(0, m, K)
  for (i in seq_len(m)) {
    dvec <- as.numeric(Gw %*% F_hat[i, ])
    sol  <- quadprog::solve.QP(Dmat = Dmat, dvec = dvec,
                               Amat = Amat, bvec = bvec, meq = 1L)
    pi_i <- pmax(sol$solution, 0)
    s <- sum(pi_i)
    if (s <= 0) {
      warning("project_to_simplex(): solver returned all-zero weights for row ",
              i, "; falling back to uniform weights.", call. = FALSE)
      pi_hat[i, ] <- 1 / K
    } else {
      pi_hat[i, ] <- pi_i / s
    }
  }
  pi_hat
}

#' Fit a weight model mapping study-level covariates to simplex weights
#'
#' Given a matrix of simplex-valued weights \eqn{\hat\pi_1,\ldots,\hat\pi_m}
#' (e.g. from [project_to_simplex()]) and associated study-level covariates
#' \eqn{\mathbf W_1,\ldots,\mathbf W_m}, fit a model
#' \eqn{\widehat{\mathcal M}:\mathbf W \mapsto \boldsymbol\pi}.
#' The default method is Dirichlet regression via the `DirichletReg` package.
#'
#' @details
#' Dirichlet regression cannot handle weights exactly at the simplex boundary
#' (`0` or `1`), which frequently arise after constrained projection. Before
#' fitting, rows of `pi_hat` are shrunk toward the barycenter via
#' \eqn{\tilde\pi = (\pi + \varepsilon) / (1 + K\varepsilon)}, with
#' \eqn{\varepsilon} set by `boundary_eps`.
#'
#' @param pi_hat An `m`-by-`K` numeric matrix of simplex weights; rows must
#'   be non-negative and sum to 1 (up to tolerance `1e-6`).
#' @param W An `m`-by-`p` matrix or data frame of study-level covariates.
#' @param method Weight-model method. Currently only `"dirichlet"` is
#'   supported.
#' @param boundary_eps Small positive scalar used to shrink weights away
#'   from the simplex boundary before Dirichlet fitting. Defaults to `1e-4`.
#' @param formula Optional RHS-only formula (e.g. `~ x1 + I(x2^2)`) describing
#'   the covariate part of the Dirichlet regression. Defaults to `~ .`
#'   (all columns of `W`).
#' @param ... Passed through to [DirichletReg::DirichReg()].
#'
#' @return An object of class `"metahunt_weight_model"`: a list with the
#'   fitted model, formula, method, `K`, and training covariate names.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- 80; K <- 3; p <- 2
#' W <- matrix(rnorm(m * p), m, p); colnames(W) <- c("w1", "w2")
#' # generate simplex weights driven by W
#' eta <- cbind(0.5 * W[, 1], -0.3 * W[, 2], rep(0, m))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' pi_hat <- pi_true + matrix(rnorm(m * K, sd = 0.01), m, K)
#' pi_hat <- pmax(pi_hat, 0); pi_hat <- pi_hat / rowSums(pi_hat)
#' model <- fit_weight_model(pi_hat, W)
#' predict(model, newdata = matrix(c(0, 0), 1, 2, dimnames = list(NULL, c("w1","w2"))))
#' }
#'
#' @export
fit_weight_model <- function(pi_hat, W, method = c("dirichlet"),
                             boundary_eps = 1e-4, formula = NULL, ...) {
  method <- match.arg(method)

  if (!is.numeric(pi_hat)) stop("`pi_hat` must be numeric.")
  pi_hat <- as.matrix(pi_hat)
  if (anyNA(pi_hat)) stop("`pi_hat` must not contain NA values.")
  if (any(pi_hat < -1e-8)) stop("`pi_hat` must have non-negative entries.")
  pi_hat <- pmax(pi_hat, 0)
  rs <- rowSums(pi_hat)
  if (any(abs(rs - 1) > 1e-6)) {
    stop("Rows of `pi_hat` must sum to 1 (tolerance 1e-6).")
  }
  m <- nrow(pi_hat); K <- ncol(pi_hat)

  if (is.matrix(W)) {
    W <- as.data.frame(W)
  }
  if (!is.data.frame(W)) {
    stop("`W` must be a matrix or data frame.")
  }
  if (nrow(W) != m) {
    stop("`W` must have the same number of rows as `pi_hat`.")
  }
  if ("..Y.." %in% names(W)) {
    stop("`W` must not contain a column named '..Y..' (reserved).")
  }

  if (!is.numeric(boundary_eps) || length(boundary_eps) != 1L ||
      !is.finite(boundary_eps) || boundary_eps < 0 || boundary_eps >= 1 / K) {
    stop("`boundary_eps` must be a scalar in [0, 1/K).")
  }

  switch(method,
    dirichlet = fit_dirichlet_(pi_hat, W, boundary_eps, formula, ...)
  )
}

fit_dirichlet_ <- function(pi_hat, W, boundary_eps, formula, ...) {
  K <- ncol(pi_hat)
  pi_adj <- (pi_hat + boundary_eps) / (1 + K * boundary_eps)

  predictor_names <- names(W)
  if (is.null(formula)) {
    rhs <- if (length(predictor_names) == 0L) "1"
           else paste(predictor_names, collapse = " + ")
  } else {
    if (!inherits(formula, "formula") || length(formula) != 2L) {
      stop("`formula` must be a one-sided formula, e.g. ~ x1 + x2.")
    }
    rhs <- paste(deparse(formula[[2]]), collapse = " ")
  }

  metahunt_Y <- DirichletReg::DR_data(pi_adj)
  fmla <- stats::as.formula(paste("metahunt_Y ~", rhs))
  environment(fmla) <- environment()

  # Inline the formula into the call so that DirichReg's internal match.call()
  # sees the formula object, not the symbol `fmla`.
  fit_call <- substitute(DirichletReg::DirichReg(FMLA, data = W),
                         list(FMLA = fmla))
  fit <- eval(fit_call)

  structure(
    list(
      method          = "dirichlet",
      fit             = fit,
      formula         = fmla,
      K               = K,
      predictor_names = predictor_names,
      boundary_eps    = boundary_eps
    ),
    class = "metahunt_weight_model"
  )
}

#' Predict simplex weights for new study-level covariates
#'
#' @param object A fitted `"metahunt_weight_model"` from [fit_weight_model()].
#' @param newdata A matrix or data frame of new study-level covariates with
#'   the same columns used at fitting.
#' @param ... Ignored.
#'
#' @return An `nrow(newdata)`-by-`K` numeric matrix of predicted simplex
#'   weights (component means); rows sum to 1.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- 40; K <- 3
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- cbind(0.5 * W$w1, -0.3 * W$w2, rep(0, m))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' pi_hat <- pi_true + matrix(rnorm(m * K, sd = 0.01), m, K)
#' pi_hat <- pmax(pi_hat, 0); pi_hat <- pi_hat / rowSums(pi_hat)
#' model <- fit_weight_model(pi_hat, W)
#' predict(model, newdata = data.frame(w1 = c(0, 1), w2 = c(0, -1)))
#' }
#'
#' @export
predict.metahunt_weight_model <- function(object, newdata, ...) {
  if (is.matrix(newdata)) newdata <- as.data.frame(newdata)
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a matrix or data frame.")
  }
  missing_cols <- setdiff(object$predictor_names, names(newdata))
  if (length(missing_cols)) {
    stop("`newdata` is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  preds <- stats::predict(object$fit, newdata = newdata)
  preds <- unname(as.matrix(preds))
  preds
}

#' @export
print.metahunt_weight_model <- function(x, ...) {
  cat("MetaHunt weight model\n")
  cat("  method:     ", x$method, "\n")
  cat("  K:          ", x$K, "\n")
  cat("  predictors: ", paste(x$predictor_names, collapse = ", "), "\n")
  invisible(x)
}

#' Extract coefficients from a MetaHunt weight model
#'
#' Returns the regression coefficients from the underlying weight-model fit.
#' For the default `"dirichlet"` method this delegates to
#' [DirichletReg::DirichReg()]'s `coef()` method.
#'
#' @param object A fitted `"metahunt_weight_model"` from [fit_weight_model()].
#' @param ... Passed through to [stats::coef()].
#'
#' @return The coefficient vector / matrix returned by
#'   `DirichletReg::DirichReg()`'s `coef()` method (numeric vector or matrix
#'   depending on the parametrisation used by the underlying fit).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' m <- 60; K <- 3
#' W <- data.frame(w1 = rnorm(m), w2 = rnorm(m))
#' eta <- cbind(0.5 * W$w1, -0.3 * W$w2, rep(0, m))
#' pi_true <- exp(eta) / rowSums(exp(eta))
#' pi_hat <- pi_true + matrix(rnorm(m * K, sd = 0.01), m, K)
#' pi_hat <- pmax(pi_hat, 0); pi_hat <- pi_hat / rowSums(pi_hat)
#' model <- fit_weight_model(pi_hat, W)
#' coef(model)
#' }
#'
#' @export
coef.metahunt_weight_model <- function(object, ...) {
  if (!inherits(object, "metahunt_weight_model")) {
    stop("`object` must be a `metahunt_weight_model` object.")
  }
  stats::coef(object$fit, ...)
}
