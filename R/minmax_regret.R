#' Minimax-regret aggregator for multisite function-valued estimands
#'
#' Implements the minimax-regret estimator of Zhang, Huang, and Imai
#' (\href{https://arxiv.org/abs/2412.11136}{arXiv:2412.11136}) for
#' aggregating site-level function estimates. Given study-level functions
#' \eqn{\hat f^{(i)}}, the estimator is
#' \deqn{\hat q = \arg\min_{q \in \Delta_{m-1}}\
#'   q^\top \hat\Gamma\, q - \hat d^\top q,
#'   \qquad
#'   \hat\Gamma_{ij} = \sum_g w_g\, \hat f^{(i)}(x_g)\hat f^{(j)}(x_g),\ \
#'   \hat d_i = \sum_g w_g\, (\hat f^{(i)}(x_g))^2,}
#' yielding the predicted target function
#' \eqn{\tilde f(x) = \sum_{i=1}^m \hat q_i\, \hat f^{(i)}(x)}.
#' Unlike [metahunt()], this method does not use study-level covariates;
#' the target is the worst-case-regret aggregator over the convex hull of
#' source functions.
#'
#' @details
#' The simplex-constrained QP is solved with [quadprog::solve.QP()].
#' A small ridge is added to \eqn{\hat\Gamma} for numerical stability when
#' source functions are highly collinear. The resulting `q` is clipped to
#' be non-negative and renormalised to sum to 1 to absorb floating-point
#' drift.
#'
#' @param F_hat An `m`-by-`G_grid` numeric matrix of source-site function
#'   evaluations on a shared grid; row `i` is \eqn{\hat f^{(i)}}.
#' @param grid_weights Optional length-`G_grid` non-negative numeric vector
#'   defining the target measure used to compute \eqn{\hat\Gamma} and
#'   \eqn{\hat d}. Defaults to uniform weights `1 / G_grid`.
#' @param ridge Non-negative scalar; replaces \eqn{\hat\Gamma} with
#'   \eqn{\hat\Gamma + \text{ridge}\cdot I} for numerical stability.
#'   Defaults to `1e-10`.
#' @param wrapper Optional reduction function (see [apply_wrapper()]). If
#'   `NULL`, `prediction` is the length-`G_grid` predicted target function;
#'   if a function, `prediction` is the scalar `wrapper(prediction)`
#'   (e.g. an ATE when `wrapper = mean` and rows of `F_hat` are CATE
#'   evaluations).
#'
#' @return An object of class `"minmax_regret"`: a list with
#'   \describe{
#'     \item{`prediction`}{Predicted target. Length-`G_grid` vector when
#'       `wrapper = NULL`; scalar otherwise.}
#'     \item{`q`}{Length-`m` simplex weights from the minimax-regret QP.}
#'     \item{`Gamma`}{The `m`-by-`m` Gram matrix used (post-ridge).}
#'     \item{`d`}{Length-`m` linear coefficient vector.}
#'     \item{`grid_weights`}{Grid weights used.}
#'     \item{`ridge`}{Ridge value used.}
#'     \item{`wrapper`}{Wrapper function or `NULL`.}
#'   }
#'
#' @references
#' Zhang, Y., Huang, M., and Imai, K. (2024). Minimax regret estimation
#' for generalizing heterogeneous treatment effects with multisite data.
#' \href{https://arxiv.org/abs/2412.11136}{arXiv:2412.11136}.
#'
#' @examples
#' set.seed(1)
#' G <- 30; m <- 6
#' x <- seq(0, 1, length.out = G)
#' F_hat <- rbind(
#'   sin(pi * x),
#'   cos(pi * x),
#'   x,
#'   0.5 * sin(pi * x) + 0.5 * x,
#'   0.3 * cos(pi * x) + 0.7 * x,
#'   0.4 * sin(pi * x) + 0.4 * cos(pi * x) + 0.2 * x
#' )
#'
#' fit <- minmax_regret(F_hat)
#' fit$q                                        # simplex weights over sources
#' length(fit$prediction)                       # G: predicted target function
#'
#' # ATE-style scalar via wrapper
#' minmax_regret(F_hat, wrapper = mean)$prediction
#'
#' @export
minmax_regret <- function(F_hat, grid_weights = NULL,
                          ridge = 1e-10, wrapper = NULL) {
  if (!is.matrix(F_hat) || !is.numeric(F_hat)) {
    stop("`F_hat` must be a numeric matrix with sources as rows.")
  }
  if (anyNA(F_hat)) stop("`F_hat` must not contain NA or NaN values.")
  m <- nrow(F_hat); G_grid <- ncol(F_hat)
  if (m < 1L || G_grid < 1L) stop("`F_hat` must be non-empty.")
  if (!is.numeric(ridge) || length(ridge) != 1L ||
      ridge < 0 || !is.finite(ridge)) {
    stop("`ridge` must be a non-negative finite scalar.")
  }

  if (is.null(grid_weights)) {
    grid_weights <- rep(1 / G_grid, G_grid)
  } else {
    if (!is.numeric(grid_weights) || length(grid_weights) != G_grid ||
        any(!is.finite(grid_weights)) || any(grid_weights < 0)) {
      stop("`grid_weights` must be a length-`ncol(F_hat)` non-negative numeric vector.")
    }
    if (sum(grid_weights) <= 0) {
      stop("`grid_weights` must have strictly positive total mass.")
    }
  }

  # Gram matrix: Gamma[i,j] = sum_g w_g F[i,g] F[j,g]
  Fw    <- sweep(F_hat, 2, grid_weights, `*`)        # m x G
  Gamma <- tcrossprod(Fw, F_hat)                     # m x m
  Gamma <- (Gamma + t(Gamma)) / 2                    # symmetrise
  d_vec <- rowSums(F_hat * Fw)                       # m

  # Single-source short circuit
  if (m == 1L) {
    q <- 1
    pred <- as.numeric(F_hat[1, ])
  } else {
    Dmat <- 2 * (Gamma + ridge * diag(m))
    Amat <- cbind(rep(1, m), diag(m))
    bvec <- c(1, rep(0, m))
    sol  <- quadprog::solve.QP(Dmat = Dmat, dvec = d_vec,
                               Amat = Amat, bvec = bvec, meq = 1L)
    q <- pmax(sol$solution, 0)
    if (sum(q) <= 0) stop("Solver returned an all-zero weight vector.")
    q <- q / sum(q)
    pred <- as.numeric(crossprod(q, F_hat))
  }

  if (!is.null(wrapper)) {
    if (!is.function(wrapper)) stop("`wrapper` must be a function or NULL.")
    pred <- apply_wrapper(matrix(pred, nrow = 1L),
                          wrapper = wrapper,
                          grid_weights = grid_weights)
    pred <- as.numeric(pred)
  }

  structure(
    list(
      prediction   = pred,
      q            = q,
      Gamma        = Gamma + ridge * diag(m),
      d            = d_vec,
      grid_weights = grid_weights,
      ridge        = ridge,
      wrapper      = wrapper
    ),
    class = "minmax_regret"
  )
}

#' @export
print.minmax_regret <- function(x, ...) {
  cat("Minimax-regret aggregator (Zhang, Huang & Imai 2024)\n")
  cat("  m (sources):  ", length(x$q), "\n")
  cat("  G (grid size):", length(x$grid_weights), "\n")
  cat("  ridge:        ", signif(x$ridge, 4), "\n")
  if (is.null(x$wrapper)) {
    cat("  prediction:    function on grid (length", length(x$prediction), ")\n")
  } else {
    cat("  prediction:    scalar =", signif(x$prediction, 4), "\n")
  }
  ord <- order(x$q, decreasing = TRUE)
  top <- min(5L, length(x$q))
  cat("  top weights:\n")
  for (k in seq_len(top)) {
    cat(sprintf("    source %d: q = %.4f\n", ord[k], x$q[ord[k]]))
  }
  invisible(x)
}
