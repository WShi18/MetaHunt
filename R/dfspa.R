#' Denoised functional Successive Projection Algorithm (d-fSPA)
#'
#' Recovers a set of `K` latent basis functions from a collection of
#' study-level function estimates under the low-rank cross-study heterogeneity
#' assumption of Shi, Imai, and Zhang. Implements Algorithm 1 of the paper
#' ("The d-fSPA Algorithm for basis hunting").
#'
#' @details
#' Each study-level function is represented by its evaluations on a shared
#' grid of `G` points. The (weighted) \eqn{L^2(\mu)} inner product is
#' \eqn{\langle f,g\rangle = \sum_{j=1}^G w_j f(x_j) g(x_j)}, where the
#' `grid_weights` `w_j` are proportional to the measure \eqn{\mu}. If not
#' supplied, uniform weights `1 / G` are used.
#'
#' Denoising follows Jin (2024): for each study `i`, let
#' \eqn{B_\Delta(\hat f^{(i)}) = \{j : \|\hat f^{(j)} - \hat f^{(i)}\| \le
#' \Delta\}}. If `|B_Δ(f̂^{(i)})| < N`, study `i` is discarded; otherwise
#' `f̂^{(i)}` is replaced by the average of the functions in `B_Δ`.
#' After denoising, the functional SPA step iteratively selects, at each of
#' the `K` iterations, the remaining function with the largest norm after
#' projecting out the span of previously selected bases.
#'
#' Default tuning parameters follow the heuristics of the paper:
#' `N = 0.5 * log(m)` and
#' `Delta = max_{ij} ||f̂^{(i)} - f̂^{(j)}|| / 10`.
#'
#' @param F_hat An `m`-by-`G` numeric matrix where row `i` is the evaluation
#'   of the estimated function \eqn{\hat f^{(i)}} at `G` grid points.
#' @param K Integer number of basis functions to recover. Must satisfy
#'   `1 <= K <= m` after denoising.
#' @param grid_weights Optional length-`G` non-negative numeric vector of
#'   grid weights defining the \eqn{L^2(\mu)} inner product. Defaults to
#'   uniform weights `1 / G`.
#' @param N,Delta Optional numeric tuning parameters controlling denoising.
#'   See Details.
#' @param denoise Logical; if `FALSE`, the denoising step is skipped and
#'   plain fSPA is run. Defaults to `TRUE`.
#'
#' @return An object of class `"dfspa"`: a list containing
#' \describe{
#'   \item{`bases`}{A `K`-by-`G` matrix whose rows are the recovered basis
#'     functions evaluated on the grid (denoised, if applicable).}
#'   \item{`selected`}{Length-`K` integer vector of the selected row indices
#'     into the post-denoising function matrix.}
#'   \item{`original_indices`}{Length-`K` integer vector of the selected
#'     study indices in the original input `F_hat` (before any rows were
#'     dropped by denoising).}
#'   \item{`kept`}{Integer vector of row indices of `F_hat` that survived
#'     denoising.}
#'   \item{`F_denoised`}{The post-denoising function matrix
#'     (`length(kept)`-by-`G`).}
#'   \item{`grid_weights`}{Grid weights used.}
#'   \item{`N`, `Delta`}{Tuning parameters actually used (or `NA` when
#'     `denoise = FALSE`).}
#'   \item{`K`}{Number of bases requested.}
#'   \item{`call`}{The matched call.}
#' }
#'
#' @examples
#' set.seed(1)
#' G <- 50
#' x <- seq(0, 1, length.out = G)
#' basis <- rbind(sin(pi * x), cos(pi * x), x)          # 3 true bases
#' pi_mat <- rbind(diag(3),                              # 3 pure studies
#'                 c(0.5, 0.3, 0.2),
#'                 c(0.2, 0.5, 0.3),
#'                 c(0.3, 0.3, 0.4))
#' F_hat <- pi_mat %*% basis                             # m = 6, G = 50
#' fit <- dfspa(F_hat, K = 3, denoise = FALSE)
#' fit$original_indices    # should be a permutation of 1, 2, 3
#'
#' @export
dfspa <- function(F_hat, K, grid_weights = NULL,
                  N = NULL, Delta = NULL, denoise = TRUE) {
  call <- match.call()

  if (!is.matrix(F_hat) || !is.numeric(F_hat)) {
    stop("`F_hat` must be a numeric matrix with studies as rows.")
  }
  if (anyNA(F_hat)) {
    stop("`F_hat` must not contain NA or NaN values.")
  }
  m <- nrow(F_hat)
  G <- ncol(F_hat)
  if (m < 2L) stop("`F_hat` must have at least two rows (studies).")

  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) ||
      K < 1 || K != as.integer(K)) {
    stop("`K` must be a positive integer.")
  }
  K <- as.integer(K)
  if (K > m) stop("`K` must not exceed the number of studies `m`.")

  if (is.null(grid_weights)) {
    grid_weights <- rep(1 / G, G)
  } else {
    if (!is.numeric(grid_weights) || length(grid_weights) != G ||
        any(!is.finite(grid_weights)) || any(grid_weights < 0)) {
      stop("`grid_weights` must be a length-`ncol(F_hat)` non-negative numeric vector.")
    }
    if (sum(grid_weights) <= 0) {
      stop("`grid_weights` must have strictly positive total mass.")
    }
  }

  sqw <- sqrt(grid_weights)
  F_tilde <- sweep(F_hat, 2, sqw, `*`)
  D <- as.matrix(stats::dist(F_tilde))

  if (denoise) {
    if (is.null(N))     N     <- 0.5 * log(m)
    if (is.null(Delta)) Delta <- max(D) / 10
    if (!is.numeric(N) || length(N) != 1L || !is.finite(N) || N <= 0) {
      stop("`N` must be a positive finite scalar.")
    }
    if (!is.numeric(Delta) || length(Delta) != 1L || !is.finite(Delta) ||
        Delta < 0) {
      stop("`Delta` must be a non-negative finite scalar.")
    }
    nbr_counts <- rowSums(D <= Delta)           # includes self
    kept <- which(nbr_counts >= N)
    if (length(kept) < K) {
      stop(sprintf(
        "Only %d studies survive denoising but K = %d; consider increasing `Delta`, decreasing `N`, or reducing `K`.",
        length(kept), K))
    }
    F_denoised <- matrix(0, nrow = length(kept), ncol = G)
    for (ii in seq_along(kept)) {
      nbrs <- which(D[kept[ii], ] <= Delta)
      F_denoised[ii, ] <- colMeans(F_hat[nbrs, , drop = FALSE])
    }
  } else {
    N <- NA_real_
    Delta <- NA_real_
    kept <- seq_len(m)
    F_denoised <- F_hat
  }

  m_kept <- length(kept)
  F_denoised_tilde <- sweep(F_denoised, 2, sqw, `*`)

  H <- F_denoised_tilde
  selected <- integer(K)
  tol <- .Machine$double.eps^0.5 * max(1, max(rowSums(H^2)))
  for (k in seq_len(K)) {
    norms_sq <- rowSums(H^2)
    s_k <- which.max(norms_sq)
    if (norms_sq[s_k] <= tol) {
      stop(sprintf(
        "d-fSPA halted at k = %d because residual norms are numerically zero; K may be too large.",
        k))
    }
    selected[k] <- s_k
    u <- H[s_k, ]
    u <- u / sqrt(sum(u^2))
    proj_coef <- as.numeric(H %*% u)
    H <- H - outer(proj_coef, u)
  }

  structure(
    list(
      bases            = F_denoised[selected, , drop = FALSE],
      selected         = selected,
      original_indices = kept[selected],
      kept             = kept,
      F_denoised       = F_denoised,
      grid_weights     = grid_weights,
      N                = N,
      Delta            = Delta,
      K                = K,
      call             = call
    ),
    class = "dfspa"
  )
}

#' @export
print.dfspa <- function(x, ...) {
  cat("d-fSPA basis hunting result\n")
  cat("  studies kept after denoising:", length(x$kept),
      "of", nrow(x$F_denoised) + max(0L, length(x$kept) - nrow(x$F_denoised)), "\n")
  cat("  K (bases recovered):         ", x$K, "\n")
  cat("  selected study indices:      ", paste(x$original_indices, collapse = ", "), "\n")
  if (is.finite(x$N))     cat("  N:                           ", signif(x$N, 4), "\n")
  if (is.finite(x$Delta)) cat("  Delta:                       ", signif(x$Delta, 4), "\n")
  invisible(x)
}
