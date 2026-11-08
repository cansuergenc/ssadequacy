#' Fleishman power-method coefficients
#'
#' Solves the Fleishman (1978) system for the coefficients of
#' \eqn{Y = a + bZ + cZ^2 + dZ^3} that yield a standardized variable with the
#' requested skewness and excess kurtosis.
#'
#' @param skew Target skewness.
#' @param exkurt Target excess kurtosis.
#' @return Named numeric vector \code{c(a, b, c, d)}.
#' @references Fleishman, A. I. (1978). A method for simulating non-normal distributions.
#'   \emph{Psychometrika, 43}(4), 521--532. \doi{10.1007/BF02293811}
#' @examples
#' fleishman_coef(2, 7)
#' @export
fleishman_coef <- function(skew = 0, exkurt = 0) {
  if (skew == 0 && exkurt == 0) return(c(a = 0, b = 1, c = 0, d = 0))
  eqs <- function(p) {
    b <- p[1]; c <- p[2]; d <- p[3]
    c(b^2 + 6 * b * d + 2 * c^2 + 15 * d^2 - 1,
      2 * c * (b^2 + 24 * b * d + 105 * d^2 + 2) - skew,
      24 * (b * d + c^2 * (1 + b^2 + 28 * b * d) +
              d^2 * (12 + 48 * b * d + 141 * c^2 + 225 * d^2)) - exkurt)
  }
  newton <- function(p) {
    for (it in 1:200) {
      f0 <- eqs(p); if (max(abs(f0)) < 1e-12) break
      J <- vapply(1:3, function(i) { h <- 1e-7; e <- replace(numeric(3), i, h); (eqs(p + e) - f0) / h }, numeric(3))
      step <- tryCatch(solve(J, f0), error = function(e) rep(NA, 3))
      if (anyNA(step)) return(NULL)
      p <- p - step
    }
    if (max(abs(eqs(p))) < 1e-9) p else NULL
  }
  p <- NULL
  for (start in list(c(.9, .3, 0), c(1, .1, 0), c(.8, .2, .05), c(.7, .4, .02))) {
    p <- newton(start); if (!is.null(p) && p[1] > 0) break
  }
  if (is.null(p)) stop("Fleishman system has no solution for these moments.")
  c(a = -p[2], b = p[1], c = p[2], d = p[3])
}

fleishman_transform <- function(z, cf) cf[1] + cf[2] * z + cf[3] * z^2 + cf[4] * z^3

# Vale & Maurelli (1983) intermediate correlation for a common Fleishman transform
vm_intermediate <- function(R, cf) {
  b <- cf[2]; c <- cf[3]; d <- cf[4]
  g <- function(r) r * (b^2 + 6 * b * d + 9 * d^2) + r^2 * 2 * c^2 + r^3 * 6 * d^2
  Ri <- diag(nrow(R))
  for (i in seq_len(nrow(R))) for (j in seq_len(i - 1)) {
    Ri[i, j] <- Ri[j, i] <- stats::uniroot(function(r) g(r) - R[i, j], c(-1, 1), tol = 1e-12)$root
  }
  Ri
}

#' Generate multivariate (non-)normal data with a given covariance matrix
#'
#' Uses the Vale-Maurelli (1983) extension of the Fleishman method so that the
#' generated variables reproduce \code{Sigma} exactly in the population while
#' every margin has the requested skewness and excess kurtosis.
#'
#' @param n Number of observations.
#' @param Sigma Population covariance matrix.
#' @param skew,exkurt Marginal skewness and excess kurtosis (common to all variables).
#' @return An \code{n x p} matrix with column names taken from \code{Sigma}.
#' @references Vale, C. D., & Maurelli, V. A. (1983). Simulating multivariate nonnormal
#'   distributions. \emph{Psychometrika, 48}(3), 465--471. \doi{10.1007/BF02293687}
#' @examples
#' S <- matrix(c(1, .5, .5, 1), 2)
#' X <- rmvnonnorm(1000, S, skew = 2, exkurt = 7)
#' cor(X)
#' @export
rmvnonnorm <- function(n, Sigma, skew = 0, exkurt = 0) {
  p <- nrow(Sigma); sds <- sqrt(diag(Sigma)); R <- stats::cov2cor(Sigma)
  if (skew == 0 && exkurt == 0) {
    X <- matrix(stats::rnorm(n * p), n) %*% chol(Sigma)
  } else {
    cf <- fleishman_coef(skew, exkurt)
    Ri <- vm_intermediate(R, cf)
    Z <- matrix(stats::rnorm(n * p), n) %*% chol(Ri)
    X <- sweep(fleishman_transform(Z, cf), 2, sds, "*")
  }
  colnames(X) <- colnames(Sigma) %||% paste0("x", seq_len(p))
  X
}

`%||%` <- function(a, b) if (is.null(a)) b else a
