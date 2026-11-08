#' Design constructors
#'
#' A design bundles a population model, a data generator, an estimator and
#' the focal parameters on which all four SSA criteria are evaluated
#' (estimand coherence: power and precision always refer to the same
#' parameters).
#'
#' @name ssa_design
#' @return An object of class \code{ssa_design}: a list holding the population model, the
#'   focal parameter values (\code{theta}) and the data-generation settings.
#' @examples
#' ssa_design_cfa(nfactors = 3, nind = 5, loadings = c(.55, .60))
#' ssa_design_regression(k = 10, R2 = .25, skew = 2, exkurt = 7)
#' ssa_design_anova(k = 3, f = .25)
#' \donttest{
#' fit <- lavaan::cfa("visual =~ x1 + x2 + x3\n textual =~ x4 + x5 + x6\n speed =~ x7 + x8 + x9",
#'                    data = lavaan::HolzingerSwineford1939, std.lv = TRUE)
#' ssa_design_from_lavaan(fit)
#' }
NULL

#' @describeIn ssa_design Confirmatory factor analysis with simple structure.
#' @param nfactors Number of correlated factors.
#' @param nind Number of indicators per factor.
#' @param loadings Standardized loadings: either one value, a range \code{c(lo, hi)}
#'   (spread evenly within each factor) or a full vector of length \code{nfactors * nind}.
#' @param phi Common inter-factor correlation, or a full factor correlation matrix.
#' @param skew,exkurt Marginal skewness and excess kurtosis of the indicators.
#' @param estimator \code{"ML"} or \code{"MLR"}; defaults to \code{"MLR"} for non-normal data.
#' @export
ssa_design_cfa <- function(nfactors = 3, nind = 5, loadings = c(.55, .60), phi = .30,
                           skew = 0, exkurt = 0, estimator = NULL) {
  check_count(nfactors, "nfactors", 1); check_count(nind, "nind", 2)
  p <- nfactors * nind
  if (!is.numeric(loadings) || !length(loadings) %in% c(1, 2, p))
    stop("'loadings' must have length 1 (common value), 2 (range c(lo, hi)) or nfactors * nind (full vector).")
  if (any(!is.finite(loadings)) || any(abs(loadings) >= 1))
    stop("Standardized loadings must be finite and lie strictly between -1 and 1.")
  if (is.matrix(phi)) {
    if (!all(dim(phi) == nfactors)) stop("'phi' must be a scalar or an nfactors x nfactors correlation matrix.")
  } else if (length(phi) != 1 || abs(phi) >= 1) stop("'phi' must be a single correlation in (-1, 1) or a matrix.")
  if (!is.null(estimator) && !estimator %in% c("ML", "MLR")) stop("'estimator' must be \"ML\" or \"MLR\".")
  lam <- if (length(loadings) == p) loadings else
    rep(if (length(loadings) == 2) seq(loadings[1], loadings[2], length.out = nind) else rep(loadings, nind), nfactors)
  L <- matrix(0, p, nfactors)
  for (j in seq_len(nfactors)) L[(j - 1) * nind + seq_len(nind), j] <- lam[(j - 1) * nind + seq_len(nind)]
  Phi <- if (is.matrix(phi)) phi else { P <- matrix(phi, nfactors, nfactors); diag(P) <- 1; P }
  Sigma <- L %*% Phi %*% t(L); diag(Sigma) <- 1
  vn <- paste0("x", seq_len(p)); dimnames(Sigma) <- list(vn, vn)
  model <- paste(vapply(seq_len(nfactors), function(j)
    paste0("f", j, " =~ ", paste(vn[(j - 1) * nind + seq_len(nind)], collapse = " + ")), ""), collapse = "\n")
  if (is.null(estimator)) estimator <- if (skew == 0 && exkurt == 0) "ML" else "MLR"
  fc <- fleishman_coef(skew, exkurt)
  Ri <- if (skew == 0 && exkurt == 0) stats::cov2cor(Sigma) else vm_intermediate(stats::cov2cor(Sigma), fc)
  structure(list(type = "cfa", fc = fc, cholR = chol(Ri), model = model, Sigma = Sigma, theta = lam, skew = skew, exkurt = exkurt,
                 estimator = estimator, nfactors = nfactors, nind = nind, phi = Phi,
                 label = sprintf("CFA %dF x %d, lambda %.2f-%.2f, %s", nfactors, nind, min(lam), max(lam),
                                 if (skew == 0 && exkurt == 0) "normal" else "non-normal")),
            class = "ssa_design")
}

#' @describeIn ssa_design OLS multiple regression with equicorrelated predictors and equal standardized coefficients.
#' @param k Number of predictors (regression) or groups (ANOVA).
#' @param R2 Population coefficient of determination.
#' @param rho Common predictor intercorrelation.
#' @export
ssa_design_regression <- function(k = 3, R2 = .10, rho = .30, skew = 0, exkurt = 0) {
  check_count(k, "k", 1)
  if (!is.numeric(R2) || length(R2) != 1 || R2 <= 0 || R2 >= 1) stop("'R2' must be a single value in (0, 1).")
  if (!is.numeric(rho) || length(rho) != 1 || rho <= -1 / max(k - 1, 1) || rho >= 1)
    stop("'rho' must be a single value giving a positive definite equicorrelation matrix.")
  Sx <- matrix(rho, k, k); diag(Sx) <- 1
  beta <- rep(sqrt(R2 / sum(Sx)), k)
  fc <- fleishman_coef(skew, exkurt)
  Ri <- if (skew == 0 && exkurt == 0) Sx else vm_intermediate(Sx, fc)
  structure(list(type = "regression", fc = fc, cholR = chol(Ri), k = k, R2 = R2, rho = rho, Sx = Sx, beta = beta, theta = beta,
                 skew = skew, exkurt = exkurt,
                 label = sprintf("Regression k = %d, R2 = %.2f, %s", k, R2,
                                 if (skew == 0 && exkurt == 0) "normal" else "non-normal")),
            class = "ssa_design")
}

#' @describeIn ssa_design One-way between-subjects ANOVA with equally spaced means; the focal
#'   parameter is the standardized difference between the two extreme groups.
#' @param f Cohen's f (SD of group means over the within-group SD).
#' @export
ssa_design_anova <- function(k = 2, f = .25, skew = 0, exkurt = 0) {
  check_count(k, "k", 2)
  if (!is.numeric(f) || length(f) != 1 || !is.finite(f) || f <= 0) stop("'f' must be a single positive value.")
  m <- seq(-1, 1, length.out = k); m <- m / sqrt(mean((m - mean(m))^2)) * f
  structure(list(type = "anova", fc = fleishman_coef(skew, exkurt), k = k, f = f, means = m, theta = max(m) - min(m),
                 skew = skew, exkurt = exkurt,
                 label = sprintf("ANOVA k = %d, f = %.2f, %s", k, f,
                                 if (skew == 0 && exkurt == 0) "normal" else "non-normal")),
            class = "ssa_design")
}

#' @describeIn ssa_design Planning population taken from a fitted 'lavaan' CFA (e.g., a pilot or
#'   published study): uses its standardized loadings and factor correlations. Requires
#'   simple structure with the same number of indicators per factor.
#' @param fit A fitted \code{lavaan} CFA object.
#' @export
ssa_design_from_lavaan <- function(fit, skew = 0, exkurt = 0, estimator = NULL) {
  ss <- lavaan::standardizedSolution(fit)
  ld <- ss[ss$op == "=~", ]
  fac <- unique(ld$lhs); nind <- table(ld$lhs)[fac]
  if (length(unique(nind)) != 1) stop("All factors must have the same number of indicators.")
  Phi <- lavaan::lavInspect(fit, "cor.lv")[fac, fac]
  d <- ssa_design_cfa(length(fac), nind[[1]], loadings = ld$est.std, phi = unname(Phi),
                      skew = skew, exkurt = exkurt, estimator = estimator)
  d$label <- sprintf("CFA from fitted model (%d factors x %d), lambda %.2f-%.2f, %s", length(fac), nind[[1]],
                     min(ld$est.std), max(ld$est.std), if (skew == 0 && exkurt == 0) "normal" else "non-normal")
  d
}

check_count <- function(x, name, min) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x != round(x) || x < min)
    stop(sprintf("'%s' must be a single integer >= %d.", name, min))
  invisible(x)
}

check_design <- function(design, type = NULL) {
  if (!inherits(design, "ssa_design")) stop("'design' must be an object created by ssa_design_cfa(), ssa_design_regression(), ssa_design_anova() or ssa_design_from_lavaan().")
  if (!is.null(type) && design$type != type) stop(sprintf("'design' must be an ssa_design of type \"%s\".", type))
  invisible(design)
}

#' @export
print.ssa_design <- function(x, ...) {
  cat("<ssa_design>", x$label, "\n")
  cat("  focal parameters:", length(x$theta), " (mean |theta| =", round(mean(abs(x$theta)), 3), ")\n")
  invisible(x)
}
