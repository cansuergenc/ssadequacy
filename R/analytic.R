#' Power-precision crossover ratio
#'
#' Under unit-information asymptotics the sample sizes required for a
#' Wald test with power \eqn{\pi} and for an expected two-sided
#' \eqn{(1-\alpha)} confidence interval of width \eqn{\delta} on the same
#' parameter \eqn{\theta} satisfy
#' \deqn{N_\delta / N_\pi = [ r^* / r ]^2, \quad r = \delta/|\theta|,\quad
#'       r^* = 2 z_{1-\alpha/2} / (z_{1-\alpha/2} + z_\pi).}
#' Precision is the binding criterion if and only if \eqn{r < r^*}.
#' The ratio does not depend on the unit variance of the estimator, so method,
#' model size and distribution cannot change which of the two criteria binds.
#'
#' @param alpha Two-sided significance level.
#' @param power Target power.
#' @return The crossover ratio \eqn{r^*}.
#' @examples
#' ssa_crossover()            # 1.399
#' ssa_crossover(power = .90) # 1.209
#' @export
ssa_crossover <- function(alpha = .05, power = .80) {
  z <- stats::qnorm(1 - alpha / 2)
  2 * z / (z + stats::qnorm(power))
}

#' Unit (per-observation) asymptotic variance of the focal estimators
#'
#' Returns \eqn{v_j} such that \eqn{\mathrm{Var}(\hat\theta_j) \approx v_j / N}.
#' Closed forms are used for regression and ANOVA. For CFA the variance is taken
#' from a 'lavaan' fit to the population covariance matrix (normal data, ML) or,
#' for non-normal data, from a robust (MLR) fit to one very large generated sample.
#'
#' @param design An \code{ssa_design}.
#' @param M Size of the population-scale sample used for non-normal CFA designs.
#' @return Numeric vector of unit variances, one per focal parameter.
#' @examples
#' ssa_unit_variance(ssa_design_regression(k = 3, R2 = .10))
#' ssa_unit_variance(ssa_design_anova(k = 2, f = .25))
#' \donttest{
#' ssa_unit_variance(ssa_design_cfa(3, 5, loadings = c(.55, .60)))
#' }
#' @export
ssa_unit_variance <- function(design, M = 2e5) {
  switch(design$type,
    regression = rep((1 - design$R2) * diag(solve(design$Sx)), length.out = design$k),
    anova = 2 * design$k,
    cfa = {
      if (normal_design(design) && design$estimator == "ML") {
        fit <- lavaan::cfa(design$model, sample.cov = design$Sigma, sample.nobs = 1e6, std.lv = TRUE,
                           sample.cov.rescale = FALSE)
        n_eff <- 1e6
      } else {
        X <- gen_mv(M, design); colnames(X) <- colnames(design$Sigma)
        fit <- lavaan::cfa(design$model, data = as.data.frame(X), std.lv = TRUE, estimator = design$estimator)
        n_eff <- M
      }
      ss <- lavaan::standardizedSolution(fit, type = "std.all")
      ss$se[ss$op == "=~"]^2 * n_eff
    })
}

#' Closed-form bottleneck analysis
#'
#' Computes the sample size required by each SSA criterion from unit-information
#' asymptotics and identifies the binding criterion.
#'
#' @param design An \code{ssa_design}.
#' @param delta Target confidence interval width (same metric as the focal parameter).
#' @param power Target power \eqn{\pi} for the focal Wald tests.
#' @param alpha Two-sided significance level.
#' @param assurance Optional probability \eqn{\gamma} that the realized width does not exceed \code{delta}.
#' @param se_cv Unit coefficient of variation of the estimated standard error,
#'   \eqn{\sqrt{N}\,\mathrm{CV}(\widehat{SE})}; \eqn{1/\sqrt{2}} is exact for regression and ANOVA.
#' @param tau Maximum relative bias.
#' @param bias_const Second-order bias constant \eqn{b} in \eqn{\mathrm{RelBias} \approx b/N} (optional).
#' @param kappa Target rate of admissible solutions (CFA designs; uses the calibrated law in \code{\link{ssa_feasibility}}).
#' @param v Optional unit variances (computed with \code{\link{ssa_unit_variance}} if omitted).
#' @return An object of class \code{ssa_analytic}.
#' @references
#' Goodman, S. N., & Berlin, J. A. (1994). The use of predicted confidence intervals when
#' planning experiments and the misuse of power when interpreting results.
#' \emph{Annals of Internal Medicine, 121}(3), 200--206. \doi{10.7326/0003-4819-121-3-199408010-00008}
#'
#' Kelley, K., & Rausch, J. R. (2006). Sample size planning for the standardized mean
#' difference: Accuracy in parameter estimation via narrow confidence intervals.
#' \emph{Psychological Methods, 11}(4), 363--385. \doi{10.1037/1082-989X.11.4.363}
#' @examples
#' # Regression: power and precision for the standardized coefficients
#' ssa_analytic(ssa_design_regression(k = 3, R2 = .10), delta = .25)
#' # Two-group comparison, d = .50, with assurance for the interval width
#' ssa_analytic(ssa_design_anova(k = 2, f = .25), delta = .50, assurance = .80)
#' \donttest{
#' # Factor model with medium loadings
#' ssa_analytic(ssa_design_cfa(3, 5, loadings = c(.55, .60)), delta = .20)
#' }
#' @export
ssa_analytic <- function(design, delta, power = .80, alpha = .05, assurance = NULL,
                         se_cv = 1 / sqrt(2), tau = .05, bias_const = NULL, kappa = .95, v = NULL) {
  check_design(design)
  check_positive(delta, "delta"); check_positive(tau, "tau")
  check_prob(power, "power"); check_prob(alpha, "alpha"); check_prob(kappa, "kappa")
  if (!is.null(assurance)) check_prob(assurance, "assurance")
  if (is.null(v)) v <- ssa_unit_variance(design)
  theta <- abs(design$theta); v <- rep(v, length.out = length(theta))
  z <- stats::qnorm(1 - alpha / 2)
  df_fun <- switch(design$type, regression = function(N) N - design$k - 1,
                   anova = function(N) N - design$k, cfa = function(N) Inf)
  crit <- function(N) if (is.finite(df_fun(N))) stats::qt(1 - alpha / 2, df_fun(N)) else z
  pow <- function(N) mean(stats::pnorm(theta * sqrt(N / v) - crit(N)))
  wid <- function(N) 2 * crit(N) * mean(sqrt(v / N))
  lo <- design_min_n(design)
  # smallest N with g(N) <= 0 for a decreasing function g
  solve_dec <- function(g) if (g(lo) <= 0) lo else stats::uniroot(g, c(lo, 1e9), tol = 1e-6)$root
  N_pi <- solve_dec(function(N) power - pow(N))
  N_delta <- solve_dec(function(N) wid(N) - delta)
  N_delta_g <- if (is.null(assurance)) NA else
    solve_dec(function(N) wid(N) * (1 + stats::qnorm(assurance) * se_cv / sqrt(N)) - delta)
  N_tau <- if (is.null(bias_const)) NA else abs(bias_const) / tau
  r <- delta / mean(theta)
  N_kappa <- if (design$type == "cfa") ssa_n_kappa(design, kappa) else NA
  req <- c(Detectability = N_pi, Precision = if (is.na(N_delta_g)) N_delta else N_delta_g, Stability = N_tau,
           Feasibility = N_kappa)
  structure(list(design = design, v = v, delta = delta, power = power, alpha = alpha, assurance = assurance,
                 r = r, r_star = ssa_crossover(alpha, power), N = req, N_delta_expected = N_delta,
                 binding = names(which.max(replace(req, is.na(req), -Inf))), N_star = max(req, na.rm = TRUE)),
            class = "ssa_analytic")
}

design_min_n <- function(design) switch(design$type, regression = design$k + 3, anova = 2 * design$k + 1, cfa = 10)

check_prob <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0 || x >= 1)
    stop(sprintf("'%s' must be a single value strictly between 0 and 1.", name))
  invisible(x)
}

#' @export
print.ssa_analytic <- function(x, ...) {
  cat("<ssa_analytic>", x$design$label, "\n")
  cat(sprintf("  r = delta/|theta| = %.3f   crossover r* = %.3f  -> %s binds between power and precision\n",
              x$r, x$r_star, if (x$r < x$r_star) "Precision" else "Detectability"))
  for (k in names(x$N)) if (!is.na(x$N[k])) cat(sprintf("  N required by %-13s %8.0f\n", k, ceiling(x$N[k])))
  cat(sprintf("  N* (analytic) = %.0f, binding criterion: %s\n", ceiling(x$N_star), x$binding))
  if (x$design$type == "cfa")
    cat("  Feasibility uses a calibrated weak-identification law (see ?ssa_feasibility); verify with ssa_simulate().\n")
  invisible(x)
}
