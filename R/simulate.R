#' Monte Carlo SSA analysis with root finding
#'
#' Simulates the design on a grid of sample sizes, smooths each criterion with
#' the functional form implied by asymptotic theory, solves for the sample size
#' at which each criterion is met, and reports \eqn{N^*} and the binding
#' criterion.
#'
#' Smoothing models: power, \eqn{\Phi^{-1}(\pi_N) = a + b\sqrt{N}};
#' width, \eqn{\log W_N = a + b \log N}; assurance and feasibility,
#' \eqn{\mathrm{logit}(p_N) = a + b \log N}; relative bias, \eqn{\mathrm{RelBias}_N = c/N}.
#'
#' @param design An \code{ssa_design}.
#' @param delta Target CI width.
#' @param power,alpha Power target and significance level.
#' @param assurance Optional assurance probability for the width criterion.
#' @param tau Maximum absolute relative bias.
#' @param kappa Minimum rate of admissible solutions.
#' @param n_grid Sample sizes to simulate; by default centred on the analytic \eqn{N^*}.
#' @param reps Replications per grid point.
#' @param cores Cores for \code{parallel::mclapply}.
#' @param verify Logical; re-simulate at the solved \eqn{N^*} and report all criteria there.
#' @return An object of class \code{ssa_sim}.
#' @examples
#' \donttest{
#' s <- ssa_simulate(ssa_design_regression(k = 3, R2 = .10), delta = .25, reps = 300)
#' s
#' ssa_report(s)
#' }
#' @export
ssa_simulate <- function(design, delta, power = .80, alpha = .05, assurance = NULL, tau = .05,
                         kappa = .95, n_grid = NULL, reps = 500, cores = 1, verify = TRUE) {
  check_design(design); check_positive(delta, "delta"); check_positive(tau, "tau")
  check_prob(power, "power"); check_prob(alpha, "alpha"); check_prob(kappa, "kappa")
  if (!is.null(assurance)) check_prob(assurance, "assurance")
  if (!is.null(n_grid)) {
    check_positive(n_grid, "n_grid")
    if (length(n_grid) < 2) stop("'n_grid' must contain at least two sample sizes.")
    n_grid <- sort(unique(round(n_grid)))
  }
  if (is.null(n_grid)) {
    a <- ssa_analytic(design, delta, power, alpha, assurance)
    centre <- max(a$N_star, 4 * design_min_n(design))
    n_grid <- unique(pmax(design_min_n(design) + 2, round(centre * c(.35, .5, .7, 1, 1.4, 2))))
  }
  cells <- lapply(n_grid, function(n) ssa_cell(design, n, reps, alpha, cores))
  tab <- t(vapply(cells, function(cl) suppressWarnings(summary(cl, delta = delta)), numeric(8)))
  tab <- data.frame(N = n_grid, tab)
  req <- solve_criteria(tab, delta, power, tau, kappa, assurance, reps)
  # NA means "met at every simulated N"; if that holds for all criteria no criterion binds
  none <- all(is.na(req))
  out <- structure(list(design = design, table = tab, N = req,
                        N_star = if (none) NA_integer_ else as.integer(ceiling(max(req, na.rm = TRUE))),
                        binding = if (none) "none" else names(req)[which.max(replace(req, is.na(req), -Inf))],
                        delta = delta, power = power, alpha = alpha,
                        assurance = assurance, tau = tau, kappa = kappa, reps = reps),
                   class = "ssa_sim")
  if (none) message("Every criterion is met at the smallest simulated N (", min(n_grid),
                    "); choose a lower n_grid to locate N*.")
  if (verify && !is.na(out$N_star) && out$N_star >= design_min_n(design)) {
    vc <- ssa_cell(design, out$N_star, reps, alpha, cores)
    out$verification <- summary(vc, delta = delta)
  }
  out
}

# Solve each criterion on the simulated grid using theory-based smoothing.
solve_criteria <- function(tab, delta, power, tau, kappa, assurance, reps) {
  N <- tab$N
  clamp <- function(p) pmin(pmax(p, .5 / reps), 1 - .5 / reps)
  root_lin <- function(x, y, target, tf) {
    fit <- stats::lm(y ~ x)
    b <- stats::coef(fit)
    if (!is.finite(b[2]) || b[2] == 0) return(NA)
    unname(tf((target - b[1]) / b[2]))
  }
  # power: probit linear in sqrt(N)
  # local root: interpolate the transformed criterion between the two grid points that bracket the target
  root_local <- function(x, y, target, tf, increasing = TRUE) {
    ok <- is.finite(y); x <- x[ok]; y <- y[ok]
    if (length(x) < 2) return(NA_real_)
    above <- if (increasing) y >= target else y <= target
    j <- which(above)[1]
    if (is.na(j)) return(root_lin(x, y, target, tf))
    if (j == 1) return(root_lin(x[1:2], y[1:2], target, tf))
    unname(tf(x[j - 1] + (target - y[j - 1]) * (x[j] - x[j - 1]) / (y[j] - y[j - 1])))
  }
  y_pi <- stats::qnorm(clamp(tab$power))
  N_pi <- if (all(tab$power >= power)) NA_real_
          else root_local(sqrt(N), y_pi, stats::qnorm(power), function(s) s^2)
  # width: log-log (NA, like the other criteria, when the target is met at every grid point)
  N_delta <- if (all(tab$width <= delta)) NA_real_ else root_lin(log(N), log(tab$width), log(delta), exp)
  if (!is.null(assurance)) {
    N_delta <- if (all(tab$assurance >= assurance)) NA_real_ else
      root_local(log(N), stats::qlogis(clamp(tab$assurance)), stats::qlogis(assurance), exp)
  }
  # feasibility
  N_kappa <- if (all(tab$feasibility >= kappa)) NA_real_ else {
    root_local(log(N), stats::qlogis(clamp(tab$feasibility)), stats::qlogis(kappa), exp)
  }
  # stability: RelBias ~ c / N through the origin
  cb <- sum(tab$rel_bias / N) / sum(1 / N^2)
  N_tau <- if (all(abs(tab$rel_bias) <= tau)) NA_real_ else abs(cb) / tau
  c(Detectability = N_pi, Precision = N_delta, Stability = N_tau, Feasibility = N_kappa)
}

#' @export
print.ssa_sim <- function(x, ...) {
  cat("<ssa_sim>", x$design$label, "\n")
  cat(sprintf("  thresholds: power >= %.2f, CI width <= %.3f%s, |rel. bias| <= %.2f, admissible >= %.2f\n",
              x$power, x$delta, if (is.null(x$assurance)) "" else sprintf(" (assurance %.2f)", x$assurance),
              x$tau, x$kappa))
  print(round(x$table, 3), row.names = FALSE)
  for (k in names(x$N)) cat(sprintf("  N required by %-13s %s\n", k,
    if (is.na(x$N[k])) sprintf("< %d (met at every grid point)", min(x$table$N)) else sprintf("%8.0f", ceiling(x$N[k]))))
  if (is.na(x$N_star)) cat(sprintf("  N* < %d: every criterion is met at the smallest simulated N; no criterion binds\n", min(x$table$N)))
  else cat(sprintf("  N* = %d, binding criterion: %s\n", x$N_star, x$binding))
  if (!is.null(x$verification)) {
    cat("  verification at N*:\n"); print(round(x$verification, 3))
  }
  invisible(x)
}
