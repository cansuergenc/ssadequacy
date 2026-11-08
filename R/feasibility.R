# Calibrated weak-identification hazard model (Study 2):
#   -log P(admissible) = exp(c0) * N^c1 * S^c2 * F^c3,
# where S is the mean per-factor sum of squared standardized loadings and F the number of factors.
FEAS_COEF <- c(c0 = 11.7520, c1 = -3.6258, c2 = -5.3081, c3 = 1.4854)

#' Identification strength and predicted admissibility for CFA designs
#'
#' Implements the calibrated weak-identification scaling law for the rate of
#' admissible (converged, Heywood-free) maximum likelihood solutions:
#' \deqn{-\log P(\mathrm{admissible}) \approx A\, F^{c_3} (N S^{c_2/c_1})^{c_1},}
#' with \eqn{S} the effective per-factor sum of squared standardized loadings (a power mean
#' over factors that is dominated by the weakest factor, because hazards add across factors) and \eqn{F}
#' the number of factors. The law was calibrated on simple-structure models with
#' 3 to 5 factors, 3 to 5 indicators per factor, loadings .35 to .60, normal data and
#' ML estimation; outside that range use \code{\link{ssa_simulate}}.
#'
#' @param design An \code{ssa_design} of type \code{"cfa"}.
#' @param N Sample size(s).
#' @param kappa Target rate of admissible solutions (for \code{ssa_n_kappa}).
#' @return \code{ssa_feasibility}: predicted admissibility rate at each \code{N};
#'   \code{ssa_n_kappa}: the smallest N with predicted rate \code{kappa};
#'   \code{ssa_identification_strength}: the index \eqn{N S^{c_2/c_1}} at each \code{N}.
#' @examples
#' d <- ssa_design_cfa(nfactors = 3, nind = 5, loadings = c(.35, .40))
#' ssa_feasibility(d, N = c(100, 150, 200))
#' ssa_n_kappa(d, kappa = .95)
#' ssa_identification_strength(d, N = 150)
#' @export
ssa_feasibility <- function(design, N) {
  check_design(design, "cfa"); check_positive(N, "N")
  S <- design_strength(design)
  exp(-exp(FEAS_COEF[["c0"]] + FEAS_COEF[["c1"]] * log(N) + FEAS_COEF[["c2"]] * log(S) +
             FEAS_COEF[["c3"]] * log(design$nfactors)))
}

#' @rdname ssa_feasibility
#' @export
ssa_n_kappa <- function(design, kappa = .95) {
  check_design(design, "cfa")
  if (!is.numeric(kappa) || length(kappa) != 1 || kappa <= 0 || kappa >= 1) stop("'kappa' must be a single value in (0, 1).")
  S <- design_strength(design)
  unname(exp((log(-log(kappa)) - FEAS_COEF[["c0"]] - FEAS_COEF[["c2"]] * log(S) -
                FEAS_COEF[["c3"]] * log(design$nfactors)) / FEAS_COEF[["c1"]]))
}

#' @rdname ssa_feasibility
#' @export
ssa_identification_strength <- function(design, N) {
  check_design(design, "cfa"); check_positive(N, "N")
  N * design_strength(design)^(FEAS_COEF[["c2"]] / FEAS_COEF[["c1"]])
}

# Effective per-factor strength: hazards add over factors, h = sum_f A N^c1 S_f^c2, so the
# effective S is the power mean of the factor strengths with exponent c2 (dominated by the weakest factor).
design_strength <- function(design) {
  lam <- matrix(design$theta, design$nind)
  Sf <- colSums(lam^2)
  mean(Sf^FEAS_COEF[["c2"]])^(1 / FEAS_COEF[["c2"]])
}

check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) < 1 || any(!is.finite(x)) || any(x <= 0))
    stop(sprintf("'%s' must be numeric and positive.", name))
  invisible(x)
}
