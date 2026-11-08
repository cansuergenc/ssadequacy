#' Sample size justification text
#'
#' Produces a paragraph suitable for a preregistration or a Method section
#' from an \code{ssa_analytic} or \code{ssa_sim} result.
#'
#' @param x An \code{ssa_analytic} or \code{ssa_sim} object.
#' @param verbose Logical; print the paragraph to the console (default) or only return it.
#' @param ... Unused.
#' @return A character string (invisibly when \code{verbose = TRUE}).
#' @examples
#' ssa_report(ssa_analytic(ssa_design_anova(k = 2, f = .25), delta = .5))
#' txt <- ssa_report(ssa_analytic(ssa_design_anova(k = 2, f = .25), delta = .5), verbose = FALSE)
#' @export
ssa_report <- function(x, verbose = TRUE, ...) {
  if (!inherits(x, c("ssa_analytic", "ssa_sim"))) stop("'x' must be an ssa_analytic or ssa_sim object.")
  d <- x$design
  crit_txt <- sprintf("power of at least %.2f for two-sided Wald tests (alpha = %.2f) of each focal parameter, an expected 95%% confidence interval width of at most %.3f%s",
                      x$power, x$alpha, x$delta,
                      if (is.null(x$assurance)) "" else sprintf(" with assurance %.2f", x$assurance))
  focal <- switch(d$type, cfa = "the standardized factor loadings", regression = "the standardized regression coefficients",
                  anova = "the standardized difference between the extreme group means")
  if (inherits(x, "ssa_analytic")) {
    txt <- sprintf(paste(
      "Sample size was planned with the Sample Size Adequacy (SSA) framework (Ergen\u00e7 \u00d6zda\u015f, in preparation) for %s.",
      "The focal parameters were %s (mean |theta| = %.3f). We required %s.",
      "Under unit-information asymptotics these requirements imply N = %.0f for detectability and N = %.0f for precision;",
      "because the ratio of the target width to the effect size (r = %.2f) is %s the crossover value r* = %.2f, %s is the binding criterion and the minimum adequate sample size is N* = %.0f."),
      d$label, focal, mean(abs(d$theta)), crit_txt, ceiling(x$N["Detectability"]), ceiling(x$N["Precision"]),
      x$r, if (x$r < x$r_star) "below" else "above", x$r_star, tolower(x$binding), ceiling(x$N_star))
    if (d$type == "cfa") txt <- paste(txt, "Because admissibility of latent variable solutions is not governed by regular asymptotics, this value should be verified by simulation (ssa_simulate()).")
  } else {
    v <- x$verification
    txt <- sprintf(paste(
      "Sample size was planned with the Sample Size Adequacy (SSA) framework for %s, using %d Monte Carlo replications per candidate N.",
      "The focal parameters were %s. We required %s, an absolute relative bias of at most %.2f, and at least %.0f%% admissible solutions.",
      "%s%s"),
      d$label, x$reps, focal, crit_txt, x$tau, 100 * x$kappa,
      if (is.na(x$N_star)) sprintf("Every criterion was met at the smallest simulated sample size (N = %d), so the minimum adequate sample size is below this value.", min(x$table$N))
      else sprintf("The minimum adequate sample size was N* = %d, with %s as the binding criterion.", x$N_star, tolower(x$binding)),
      if (is.null(v)) "" else sprintf(" At N* the simulated power was %.2f, the mean interval width %.3f, the relative bias %.3f, and the rate of admissible solutions %.2f.",
                                      v["power"], v["width"], v["rel_bias"], v["feasibility"]))
  }
  txt <- paste(strwrap(txt, 90), collapse = "\n")
  if (!verbose) return(txt)
  cat(txt, "\n")
  invisible(txt)
}

#' Launch the interactive SSA planner
#'
#' @param ... Passed to \code{shiny::runApp}.
#' @return Called for its side effect of launching the app; returns the value of \code{shiny::runApp()}.
#' @examples
#' if (interactive()) run_ssa_app()
#' @export
run_ssa_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) stop("Install the 'shiny' package to use the app.")
  shiny::runApp(system.file("shiny", package = "ssadequacy"), ...)
}
