# One Monte Carlo replication for a design at sample size n.
# Returns a list: adm (0/1), est, se, crit (critical value used for Wald CI/test), omni (0/1 or NA)

ssa_replicate <- function(design, n, alpha = .05) {
  switch(design$type,
         cfa = rep_cfa(design, n, alpha),
         regression = rep_reg(design, n, alpha),
         anova = rep_anova(design, n, alpha),
         stop("Unknown design type"))
}

rep_cfa <- function(d, n, alpha) {
  X <- gen_mv(n, d); colnames(X) <- colnames(d$Sigma)
  fail <- list(adm = 0, est = rep(NA_real_, length(d$theta)), se = rep(NA_real_, length(d$theta)),
               crit = stats::qnorm(1 - alpha / 2), omni = NA_real_)
  fit <- tryCatch(suppressWarnings(lavaan::cfa(d$model, data = as.data.frame(X), std.lv = TRUE,
                                                estimator = d$estimator, warn = FALSE,
                                                control = list(iter.max = 1000))),
                  error = function(e) NULL)
  if (is.null(fit) || !lavaan::lavInspect(fit, "converged")) return(fail)
  pe <- lavaan::parameterEstimates(fit, se = FALSE, zstat = FALSE, pvalue = FALSE, ci = FALSE)
  resid <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% colnames(d$Sigma)]
  lvcov <- lavaan::lavInspect(fit, "cov.lv")
  if (any(resid <= 0) || min(eigen(lvcov, symmetric = TRUE, only.values = TRUE)$values) <= 0) return(fail)
  ss <- tryCatch(lavaan::standardizedSolution(fit, type = "std.all", zstat = FALSE, pvalue = FALSE, ci = FALSE),
                 error = function(e) NULL)
  if (is.null(ss)) return(fail)
  ld <- ss[ss$op == "=~", ]
  if (any(abs(ld$est.std) >= 1) || any(!is.finite(ld$se))) return(fail)
  list(adm = 1, est = ld$est.std, se = ld$se, crit = stats::qnorm(1 - alpha / 2), omni = NA_real_)
}

normal_design <- function(d) d$skew == 0 && d$exkurt == 0

err_draw <- function(n, d) {
  z <- stats::rnorm(n)
  if (normal_design(d)) z else fleishman_transform(z, d$fc)
}

# correlated standardized variables (all designs here have unit variances)
gen_mv <- function(n, d) {
  Z <- matrix(stats::rnorm(n * ncol(d$cholR)), n) %*% d$cholR
  if (normal_design(d)) Z else fleishman_transform(Z, d$fc)
}

rep_reg <- function(d, n, alpha) {
  X <- gen_mv(n, d)
  y <- drop(X %*% d$beta) + sqrt(1 - d$R2) * err_draw(n, d)
  Xs <- scale(X); ys <- drop(scale(y))
  D <- cbind(1, Xs)
  qrD <- qr(D); b <- qr.coef(qrD, ys); res <- ys - D %*% b
  df <- n - d$k - 1; s2 <- sum(res^2) / df
  XtXi <- chol2inv(qr.R(qrD))
  se <- sqrt(diag(XtXi) * s2)[-1]
  r2 <- 1 - sum(res^2) / sum(ys^2)
  Fst <- (r2 / d$k) / ((1 - r2) / df)
  list(adm = 1, est = b[-1], se = se, crit = stats::qt(1 - alpha / 2, df),
       omni = as.numeric(stats::pf(Fst, d$k, df, lower.tail = FALSE) < alpha))
}

rep_anova <- function(d, n, alpha) {
  ng <- floor(n / d$k)
  Y <- matrix(d$means, d$k, ng) + matrix(err_draw(d$k * ng, d), d$k, ng)
  gm <- rowMeans(Y); dfw <- d$k * ng - d$k
  msw <- sum((Y - gm)^2) / dfw
  ssb <- ng * sum((gm - mean(gm))^2)
  Fst <- (ssb / (d$k - 1)) / msw
  list(adm = 1, est = gm[d$k] - gm[1], se = sqrt(2 * msw / ng), crit = stats::qt(1 - alpha / 2, dfw),
       omni = as.numeric(stats::pf(Fst, d$k - 1, dfw, lower.tail = FALSE) < alpha))
}

#' Simulate one design cell
#'
#' Runs \code{reps} Monte Carlo replications of \code{design} at sample size
#' \code{n} and returns the raw per-replication quantities needed to evaluate
#' every SSA criterion for any choice of thresholds.
#'
#' @param design An \code{ssa_design}.
#' @param n Sample size.
#' @param reps Number of replications.
#' @param alpha Two-sided significance level.
#' @param cores Number of cores for \code{\link[parallel]{mclapply}} (forked processes).
#'   Forking is not available on Windows, where the value is reset to 1 with a warning.
#' @return An object of class \code{ssa_cell}.
#' @examples
#' cl <- ssa_cell(ssa_design_anova(k = 2, f = .25), n = 128, reps = 200)
#' summary(cl, delta = .5)
#' @export
ssa_cell <- function(design, n, reps = 500, alpha = .05, cores = 1) {
  check_design(design); check_count(n, "n", design_min_n(design)); check_count(reps, "reps", 1)
  check_prob(alpha, "alpha"); check_count(cores, "cores", 1)
  if (cores > 1 && .Platform$OS.type == "windows") {
    warning("Forked parallel processing is not available on Windows; using cores = 1.")
    cores <- 1
  }
  run <- function(i) ssa_replicate(design, n, alpha)
  out <- if (cores > 1) parallel::mclapply(seq_len(reps), run, mc.cores = cores, mc.set.seed = TRUE)
         else lapply(seq_len(reps), run)
  adm <- vapply(out, `[[`, 0, "adm")
  est <- do.call(rbind, lapply(out, `[[`, "est"))
  se <- do.call(rbind, lapply(out, `[[`, "se"))
  crit <- vapply(out, `[[`, 0, "crit")
  omni <- vapply(out, `[[`, 0, "omni")
  structure(list(design = design, n = n, reps = reps, alpha = alpha, adm = adm, est = est, se = se,
                 crit = crit, omni = omni), class = "ssa_cell")
}

#' Summarize a simulated cell
#'
#' @param object An \code{ssa_cell}.
#' @param delta Precision threshold (target CI width); used for the assurance probability.
#' @param ... Unused.
#' @return Named vector: feasibility (rate of admissible solutions), power (focal, mean over
#'   focal parameters), omnibus power, mean CI width, assurance P(width <= delta),
#'   relative bias, coverage, and Monte Carlo SD-based width.
#' @examples
#' summary(ssa_cell(ssa_design_regression(3, .25), n = 150, reps = 100), delta = .3)
#' @export
summary.ssa_cell <- function(object, delta = NA, ...) {
  ok <- object$adm == 1
  if (!any(ok)) {
    warning("No admissible solutions in this cell; all criteria other than feasibility are NA.")
    return(c(feasibility = 0, power = NA_real_, power_omnibus = NA_real_, width = NA_real_,
             assurance = NA_real_, rel_bias = NA_real_, coverage = NA_real_, width_mc = NA_real_))
  }
  est <- object$est[ok, , drop = FALSE]; se <- object$se[ok, , drop = FALSE]; crit <- object$crit[ok]
  theta <- object$design$theta
  width_rep <- rowMeans(2 * crit * se)
  c(feasibility = mean(object$adm),
    power = mean(abs(est / se) > crit),
    power_omnibus = if (all(is.na(object$omni))) NA else mean(object$omni[ok]),
    width = mean(width_rep),
    assurance = if (is.na(delta)) NA else mean(width_rep <= delta),
    rel_bias = mean(colMeans(est) / theta - 1),
    coverage = mean(abs(sweep(est, 2, theta)) <= crit * se),
    width_mc = mean(2 * stats::qnorm(1 - object$alpha / 2) * apply(est, 2, stats::sd)))
}
