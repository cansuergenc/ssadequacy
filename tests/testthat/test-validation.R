test_that("design constructors reject invalid arguments", {
  expect_error(ssa_design_cfa(3, 5, loadings = c(.4, .5, .6)), "length 1")
  expect_error(ssa_design_cfa(3, 5, loadings = 1.2), "between -1 and 1")
  expect_error(ssa_design_cfa(3, 5, phi = matrix(.3, 2, 2)), "nfactors x nfactors")
  expect_error(ssa_design_cfa(3, 5, estimator = "WLSMV"), "ML")
  expect_error(ssa_design_cfa(0, 5), "nfactors")
  expect_error(ssa_design_regression(3, R2 = 1.5), "R2")
  expect_error(ssa_design_regression(k = 2.5), "'k'")
  expect_error(ssa_design_anova(k = 1), "'k'")
  expect_error(ssa_design_anova(2, f = -.1), "'f'")
  # a full loading vector and a single value are still accepted
  expect_s3_class(ssa_design_cfa(2, 3, loadings = seq(.4, .9, length.out = 6)), "ssa_design")
  expect_equal(ssa_design_cfa(2, 3, loadings = .6)$theta, rep(.6, 6))
})

test_that("feasibility functions require a CFA design", {
  d_reg <- ssa_design_regression(3, .10)
  expect_error(ssa_feasibility(d_reg, 100), "type \"cfa\"")
  expect_error(ssa_n_kappa(d_reg), "type \"cfa\"")
  expect_error(ssa_identification_strength(d_reg, 100), "type \"cfa\"")
  expect_error(ssa_feasibility(list(type = "cfa"), 100), "ssa_design")
  d <- ssa_design_cfa(3, 5, loadings = c(.35, .40))
  expect_error(ssa_feasibility(d, N = -5), "'N'")
  expect_error(ssa_n_kappa(d, kappa = 1), "'kappa'")
  # the calibrated law reproduces the values reported for Study 2
  expect_equal(round(ssa_n_kappa(d, .95)), 152)
  expect_equal(round(ssa_feasibility(d, 100), 2), .79)
})

test_that("ssa_analytic and ssa_cell validate thresholds", {
  d <- ssa_design_anova(2, .25)
  expect_error(ssa_analytic(d, delta = 0), "'delta'")
  expect_error(ssa_analytic(d, delta = .5, power = 1), "'power'")
  expect_error(ssa_analytic(d, delta = .5, alpha = 0), "'alpha'")
  expect_error(ssa_analytic(d, delta = .5, assurance = 1.2), "'assurance'")
  expect_error(ssa_analytic(list(), delta = .5), "ssa_design")
  expect_error(ssa_cell(d, n = 3, reps = 10), "'n'")
  expect_error(ssa_cell(d, n = 50, reps = 0), "'reps'")
})

test_that("summary.ssa_cell warns when no solution is admissible", {
  cl <- ssa_cell(ssa_design_anova(2, .25), n = 20, reps = 3)
  cl$adm[] <- 0
  expect_warning(s <- summary(cl, delta = .5), "No admissible")
  expect_equal(unname(s["feasibility"]), 0)
  expect_true(all(is.na(s[-1])))
})

test_that("ssa_simulate treats precision like the other criteria when met everywhere", {
  set.seed(4)
  expect_message(s <- ssa_simulate(ssa_design_anova(2, .8), delta = 2, n_grid = c(200, 300, 400),
                                   reps = 30, verify = FALSE), "Every criterion is met")
  expect_true(all(is.na(s$N)))
  expect_true(is.na(s$N_star))
  expect_equal(s$binding, "none")
  expect_output(print(s), "no criterion binds")
  expect_match(ssa_report(s, verbose = FALSE), "below this value")
  expect_error(ssa_simulate(ssa_design_anova(2, .8), delta = 2, n_grid = 100, reps = 10), "at least two")
})

test_that("ssa_report returns the text and can be silent", {
  a <- ssa_analytic(ssa_design_anova(2, .25), delta = .5)
  expect_silent(txt <- ssa_report(a, verbose = FALSE))
  expect_type(txt, "character")
  expect_match(txt, "N\\* = 249")
  expect_output(ssa_report(a), "binding criterion")
  expect_error(ssa_report(list()), "ssa_analytic or ssa_sim")
})
