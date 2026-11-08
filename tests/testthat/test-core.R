test_that("crossover ratio matches closed form", {
  expect_equal(ssa_crossover(), 2 * qnorm(.975) / (qnorm(.975) + qnorm(.8)))
  expect_equal(round(ssa_crossover(), 3), 1.399)
  expect_lt(ssa_crossover(power = .9), ssa_crossover(power = .8))
})

test_that("Fleishman coefficients reproduce the target moments", {
  cf <- fleishman_coef(2, 7)
  expect_equal(unname(cf), c(-0.2600, 0.7616, 0.2600, 0.0531), tolerance = 1e-3)
  set.seed(1); x <- ssadequacy:::fleishman_transform(rnorm(2e5), cf)
  expect_equal(mean(x^3), 2, tolerance = .15); expect_equal(var(x), 1, tolerance = .02)
})

test_that("non-normal generator reproduces the covariance matrix", {
  S <- matrix(c(1, .4, .2, .4, 1, .3, .2, .3, 1), 3)
  set.seed(2); X <- rmvnonnorm(1e5, S, 2, 7)
  expect_equal(cor(X), S, tolerance = .02, ignore_attr = TRUE)
})

test_that("analytic two-group ANOVA matches classical results", {
  a <- ssa_analytic(ssa_design_anova(2, .25), delta = .5)
  expect_equal(unname(ceiling(a$N["Detectability"])), 128, tolerance = 2)   # d = .5, 64 per group
  expect_equal(a$binding, "Precision")
  b <- ssa_analytic(ssa_design_anova(2, .10), delta = .5)
  expect_equal(b$binding, "Detectability")
})

test_that("binding identity depends only on r (variance invariance)", {
  d1 <- ssa_design_regression(3, .10); d2 <- ssa_design_regression(10, .10)
  for (d in list(d1, d2)) {
    th <- mean(d$theta)
    expect_equal(ssa_analytic(d, delta = 1.2 * th)$binding, "Precision")
    expect_equal(ssa_analytic(d, delta = 1.6 * th)$binding, "Detectability")
  }
})

test_that("simulation cell returns sensible criteria", {
  set.seed(3)
  s <- summary(ssa_cell(ssa_design_regression(3, .25), 200, reps = 200), delta = .3)
  expect_equal(unname(s["feasibility"]), 1)
  expect_gt(s["coverage"], .9)
  expect_lt(abs(s["rel_bias"]), .05)
})
