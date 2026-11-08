## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(ssadequacy)

## -----------------------------------------------------------------------------
ssa_crossover()              # alpha = .05, power = .80
ssa_crossover(power = .90)

## -----------------------------------------------------------------------------
d_reg <- ssa_design_regression(k = 3, R2 = .10, rho = .30)
ssa_analytic(d_reg, delta = .25)

## -----------------------------------------------------------------------------
d_aov <- ssa_design_anova(k = 2, f = .25)          # d = .50
ssa_analytic(d_aov, delta = .50, assurance = .80)

## -----------------------------------------------------------------------------
d_cfa <- ssa_design_cfa(nfactors = 3, nind = 5, loadings = c(.35, .40))
ssa_analytic(d_cfa, delta = .20)
ssa_feasibility(d_cfa, N = c(100, 150, 200))

## -----------------------------------------------------------------------------
fit <- lavaan::cfa("visual  =~ x1 + x2 + x3
                    textual =~ x4 + x5 + x6
                    speed   =~ x7 + x8 + x9",
                   data = lavaan::HolzingerSwineford1939, std.lv = TRUE)
d_hs <- ssa_design_from_lavaan(fit)
a_hs <- ssa_analytic(d_hs, delta = .20, assurance = .80, se_cv = 1)
a_hs

## ----eval = FALSE-------------------------------------------------------------
#  s_hs <- ssa_simulate(d_hs, delta = .20, assurance = .80, reps = 1000, cores = 4)
#  s_hs

## -----------------------------------------------------------------------------
ssa_report(a_hs)

