# ssadequacy 0.1.0

Initial CRAN submission.

* Design constructors for confirmatory factor analysis (`ssa_design_cfa()`,
  `ssa_design_from_lavaan()`), multiple regression (`ssa_design_regression()`)
  and one-way ANOVA (`ssa_design_anova()`), with normal or non-normal
  (Vale-Maurelli) data generation.
* `ssa_analytic()`: closed-form sample size requirements for detectability,
  precision (optionally with assurance), stability and feasibility, with the
  power-precision crossover ratio `ssa_crossover()`.
* `ssa_feasibility()`, `ssa_n_kappa()` and `ssa_identification_strength()`:
  calibrated weak-identification scaling law for the admissibility of CFA
  solutions.
* `ssa_cell()` and `ssa_simulate()`: Monte Carlo verification with
  theory-based smoothing and root finding.
* `ssa_report()`: justification paragraph for preregistrations and Method
  sections; `run_ssa_app()`: interactive 'shiny' planner.
* All user-facing functions validate their arguments and give informative
  errors; `ssa_simulate()` reports "met at every grid point" consistently for
  all four criteria and handles the case in which no criterion binds.
