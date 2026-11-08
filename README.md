# ssadequacy

Sample Size Adequacy (SSA): joint power, precision, stability and feasibility
planning for confirmatory factor analysis, multiple regression and one-way
ANOVA.

A study is *adequate* when four requirements hold at the same time for the
parameters that carry its claims (the focal parameters):

| Criterion | Requirement |
|---|---|
| Detectability | power of the focal Wald tests is at least π |
| Precision | the 95% confidence interval width is at most δ (in expectation, or with assurance γ) |
| Stability | the absolute relative bias is at most τ |
| Feasibility | the rate of admissible (converged, Heywood-free) solutions is at least κ |

The minimum adequate sample size N\* is the largest of the four
criterion-specific requirements, and the criterion that attains it is the
*bottleneck*. `ssadequacy` predicts the bottleneck in closed form and verifies
it by Monte Carlo simulation.

## Installation

```r
# from CRAN (once accepted)
install.packages("ssadequacy")

# development version
# install.packages("remotes")
remotes::install_github("cansuergenc/ssadequacy")
```

## Workflow

```r
library(ssadequacy)

# 1. Describe the design (population model, focal parameters, data distribution)
d <- ssa_design_cfa(nfactors = 3, nind = 5, loadings = c(.35, .40))

# 2. Closed-form bottleneck analysis
a <- ssa_analytic(d, delta = .20, power = .80)
a

# Which criterion binds between power and precision depends only on
# r = delta / |theta| relative to the crossover ratio r*
ssa_crossover()            # 1.399 for alpha = .05, power = .80

# Predicted admissibility of CFA solutions (calibrated scaling law)
ssa_feasibility(d, N = c(100, 150, 200))

# 3. Monte Carlo verification with root finding (slower)
s <- ssa_simulate(d, delta = .20, reps = 500)
s

# 4. Justification paragraph for a preregistration or Method section
ssa_report(a)

# Interactive planner
run_ssa_app()
```

The planning population can also be taken from a fitted `lavaan` model, for
example a pilot study:

```r
fit <- lavaan::cfa("visual  =~ x1 + x2 + x3
                    textual =~ x4 + x5 + x6
                    speed   =~ x7 + x8 + x9",
                   data = lavaan::HolzingerSwineford1939, std.lv = TRUE)
ssa_analytic(ssa_design_from_lavaan(fit), delta = .20, assurance = .80, se_cv = 1)
```

See `vignette("ssadequacy")` for a full walk-through.

## Citation

Ergenç Özdaş, C. (2026). *Beyond rules of thumb: A bottleneck theory of sample
size adequacy.* Manuscript in preparation. See `citation("ssadequacy")`.

## License

GPL (>= 3)
