# robustExpDesign

Robust experimental design for combining experimental and external evidence.

`robustExpDesign` provides tools for planning experiments when observational or
other external estimates may be biased. It can choose experiments or moments,
allocate a cost-constrained experimental budget, and compute estimator weights
that balance variance against sensitivity to misspecification.

The package implements methods from:

> Epanomeritakis, A., and Viviano, D. (2026). *Learning What to Learn:
> Experimental Design when Combining Experimental with Observational Evidence*.
> Manuscript, August 4, 2026.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("dviviano/robustExpDesign")
```

`Matrix` and `rlang` are required dependencies and are installed with the
package. Install the optional packages needed by the features you plan to use:

```r
install.packages(c("quadprog", "CVXR", "ggplot2"))
```

The optional Gurobi R interface is required for:

- audience-regret designs;
- designs with `n_min > 0`;
- the direct experimental-weight formulation obtained with
  `force_gamma_unit_interval = FALSE`;
- large mixed-integer problems for which enumerating feasible experiment sets
  is impractical.

Install the Gurobi optimizer and its R package using the instructions for your
local Gurobi installation.

## Quick start

The following example uses the non-Gurobi backend. It selects at most two of
three experiments and allocates a total cost budget of 500.

```r
library(robustExpDesign)

Sigma_obs <- matrix(
  c(
    0.08, 0.01, 0.00,
    0.01, 0.05, 0.01,
    0.00, 0.01, 0.09
  ),
  nrow = 3,
  byrow = TRUE
)

v2 <- c(direct = 1.00, income = 1.40, wage = 1.10)
omega <- c(direct = 0.40, income = 1.00, wage = -0.70)

fit <- solve_minimax_design(
  Sigma_obs = Sigma_obs,
  v2 = v2,
  n_total = 500,
  omega = omega,
  costs = c(1, 1, 1.5),
  bias_weights = c(1, 1, 1),
  h = 2,
  min_experiments = 1,
  solver = "quadprog"
)

fit
fit$selected
fit$gamma_opt
fit$n_opt
fit$alpha_ratio
fit$beta_ratio
fit$regret
```

The proportional-regret objective is

```text
max(variance regret, bias regret).
```

The result includes the selected experiments, effective shrinkage weights,
cost-aware allocation, oracle normalizations, variance and bias ratios, and the
maximum of those two ratios.

## Evaluate a supplied design

Use `evaluate_design()` when the selection and shrinkage weights are already
specified:

```r
evaluated <- evaluate_design(
  Sigma_obs = Sigma_obs,
  v2 = v2,
  n_total = 500,
  omega = omega,
  x = c(1, 1, 0),
  gamma = c(0.8, 0.6, 0),
  costs = c(1, 1, 1.5)
)

evaluated$n_opt
evaluated$alpha
evaluated$beta
```

## Weighted bias radii

`bias_weights` implements coordinate-specific ambiguity bounds of the form

```text
|b_j| <= bias_weights[j] * B.
```

A weight of zero treats the corresponding coordinate as not misspecified. The
legacy alias `k` remains available for compatibility.

```r
fit_weighted <- solve_minimax_design(
  Sigma_obs = Sigma_obs,
  v2 = v2,
  n_total = 500,
  omega = omega,
  costs = c(1, 1, 1.5),
  bias_weights = c(1.0, 0.5, 2.0),
  h = 2,
  min_experiments = 1,
  solver = "quadprog"
)
```

## Feasibility restrictions

Experiment menus can be restricted with `h`, `min_experiments`, `x_max`,
`feasible_sets`, or linear constraints on the binary selection vector.

```r
selection_constraints <- list(
  A = matrix(c(1, 1, 0), nrow = 1),
  sense = "<=",
  rhs = 1
)

fit_restricted <- solve_minimax_design(
  Sigma_obs = Sigma_obs,
  v2 = v2,
  n_total = 500,
  omega = omega,
  costs = c(1, 1, 1.5),
  h = 2,
  selection_constraints = selection_constraints,
  min_experiments = 1,
  solver = "quadprog"
)
```

Use `n_min` to impose a minimum sample allocation on every selected experiment.
Positive `n_min` currently requires Gurobi, and the same floor is imposed on
candidate designs and oracle benchmarks.

## Budget sweeps

`sweep_n_total()` runs the design solvers over budgets and experiment-count
limits. It returns arm-level allocation rows, regret rows, fitted objects, and
optional `ggplot2` plots.

```r
sweep <- sweep_n_total(
  Sigma_obs = Sigma_obs,
  v2 = v2,
  omega = omega,
  n_grid = c(200, 300, 400, 500, 750, 1000),
  h_values = c(1, 2),
  costs = c(1, 1, 1.5),
  experiment_names = c("Direct", "Income", "Wage"),
  include_variance_design = TRUE,
  make_plots = TRUE,
  solver = "quadprog"
)

head(sweep$data$arms)
head(sweep$data$regret)
sweep$plots$allocation
sweep$plots$regret
```

## Moment selection and GMM weights

`solve_moment_design()` handles general moment-loading problems. With a supplied
weighting matrix and `optimize_W = FALSE`, it evaluates the implied linear GMM
estimator. With `optimize_W = TRUE`, it uses `CVXR` to optimize over admissible
moment sets or weighting-matrix masks.

```r
Lambda <- matrix(
  c(
    1.0, 0.0,
    0.0, 1.0,
    1.0, 0.5,
    0.5, 1.0
  ),
  nrow = 4,
  byrow = TRUE
)

Sigma_mom <- diag(c(0.05, 0.06, 0.10, 0.12))
Omega <- matrix(c(0.4, 1.0), nrow = 1)

fixed_moments <- solve_moment_design(
  Lambda = Lambda,
  Sigma = Sigma_mom,
  Omega = Omega,
  W = diag(4),
  biased_moments = c(3, 4),
  norm = "l2",
  optimize_W = FALSE
)

fixed_moments$alpha
fixed_moments$beta
```

## Audience-regret designs

`solve_audience_regret_design()` solves the Appendix C.2 objective over a grid
of scalarization weights. The transformation is

```text
lambda = B^2 / (1 + B^2),
```

and the grid must include 0 and 1. The default uses 50 equally spaced points.
This solver currently requires Gurobi.

```r
fit_audience <- solve_audience_regret_design(
  Sigma_obs = Sigma_obs,
  v2 = v2,
  n_total = 500,
  omega = omega,
  costs = c(1, 1, 1.5),
  bias_weights = c(1, 1, 1),
  h = 2,
  min_experiments = 1
)

fit_audience$r_opt
fit_audience$risk_by_lambda
```