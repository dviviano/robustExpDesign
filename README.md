# robustExpDesign

Robust experimental design for combining experimental and external evidence.

`robustExpDesign` implements tools for planning experiments when the researcher wants to combine experimental evidence with observational or external evidence that may be biased. The package chooses which experiments or moments to use, how to allocate the experimental budget, and how to combine the resulting estimates.

The package is designed for settings with

- a target estimand that depends on an underlying parameter vector;
- observational or external estimates that may be misspecified;
- a menu of feasible experiments, treatment arms, sites, moments, or data sources;
- a budget or feasibility constraint on which experiments can be run;
- design-stage covariance and cost information for power calculations.

The typical workflow is:

1. specify the target sensitivity vector `omega`, or a target-sensitivity matrix `Omega`;
2. provide observational covariance information and experimental variance inputs;
3. specify costs, experiment-count limits, and feasibility restrictions;
4. solve the minimax proportional-regret design problem;
5. inspect selected experiments, estimator weights, sample allocation, and regret components.

## Reference

This package implements methods from:

Epanomeritakis, A., & Viviano, D. (2026). *Learning What to Learn: Experimental Design when Combining Experimental with Observational Evidence*. Manuscript, July 9, 2026.

## Installation

Install the development version from GitHub with:

```r
install.packages("remotes")
remotes::install_github("ostasovskyi/robustExpDesign")
```

Load the package with:

```r
library(robustExpDesign)
```

Core utilities use `Matrix`, `rlang`, and `utils`.

```r
install.packages(c("Matrix", "rlang"))
```

The non-Gurobi baseline solver uses `quadprog`.

```r
install.packages("quadprog")
```

The moment-selection and audience-regret solvers use `CVXR`.

```r
install.packages("CVXR")
```

Gurobi is optional. To use the mixed-integer backend, install the Gurobi Optimizer and the Gurobi R package. The Gurobi R package is not installed from CRAN; it is installed from the Gurobi installation directory after Gurobi has been installed and licensed.

On Linux or WSL, the command usually has the following form, with the path adjusted to your installed Gurobi version:

```sh
R CMD INSTALL /opt/gurobi/linux64/R/gurobi_*.tar.gz
```

Check that the R interface is available:

```r
library(gurobi)
```

## Quick start

The following example solves a small three-parameter design problem. The researcher can select at most two experiments and allocate a total sample size across selected experimental arms. The example uses the non-Gurobi backend.

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
fit$x_opt
fit$gamma_opt
fit$n_opt
fit$alpha_ratio
fit$beta_ratio
fit$regret
```

The returned object contains the selected experiments, shrinkage weights, sample allocation, oracle normalizations, and the two components of the minimax proportional-regret objective.

## Objective components

For a candidate design, the package evaluates

```text
max(variance regret, bias regret).
```

The variance-regret component compares the design's variance index with the smallest feasible variance index. The bias-regret component compares the design's worst-case bias exposure with the smallest feasible bias exposure. The minimax design minimizes the larger of these two ratios.

```r
fit$alpha_opt
fit$beta_opt
fit$alpha_star
fit$beta_star
fit$alpha_ratio
fit$beta_ratio
fit$regret
```

## Weighted bias radii

The package supports weighted ambiguity sets for observational bias. Passing `bias_weights` implements coordinate-specific bounds of the form

```text
|b_j| <= bias_weights[j] * B.
```

The aliases `kappa` and `k` are also accepted.

```r
fit_weighted <- solve_minimax_design(
  Sigma_obs = Sigma_obs,
  v2 = v2,
  n_total = 500,
  omega = omega,
  costs = c(1, 1, 1.5),
  kappa = c(1.0, 0.5, 2.0),
  h = 2,
  min_experiments = 1,
  solver = "quadprog"
)

fit_weighted$beta_opt
fit_weighted$regret
```

Weights equal to zero can be used for coordinates treated as not misspecified.

## Feasibility restrictions

Experiment menus can be restricted with experiment-count bounds, eligibility indicators, explicit feasible sets, or linear constraints on the selected experiments.

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
  x_max = c(1, 1, 1),
  selection_constraints = selection_constraints,
  min_experiments = 1,
  solver = "quadprog"
)

fit_restricted$selected
fit_restricted$n_opt
```

## Budget sweeps

`sweep_n_total()` runs the design solvers over a grid of total budgets and experiment-count limits. It returns arm-level allocation rows, regret rows, and optional plots.

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

The moment-selection routines are for settings where the researcher wants to combine several moments, instruments, sites, outcomes, or data sources rather than choosing experiments that each estimate one parameter directly. Some moments may be reliable, while others may be useful but potentially biased.

The main function is `solve_moment_design()`. The user provides `Lambda`, the moment Jacobian; `Sigma`, the covariance matrix of the moments; and `omega` or `Omega`, which describes how the target estimand depends on the structural parameters. The user also specifies which moments may be biased and either supplies a fixed GMM weighting matrix `W` or asks the package to search over allowed moment combinations using masks.

### Fixed weighting matrix

Use `solve_moment_design(..., optimize_W = FALSE)` when the researcher wants to evaluate a user-supplied GMM weighting matrix.

```r
Lambda <- matrix(
  c(
    1.0, 0.0,
    0.0, 1.0,
    1.0, 0.5,
    0.5, 1.0
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(
    c("exp_direct", "exp_income", "obs_joint", "obs_aux"),
    c("theta_direct", "theta_income")
  )
)

Sigma_mom <- diag(c(0.05, 0.06, 0.10, 0.12))
Omega <- matrix(c(0.4, 1.0), nrow = 1)
W <- diag(4)

eval_fixed <- solve_moment_design(
  Lambda = Lambda,
  Sigma = Sigma_mom,
  Omega = Omega,
  W = W,
  biased_moments = c(3, 4),
  norm = "l2",
  optimize_W = FALSE
)

eval_fixed$alpha
eval_fixed$beta
```

### Optimizing over moment sets

Use `solve_moment_design()` with `optimize_W = TRUE` to optimize over the first-order linear estimator implied by selected moments. `W_mask`, `W_masks`, or `moment_sets` specify which moments are available. Moments excluded by the mask are assigned zero weight in the optimized linear estimator.

```r
fit_moments <- solve_moment_design(
  Lambda = Lambda,
  Sigma = Sigma_mom,
  Omega = Omega,
  biased_moments = c(3, 4),
  bias_weights = c(0, 0, 1, 1.5),
  norm = "linf",
  moment_sets = list(c(1, 2, 3), c(1, 2, 4), c(1, 2, 3, 4)),
  optimize_W = TRUE,
  cvxr_solver = "CLARABEL"
)

fit_moments$selected_candidate
fit_moments$alpha_ratio
fit_moments$beta_ratio
fit_moments$regret
```

## Audience-regret designs

`solve_audience_design()` solves the audience-regret objective over a supplied grid of bias radii `B_grid`. The function computes oracle benchmarks at each grid value and selects one experiment set and one sample allocation that control regret over the grid.

```r
fit_audience <- solve_audience_design(
  Sigma_obs = Sigma_obs,
  v2 = v2,
  n_total = 500,
  omega = omega,
  B_grid = c(0, 0.25, 0.5, 1, 2),
  costs = c(1, 1, 1.5),
  kappa = c(1, 1, 1),
  h = 2,
  min_experiments = 1,
  cvxr_solver = "CLARABEL"
)

fit_audience
fit_audience$n_opt
fit_audience$posterior_ratio
fit_audience$beta_ratio_large_B
```

For general linear reported statistics, pass `Lambda_obs` and `Lambda_exp`. The function then imposes the calibration constraint implied by those loading matrices.

## Core functions

### `solve_minimax_design()`

Main experiment-selection function. It chooses experiment indicators, shrinkage weights, and cost-aware sample allocation to minimize proportional regret.

### `solve_variance_design()`

Variance-optimal benchmark under the same feasibility restrictions.

### `evaluate_design()`

Evaluates a user-supplied design and returns sample allocation, variance, and bias diagnostics.

### `sweep_n_total()`

Runs the design solvers over budget and experiment-count grids and collects allocation, shrinkage, and regret diagnostics.

### `solve_moment_design()`

Solves fixed-weight and optimized moment-design problems using the first-order linear representation of GMM estimators.

### `solve_audience_design()`

Solves the audience-regret design problem over a user-supplied `B_grid`.
