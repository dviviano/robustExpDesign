#' robustExpDesign: Robust Experimental Design with Observational Evidence
#'
#' Tools for choosing experiments, allocating experimental budgets, and
#' combining experimental estimates with observational or external estimates
#' that may be misspecified.
#'
#' @section Main functions:
#' * [solve_minimax_design()] solves the baseline proportional-regret design.
#' * [solve_variance_design()] computes the variance-optimal benchmark.
#' * [evaluate_design()] evaluates a user-supplied design.
#' * [sweep_n_total()] runs budget and experiment-count sweeps.
#' * [solve_moment_design()] handles general moment and GMM designs.
#' * [solve_audience_regret_design()] solves the finite-grid audience-regret
#'   design.
#'
#' @references
#' Epanomeritakis, A., and Viviano, D. (2026). *Learning What to Learn:
#' Experimental Design when Combining Experimental with Observational
#' Evidence*. Manuscript, August 4, 2026.
#'
#' @keywords internal
"_PACKAGE"
