#' Compatibility wrapper for the minimax design problem
#'
#' @description
#' Compatibility wrapper for the stand-alone script name supplied with early
#' versions of the project. This calls [solve_minimax_design()] with
#' `solver = "quadprog"` for the shrinkage formulation and Gurobi only for
#' features that require explicit allocation or direct bounded linear weights.
#'
#' @inheritParams solve_minimax_design
#' @param max_sets Maximum number of feasible experiment sets to enumerate.
#'   Retained for compatibility with the stand-alone implementation. The wrapper
#'   passes this value to [solve_minimax_design()] when supported by the
#'   non-Gurobi backend.
#' @param tol_bisect Numeric tolerance for bisection steps. Retained for
#'   compatibility with the stand-alone implementation.
#' @param bisect_iter Integer maximum number of bisection iterations. Retained
#'   for compatibility with the stand-alone implementation.
#' @param qp_ridge Small nonnegative ridge added to quadratic programs for
#'   numerical stability. Retained for compatibility with the stand-alone
#'   implementation.
#' @param tol_ratio Numeric tolerance used when checking whether normalized
#'   regret ratios are below one because of numerical error.
#' @param verbose Logical; retained for compatibility. The current implementation
#'   is silent apart from warnings.
#'
#' @return An object of class `minimax_design`.
#' @export
solve_minimax_experiment <- function(Sigma_obs,
                                     v2,
                                     n_total,
                                     omega,
                                     c = NULL,
                                     k = NULL,
                                     h = length(omega),
                                     x_max = NULL,
                                     feasible_sets = NULL,
                                     selection_constraints = NULL,
                                     min_experiments = 0L,
                                     n_min = 0,
                                     gamma_lower = NULL,
                                     gamma_upper = NULL,
                                     a_exp_lower = NULL,
                                     a_exp_upper = NULL,
                                     force_gamma_unit_interval = TRUE,
                                     gurobi_params = list(),
                                     max_sets = 200000L,
                                     tol_bisect = 1e-7,
                                     bisect_iter = 80L,
                                     qp_ridge = 1e-10,
                                     zero_tol = 1e-9,
                                     psd_tol = 1e-10,
                                     psd_action = base::c("error", "project"),
                                     tol_ratio = 1e-6,
                                     verbose = FALSE) {
  use_quadprog <- isTRUE(force_gamma_unit_interval) &&
    is.numeric(n_min) && length(n_min) == 1L && is.finite(n_min) && n_min == 0
  backend <- if (use_quadprog) "quadprog" else "gurobi"
  if (isTRUE(verbose)) {
    if (identical(backend, "quadprog")) {
      message("Using non-Gurobi solver: finite-set enumeration + quadprog.")
    } else if (isTRUE(force_gamma_unit_interval)) {
      message("Using Gurobi because n_min > 0 requires explicit allocation.")
    } else {
      message("Using Gurobi for the direct a-weight formulation.")
    }
  }
  solve_minimax_design(
    Sigma_obs = Sigma_obs,
    v2 = v2,
    n_total = n_total,
    omega = omega,
    costs = c,
    bias_weights = k,
    h = h,
    x_max = x_max,
    feasible_sets = feasible_sets,
    selection_constraints = selection_constraints,
    min_experiments = min_experiments,
    n_min = n_min,
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    a_exp_lower = a_exp_lower,
    a_exp_upper = a_exp_upper,
    force_gamma_unit_interval = force_gamma_unit_interval,
    solver = backend,
    gurobi_params = gurobi_params,
    max_sets = max_sets,
    tol_bisect = tol_bisect,
    bisect_iter = bisect_iter,
    qp_ridge = qp_ridge,
    zero_tol = zero_tol,
    psd_tol = psd_tol,
    psd_action = psd_action,
    tol = tol_ratio
  )
}
