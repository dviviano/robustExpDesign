#' Solve the minimax design problem without Gurobi
#'
#' @description
#' Compatibility wrapper for the stand-alone script name supplied with early
#' versions of the project. This calls [solve_minimax_design()] with
#' `solver = "quadprog"`.
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
                                     kappa = NULL,
                                     h = length(omega),
                                     x_max = NULL,
                                     feasible_sets = NULL,
                                     selection_constraints = NULL,
                                     min_experiments = 0L,
                                     gamma_lower = NULL,
                                     gamma_upper = NULL,
                                     max_sets = 200000L,
                                     tol_bisect = 1e-7,
                                     bisect_iter = 80L,
                                     qp_ridge = 1e-10,
                                     tol_ratio = 1e-6,
                                     verbose = FALSE) {
  if (isTRUE(verbose)) {
    message("Using non-Gurobi solver: finite-set enumeration + quadprog.")
  }
  solve_minimax_design(
    Sigma_obs = Sigma_obs,
    v2 = v2,
    n_total = n_total,
    omega = omega,
    costs = c,
    bias_weights = k,
    kappa = kappa,
    h = h,
    x_max = x_max,
    feasible_sets = feasible_sets,
    selection_constraints = selection_constraints,
    min_experiments = min_experiments,
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    solver = "quadprog",
    max_sets = max_sets,
    tol_bisect = tol_bisect,
    bisect_iter = bisect_iter,
    qp_ridge = qp_ridge,
    tol = tol_ratio
  )
}
