#' Solve the minimax proportional-regret design problem
#'
#' @description
#' `solve_minimax_design()` implements the baseline optimization routine in
#' Algorithm 1 of Epanomeritakis and Viviano (2025). It chooses experiment
#' indicators `x`, shrinkage weights `gamma`, and cost-aware sample allocation
#' `n_opt` that minimize the maximum of the normalized variance and bias
#' components.
#'
#' The default `solver = "auto"` uses Gurobi when available and otherwise falls
#' back to an exact finite-set enumeration plus convex quadratic programs solved
#' by `quadprog`. The non-Gurobi solver is intended for moderate numbers of
#' feasible experiment sets; for very large menus, use Gurobi or pass
#' `feasible_sets` explicitly.
#'
#' @param Sigma_obs Numeric `p` by `p` covariance matrix for the observational
#'   estimators.
#' @param v2 Numeric length-`p` vector. The experimental estimator for coordinate
#'   `j` has variance `v2[j] / n_j`.
#' @param n_total Positive scalar budget `n`. If all entries of `costs` are one,
#'   this is the total experimental sample size. Otherwise it is the budget in
#'   the constraint `sum(costs * n_opt) = n_total`.
#' @param omega Numeric length-`p` sensitivity vector. For a nonlinear estimand,
#'   use the gradient evaluated at the design-stage estimate.
#' @param costs Optional positive length-`p` vector of per-unit costs. Defaults to
#'   one for each coordinate.
#' @param bias_weights Optional nonnegative length-`p` vector. Defaults to one for
#'   each coordinate. Values other than one implement the weighted ambiguity set
#'   `|b_j| <= bias_weights[j] * B`.
#' @param h Integer upper bound on the number of selected experiments. Defaults
#'   to `length(omega)`.
#' @param x_max Optional length-`p` vector of zeros and ones indicating which
#'   experiments are eligible. Defaults to one for each coordinate.
#' @param feasible_sets Optional list of integer vectors. Each vector gives a
#'   feasible selected set using 1-based indices. When supplied, the solver
#'   selects one member of this list after applying `h`, `x_max`,
#'   `min_experiments`, and `selection_constraints`.
#' @param selection_constraints Optional list with entries `A`, `sense`, and
#'   `rhs`, defining additional linear constraints on `x`. Senses may be `"<="`,
#'   `">="`, `"<"`, `">"`, or `"="`.
#' @param min_experiments Integer lower bound on the number of selected
#'   experiments. Defaults to zero.
#' @param gamma_lower,gamma_upper Optional length-`p` vectors giving coordinate
#'   bounds on `gamma`. Defaults are zero and one. Set both to one to use only
#'   experimental estimates when selected.
#' @param solver One of `"auto"`, `"gurobi"`, or `"quadprog"`. `"quadprog"` avoids
#'   Gurobi by enumerating feasible experiment sets and solving the convex
#'   continuous subproblems with `quadprog`.
#' @param gurobi_params Optional named list of parameters passed to the Gurobi R
#'   optimizer when `solver = "gurobi"`.
#' @param max_sets Maximum number of feasible experiment sets to enumerate when
#'   using the non-Gurobi solver.
#' @param tol_bisect,bisect_iter Controls for the bisection used by the non-Gurobi
#'   minimax solver.
#' @param qp_ridge Small ridge added to quadratic-program Hessians for numerical
#'   positive-definiteness in `quadprog`.
#' @param tol Numeric tolerance used for diagnostic warnings.
#' @param c Alias for `costs`, included for script compatibility.
#' @param k,kappa Aliases for `bias_weights`, included for script compatibility.
#'
#' @return An object of class `minimax_design`.
#' @export
solve_minimax_design <- function(Sigma_obs,
                                 v2,
                                 n_total,
                                 omega,
                                 costs = NULL,
                                 bias_weights = NULL,
                                 h = length(omega),
                                 x_max = NULL,
                                 feasible_sets = NULL,
                                 selection_constraints = NULL,
                                 min_experiments = 0L,
                                 gamma_lower = NULL,
                                 gamma_upper = NULL,
                                 solver = c("auto", "gurobi", "quadprog"),
                                 gurobi_params = list(),
                                 max_sets = 200000L,
                                 tol_bisect = 1e-7,
                                 bisect_iter = 80L,
                                 qp_ridge = 1e-10,
                                 tol = 1e-6,
                                 c = NULL,
                                 k = NULL,
                                 kappa = NULL) {
  prep <- .prepare_minimax_problem(
    Sigma_obs = Sigma_obs,
    v2 = v2,
    n_total = n_total,
    omega = omega,
    costs = costs,
    bias_weights = bias_weights,
    h = h,
    x_max = x_max,
    feasible_sets = feasible_sets,
    selection_constraints = selection_constraints,
    min_experiments = min_experiments,
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    costs_alias = c,
    bias_weights_alias = k,
    kappa_alias = kappa
  )

  solver <- .choose_baseline_solver(solver)

  if (identical(solver, "gurobi")) {
    alpha_oracle <- .solve_alpha_oracle(prep, gurobi_params)
    beta_oracle <- .solve_beta_oracle(prep, gurobi_params)
  } else {
    alpha_oracle <- .solve_alpha_oracle_enum(prep, max_sets = max_sets,
                                             qp_ridge = qp_ridge)
    beta_oracle <- .solve_beta_oracle_enum(prep, max_sets = max_sets)
  }
  alpha_star <- alpha_oracle$alpha
  beta_star <- beta_oracle$beta

  if (!is.finite(alpha_star) || alpha_star <= 0) {
    stop("The variance oracle returned a nonpositive value. Check Sigma_obs and v2.",
         call. = FALSE)
  }

  minimax <- if (identical(solver, "gurobi")) {
    .solve_minimax_miqcp(prep, alpha_star, beta_star, gurobi_params)
  } else {
    .solve_minimax_enum(
      prep, alpha_star, beta_star, max_sets = max_sets,
      tol_bisect = tol_bisect, bisect_iter = bisect_iter, qp_ridge = qp_ridge
    )
  }

  if (minimax$alpha_ratio < 1 - tol) {
    warning(sprintf("alpha_ratio = %.8f is below one by more than tol.", minimax$alpha_ratio),
            call. = FALSE)
  }
  if (is.finite(minimax$beta_ratio) && minimax$beta_ratio < 1 - tol) {
    warning(sprintf("beta_ratio = %.8f is below one by more than tol.", minimax$beta_ratio),
            call. = FALSE)
  }

  n_opt <- prep$compute_allocation(minimax$s)
  names_out <- prep$arm_names

  out <- list(
    selected = prep$arm_names[which(minimax$x == 1L)],
    x_opt = .make_named(minimax$x, names_out),
    gamma_opt = .make_named(minimax$gamma, names_out),
    s_opt = .make_named(minimax$s, names_out),
    n_opt = .make_named(n_opt, names_out),
    alpha_opt = minimax$alpha,
    beta_opt = minimax$beta,
    alpha_star = alpha_star,
    beta_star = beta_star,
    alpha_ratio = minimax$alpha_ratio,
    beta_ratio = minimax$beta_ratio,
    regret = minimax$regret,
    t_opt = minimax$t,
    solver = solver,
    alpha_oracle = list(
      x = .make_named(alpha_oracle$x, names_out),
      gamma = .make_named(alpha_oracle$gamma, names_out),
      s = .make_named(alpha_oracle$s, names_out),
      n = .make_named(prep$compute_allocation(alpha_oracle$s), names_out),
      alpha = alpha_oracle$alpha,
      beta = alpha_oracle$beta
    ),
    beta_oracle = list(
      x = .make_named(beta_oracle$x, names_out),
      gamma = .make_named(beta_oracle$gamma, names_out),
      s = .make_named(beta_oracle$s, names_out),
      n = .make_named(prep$compute_allocation(beta_oracle$s), names_out),
      alpha = beta_oracle$alpha,
      beta = beta_oracle$beta
    ),
    feasible_sets_used = if (prep$use_sets) prep$sets else NULL,
    optimizer_result = minimax$result,
    gurobi_result = if (identical(solver, "gurobi")) minimax$result else NULL,
    inputs = list(
      Sigma_obs = prep$Sigma_obs,
      v2 = prep$v2,
      n_total = prep$n_total,
      omega = .make_named(prep$omega, names_out),
      costs = .make_named(prep$costs, names_out),
      bias_weights = .make_named(prep$bias_weights, names_out),
      h = prep$h,
      min_experiments = prep$min_experiments,
      gamma_lower = .make_named(prep$gamma_lower, names_out),
      gamma_upper = .make_named(prep$gamma_upper, names_out)
    )
  )
  class(out) <- base::c("minimax_design", "list")
  out
}

#' @export
print.minimax_design <- function(x, ...) {
  cat("Regret-optimal design\n")
  cat("Selected experiments: ", .format_selected(x$x_opt, names(x$x_opt)), "\n", sep = "")
  cat("Solver: ", if (is.null(x$solver)) "unknown" else x$solver, "\n", sep = "")
  cat("Regret: ", format(x$regret, digits = 6), "\n", sep = "")
  cat("Variance ratio: ", format(x$alpha_ratio, digits = 6), "\n", sep = "")
  cat("Bias ratio: ", format(x$beta_ratio, digits = 6), "\n", sep = "")
  invisible(x)
}
