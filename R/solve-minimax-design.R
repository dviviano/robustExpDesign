#' Solve the minimax proportional-regret design problem
#'
#' @description
#' `solve_minimax_design()` implements the baseline optimization routine in
#' Algorithm 1 of Epanomeritakis and Viviano (2026) for the identity-loading
#' specialization, where each experimental estimate targets the corresponding
#' coordinate observed externally. It chooses experiment
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
#' @param n_min Optional nonnegative lower bound on the sample size allocated to
#'   each selected experiment. Defaults to zero, which preserves the profiled
#'   Neyman-allocation implementation.
#' @param gamma_lower,gamma_upper Optional length-`p` vectors giving coordinate
#'   bounds on `gamma`. Defaults are zero and one. Set both to one to use only
#'   experimental estimates when selected. When `force_gamma_unit_interval =
#'   FALSE` and `a_exp_lower`/`a_exp_upper` are not supplied, these bounds are
#'   translated into direct experimental weight bounds through `a_exp = omega *
#'   gamma`.
#' @param a_exp_lower,a_exp_upper Optional length-`p` vectors giving direct lower
#'   and upper bounds for experimental linear weights in the general
#'   identity-Lambda formulation used when `force_gamma_unit_interval = FALSE`.
#' @param force_gamma_unit_interval Logical. If `TRUE`, require
#'   `0 <= gamma_lower <= gamma_upper <= 1`, preserving the original shrinkage
#'   interpretation. If `FALSE`, solve the general identity-Lambda linear-weight
#'   formulation over `a_exp`; this currently requires Gurobi.
#' @param solver One of `"auto"`, `"gurobi"`, or `"quadprog"`. `"quadprog"` avoids
#'   Gurobi by enumerating feasible experiment sets and solving the convex
#'   continuous subproblems with `quadprog`.
#' @param gurobi_params Optional named list of parameters passed to the Gurobi R
#'   optimizer when `solver = "gurobi"`.
#' @param max_sets Maximum number of feasible experiment sets to enumerate when
#'   using the non-Gurobi solver.
#' @param tol_bisect,bisect_iter Controls for the bisection used by the non-Gurobi
#'   minimax solver.
#' @param qp_ridge Relative, dimensionless ridge added after normalizing the
#'   quadratic-program objective. This preserves target-unit invariance while
#'   ensuring positive-definiteness for `quadprog`.
#' @param zero_tol Relative tolerance used only to recognize a mathematically
#'   zero bias oracle. It is applied against a homogeneous bias scale.
#' @param psd_tol Relative eigenvalue tolerance for covariance validation.
#' @param psd_action Either `"error"` or `"project"`. Corrections beyond
#'   floating-point roundoff require explicit `"project"`; every projection is
#'   reported in the fitted object.
#' @param tol Numeric tolerance used for diagnostic warnings.
#' @param c Alias for `costs`, included for script compatibility.
#' @param k Alias for `bias_weights`, included for script compatibility.
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
                                 n_min = 0,
                                 gamma_lower = NULL,
                                 gamma_upper = NULL,
                                 a_exp_lower = NULL,
                                 a_exp_upper = NULL,
                                 force_gamma_unit_interval = TRUE,
                                 solver = base::c("auto", "gurobi", "quadprog"),
                                 gurobi_params = list(),
                                 max_sets = 200000L,
                                 tol_bisect = 1e-7,
                                 bisect_iter = 80L,
                                 qp_ridge = 1e-10,
                                 zero_tol = 1e-9,
                                 psd_tol = 1e-10,
                                 psd_action = base::c("error", "project"),
                                 tol = 1e-6,
                                 c = NULL,
                                 k = NULL) {
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
    n_min = n_min,
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    a_exp_lower = a_exp_lower,
    a_exp_upper = a_exp_upper,
    force_gamma_unit_interval = force_gamma_unit_interval,
    costs_alias = c,
    bias_weights_alias = k,
    zero_tol = zero_tol,
    psd_tol = psd_tol,
    psd_action = psd_action
  )

  solver <- .choose_baseline_solver(solver)
  if (identical(prep$weight_mode, "a") && !identical(solver, "gurobi")) {
    stop("The direct a-weight formulation with force_gamma_unit_interval = FALSE currently requires solver = 'gurobi'.",
         call. = FALSE)
  }
  if (prep$n_min > 0 && !identical(solver, "gurobi")) {
    stop("n_min > 0 currently requires solver = 'gurobi'.", call. = FALSE)
  }

  if (prep$n_min > 0) {
    alpha_oracle <- .solve_alpha_oracle_explicit_alloc(prep, gurobi_params)
    beta_oracle <- .solve_beta_oracle_explicit_alloc(prep, gurobi_params)
  } else if (identical(prep$weight_mode, "a")) {
    alpha_oracle <- if (prep$use_sets) {
      .solve_alpha_oracle_a_sets(prep, gurobi_params)
    } else {
      .solve_alpha_oracle_a(prep, gurobi_params)
    }
    beta_oracle <- .solve_beta_oracle_a(prep, gurobi_params)
  } else if (identical(solver, "gurobi")) {
    alpha_oracle <- .solve_alpha_oracle(prep, gurobi_params)
    beta_oracle <- .solve_beta_oracle(prep, gurobi_params)
  } else {
    alpha_oracle <- .solve_alpha_oracle_enum(prep, max_sets = max_sets,
                                             qp_ridge = qp_ridge)
    beta_oracle <- .solve_beta_oracle_enum(prep, max_sets = max_sets)
  }
  alpha_star <- alpha_oracle$alpha
  beta_star <- beta_oracle$beta
  beta_zero <- .beta_zero_feasible(prep, beta_oracle$x)
  prep$beta_zero <- beta_zero
  if (beta_zero) {
    beta_star <- 0
    beta_oracle$beta <- 0
  }

  if (!is.finite(alpha_star) || alpha_star <= 0) {
    stop("The variance oracle returned a nonpositive value. Check Sigma_obs and v2.",
         call. = FALSE)
  }

  minimax <- if (prep$n_min > 0) {
    .solve_minimax_explicit_alloc(prep, alpha_star, beta_star, gurobi_params)
  } else if (identical(prep$weight_mode, "a")) {
    if (prep$use_sets) {
      .solve_minimax_a_sets(prep, alpha_star, beta_star, gurobi_params)
    } else {
      .solve_minimax_a_miqcp(prep, alpha_star, beta_star, gurobi_params)
    }
  } else if (identical(solver, "gurobi")) {
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

  solution_weight <- if (identical(prep$weight_mode, "a")) minimax$a_exp else minimax$s
  n_opt <- if (!is.null(minimax$n)) {
    minimax$n
  } else {
    prep$compute_allocation(solution_weight, x = minimax$x)
  }
  names_out <- prep$arm_names
  minimax_a_exp <- if (is.null(minimax$a_exp)) prep$omega * minimax$s else minimax$a_exp
  minimax_a_obs <- if (is.null(minimax$a_obs)) prep$omega * (1 - minimax$s) else minimax$a_obs
  alpha_a_exp <- if (is.null(alpha_oracle$a_exp)) prep$omega * alpha_oracle$s else alpha_oracle$a_exp
  alpha_a_obs <- if (is.null(alpha_oracle$a_obs)) prep$omega * (1 - alpha_oracle$s) else alpha_oracle$a_obs
  beta_a_exp <- if (is.null(beta_oracle$a_exp)) prep$omega * beta_oracle$s else beta_oracle$a_exp
  beta_a_obs <- if (is.null(beta_oracle$a_obs)) prep$omega * (1 - beta_oracle$s) else beta_oracle$a_obs
  target_scale <- prep$target_scale
  risk_rescale <- prep$risk_rescale

  out <- list(
    selected = prep$arm_names[which(minimax$x == 1L)],
    x_opt = .make_named(minimax$x, names_out),
    gamma_opt = .make_named(minimax$gamma, names_out),
    s_opt = .make_named(minimax$s, names_out),
    a_exp_opt = .make_named(minimax_a_exp * target_scale, names_out),
    a_obs_opt = .make_named(minimax_a_obs * target_scale, names_out),
    n_opt = .make_named(n_opt, names_out),
    alpha_opt = minimax$alpha * risk_rescale,
    beta_opt = minimax$beta * risk_rescale,
    alpha_star = alpha_star * risk_rescale,
    beta_star = beta_star * risk_rescale,
    alpha_ratio = minimax$alpha_ratio,
    beta_ratio = minimax$beta_ratio,
    regret = minimax$regret,
    t_opt = minimax$t,
    solver_t_opt = if (is.null(minimax$solver_t)) minimax$t else minimax$solver_t,
    solver = solver,
    weight_mode = prep$weight_mode,
    beta_zero = beta_zero,
    psd_diagnostics = prep$psd_diagnostics,
    alpha_oracle = list(
      x = .make_named(alpha_oracle$x, names_out),
      gamma = .make_named(alpha_oracle$gamma, names_out),
      s = .make_named(alpha_oracle$s, names_out),
      a_exp = .make_named(alpha_a_exp * target_scale, names_out),
      a_obs = .make_named(alpha_a_obs * target_scale, names_out),
      n = .make_named(if (!is.null(alpha_oracle$n)) alpha_oracle$n else prep$compute_allocation(
        if (identical(prep$weight_mode, "a")) alpha_a_exp else alpha_oracle$s,
        x = alpha_oracle$x
      ), names_out),
      alpha = alpha_oracle$alpha * risk_rescale,
      beta = alpha_oracle$beta * risk_rescale
    ),
    beta_oracle = list(
      x = .make_named(beta_oracle$x, names_out),
      gamma = .make_named(beta_oracle$gamma, names_out),
      s = .make_named(beta_oracle$s, names_out),
      a_exp = .make_named(beta_a_exp * target_scale, names_out),
      a_obs = .make_named(beta_a_obs * target_scale, names_out),
      n = .make_named(if (!is.null(beta_oracle$n)) beta_oracle$n else prep$compute_allocation(
        if (identical(prep$weight_mode, "a")) beta_a_exp else beta_oracle$s,
        x = beta_oracle$x
      ), names_out),
      alpha = beta_oracle$alpha * risk_rescale,
      beta = beta_oracle$beta * risk_rescale
    ),
    feasible_sets_used = if (prep$use_sets) prep$sets else NULL,
    optimizer_result = minimax$result,
    gurobi_result = if (identical(solver, "gurobi")) minimax$result else NULL,
    inputs = list(
      Sigma_obs = prep$Sigma_obs,
      v2 = prep$v2,
      n_total = prep$n_total,
      omega = .make_named(prep$omega_original, names_out),
      costs = .make_named(prep$costs, names_out),
      bias_weights = .make_named(prep$bias_weights, names_out),
      h = prep$h,
      min_experiments = prep$min_experiments,
      n_min = prep$n_min,
      gamma_lower = .make_named(prep$gamma_lower, names_out),
      gamma_upper = .make_named(prep$gamma_upper, names_out),
      a_exp_lower = .make_named(prep$a_exp_lower * target_scale, names_out),
      a_exp_upper = .make_named(prep$a_exp_upper * target_scale, names_out),
      weight_mode = prep$weight_mode,
      force_gamma_unit_interval = prep$force_gamma_unit_interval,
      zero_tol = prep$zero_tol,
      psd_tol = prep$psd_tol,
      psd_action = prep$psd_action,
      psd_diagnostics = prep$psd_diagnostics,
      loading_structure = "identity",
      internal_target_scale = target_scale
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
