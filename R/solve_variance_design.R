#' Solve the variance-optimal design problem
#'
#' @description
#' `solve_variance_design()` computes the design that minimizes the variance
#' index `alpha(s)` in equation (9), subject to the same feasibility restrictions
#' used by [solve_minimax_design()]. This is the variance-optimal benchmark used
#' for regret diagnostics.
#'
#' @inheritParams solve_minimax_design
#'
#' @return An object of class `variance_design`.
#' @export
solve_variance_design <- function(Sigma_obs,
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

  alpha_ratio <- .safe_ratio(alpha_oracle$alpha, alpha_star)
  beta_ratio <- .safe_ratio(alpha_oracle$beta, beta_star)
  regret <- max(alpha_ratio, beta_ratio)

  if (alpha_ratio < 1 - tol) {
    warning(sprintf("alpha_ratio = %.8f is below one by more than tol.", alpha_ratio),
            call. = FALSE)
  }

  names_out <- prep$arm_names
  out <- list(
    selected = prep$arm_names[which(alpha_oracle$x == 1L)],
    x_opt = .make_named(alpha_oracle$x, names_out),
    gamma_opt = .make_named(alpha_oracle$gamma, names_out),
    s_opt = .make_named(alpha_oracle$s, names_out),
    n_opt = .make_named(prep$compute_allocation(alpha_oracle$s), names_out),
    alpha_opt = alpha_oracle$alpha,
    beta_opt = alpha_oracle$beta,
    alpha_star = alpha_star,
    beta_star = beta_star,
    alpha_ratio = alpha_ratio,
    beta_ratio = beta_ratio,
    regret = regret,
    solver = solver,
    beta_oracle = list(
      x = .make_named(beta_oracle$x, names_out),
      gamma = .make_named(beta_oracle$gamma, names_out),
      s = .make_named(beta_oracle$s, names_out),
      n = .make_named(prep$compute_allocation(beta_oracle$s), names_out),
      alpha = beta_oracle$alpha,
      beta = beta_oracle$beta
    ),
    feasible_sets_used = if (prep$use_sets) prep$sets else NULL,
    optimizer_result = alpha_oracle$result,
    gurobi_result = if (identical(solver, "gurobi")) alpha_oracle$result else NULL,
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
  class(out) <- base::c("variance_design", "list")
  out
}

#' @export
print.variance_design <- function(x, ...) {
  cat("Variance-optimal design\n")
  cat("Selected experiments: ", .format_selected(x$x_opt, names(x$x_opt)), "\n", sep = "")
  cat("Solver: ", if (is.null(x$solver)) "unknown" else x$solver, "\n", sep = "")
  cat("Variance ratio: ", format(x$alpha_ratio, digits = 6), "\n", sep = "")
  cat("Bias ratio: ", format(x$beta_ratio, digits = 6), "\n", sep = "")
  invisible(x)
}

#' Evaluate variance, bias, and allocation for a supplied design
#'
#' @description
#' Computes the variance index `alpha(s)`, the bias index `beta(s)`, and the
#' cost-aware allocation for user-supplied selection and shrinkage vectors.
#'
#' @inheritParams solve_minimax_design
#' @param x Numeric or logical length-`p` selection vector.
#' @param gamma Numeric length-`p` shrinkage vector. Entries for unselected
#'   experiments are multiplied by zero through `s = x * gamma`.
#'
#' @return A list with `s`, `n_opt`, `alpha`, and `beta`.
#' @export
evaluate_design <- function(Sigma_obs,
                            v2,
                            n_total,
                            omega,
                            x,
                            gamma,
                            costs = NULL,
                            bias_weights = NULL,
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
    h = length(omega),
    min_experiments = 0L,
    costs_alias = c,
    bias_weights_alias = k,
    kappa_alias = kappa
  )
  p <- prep$p
  x <- as.numeric(x)
  gamma <- as.numeric(gamma)
  if (length(x) != p || length(gamma) != p) {
    stop("x and gamma must have length equal to length(omega).", call. = FALSE)
  }
  if (any(!is.finite(x)) || any(!is.finite(gamma))) {
    stop("x and gamma must contain finite values.", call. = FALSE)
  }
  if (any(x < 0 | x > 1) || any(gamma < 0 | gamma > 1)) {
    stop("x and gamma entries must lie between zero and one.", call. = FALSE)
  }
  s <- x * gamma
  names_out <- prep$arm_names
  list(
    s = .make_named(s, names_out),
    n_opt = .make_named(prep$compute_allocation(s), names_out),
    alpha = prep$compute_alpha(s),
    beta = prep$compute_beta(s)
  )
}
