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
                                  n_min = 0,
                                  gamma_lower = NULL,
                                  gamma_upper = NULL,
                                  a_exp_lower = NULL,
                                  a_exp_upper = NULL,
                                  force_gamma_unit_interval = TRUE,
                                  solver = base::c("auto", "gurobi", "quadprog"),
                                  gurobi_params = list(),
                                  max_sets = 200000L,
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

  alpha_ratio <- .safe_ratio(alpha_oracle$alpha, alpha_star)
  beta_ratio <- .safe_ratio(alpha_oracle$beta, beta_star)
  regret <- max(alpha_ratio, beta_ratio)

  if (alpha_ratio < 1 - tol) {
    warning(sprintf("alpha_ratio = %.8f is below one by more than tol.", alpha_ratio),
            call. = FALSE)
  }

  names_out <- prep$arm_names
  alpha_a_exp <- if (is.null(alpha_oracle$a_exp)) prep$omega * alpha_oracle$s else alpha_oracle$a_exp
  alpha_a_obs <- if (is.null(alpha_oracle$a_obs)) prep$omega * (1 - alpha_oracle$s) else alpha_oracle$a_obs
  beta_a_exp <- if (is.null(beta_oracle$a_exp)) prep$omega * beta_oracle$s else beta_oracle$a_exp
  beta_a_obs <- if (is.null(beta_oracle$a_obs)) prep$omega * (1 - beta_oracle$s) else beta_oracle$a_obs
  target_scale <- prep$target_scale
  risk_rescale <- prep$risk_rescale
  out <- list(
    selected = prep$arm_names[which(alpha_oracle$x == 1L)],
    x_opt = .make_named(alpha_oracle$x, names_out),
    gamma_opt = .make_named(alpha_oracle$gamma, names_out),
    s_opt = .make_named(alpha_oracle$s, names_out),
    a_exp_opt = .make_named(alpha_a_exp * target_scale, names_out),
    a_obs_opt = .make_named(alpha_a_obs * target_scale, names_out),
    n_opt = .make_named(if (!is.null(alpha_oracle$n)) alpha_oracle$n else prep$compute_allocation(
      if (identical(prep$weight_mode, "a")) alpha_a_exp else alpha_oracle$s,
      x = alpha_oracle$x
    ), names_out),
    alpha_opt = alpha_oracle$alpha * risk_rescale,
    beta_opt = alpha_oracle$beta * risk_rescale,
    alpha_star = alpha_star * risk_rescale,
    beta_star = beta_star * risk_rescale,
    alpha_ratio = alpha_ratio,
    beta_ratio = beta_ratio,
    regret = regret,
    solver = solver,
    weight_mode = prep$weight_mode,
    beta_zero = beta_zero,
    psd_diagnostics = prep$psd_diagnostics,
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
    optimizer_result = alpha_oracle$result,
    gurobi_result = if (identical(solver, "gurobi")) alpha_oracle$result else NULL,
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
#'   experiments are multiplied by zero through `s = x * gamma`. In general-a
#'   mode, this can be supplied as a convenience and is translated into
#'   `a_exp = x * omega * gamma`.
#' @param a_exp Optional length-`p` vector of direct experimental linear weights
#'   for use when `force_gamma_unit_interval = FALSE`.
#' @param n_min Optional nonnegative lower bound on the sample size allocated to
#'   each selected experiment when `n_opt` is not supplied.
#' @param n_opt Optional explicit length-`p` sample allocation. If supplied,
#'   `evaluate_design()` uses this allocation directly instead of profiling over
#'   sample sizes.
#'
#' @return A list with `s`, `gamma`, `a_exp`, `a_obs`, `n_opt`, `alpha`, and
#'   `beta`.
#' @export
evaluate_design <- function(Sigma_obs,
                            v2,
                            n_total,
                            omega,
                            x,
                            gamma = NULL,
                            a_exp = NULL,
                            costs = NULL,
                            bias_weights = NULL,
                            gamma_lower = NULL,
                            gamma_upper = NULL,
                            a_exp_lower = NULL,
                            a_exp_upper = NULL,
                            force_gamma_unit_interval = TRUE,
                            n_min = 0,
                            n_opt = NULL,
                            zero_tol = 1e-9,
                            psd_tol = 1e-10,
                            psd_action = base::c("error", "project"),
                            tol = 1e-10,
                            c = NULL,
                            k = NULL) {
  prep <- .prepare_minimax_problem(
    Sigma_obs = Sigma_obs,
    v2 = v2,
    n_total = n_total,
    omega = omega,
    costs = costs,
    bias_weights = bias_weights,
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    a_exp_lower = a_exp_lower,
    a_exp_upper = a_exp_upper,
    force_gamma_unit_interval = force_gamma_unit_interval,
    n_min = n_min,
    h = length(omega),
    min_experiments = 0L,
    costs_alias = c,
    bias_weights_alias = k,
    zero_tol = zero_tol,
    psd_tol = psd_tol,
    psd_action = psd_action
  )
  p <- prep$p
  x <- as.numeric(x)
  if (length(x) != p) {
    stop("x must have length equal to length(omega).", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("x must contain finite values.", call. = FALSE)
  }
  if (any(!(x %in% base::c(0, 1)))) {
    stop("x must contain only 0 or 1.", call. = FALSE)
  }

  if (identical(prep$weight_mode, "a")) {
    if (is.null(a_exp)) {
      if (is.null(gamma)) {
        stop("Supply a_exp, or gamma to translate into a_exp, when force_gamma_unit_interval = FALSE.",
             call. = FALSE)
      }
      gamma <- as.numeric(gamma)
      if (length(gamma) != p || any(!is.finite(gamma))) {
        stop("gamma must be a finite numeric vector with length equal to length(omega).",
             call. = FALSE)
      }
      a_exp <- x * prep$omega * gamma
    } else {
      a_exp <- as.numeric(a_exp)
      if (length(a_exp) != p || any(!is.finite(a_exp))) {
        stop("a_exp must be a finite numeric vector with length equal to length(omega).",
             call. = FALSE)
      }
      # Direct weights are supplied and returned in the user's target units.
      a_exp <- a_exp / prep$target_scale
    }
    tol <- .validate_tolerance(tol)
    weight_scale <- max(base::c(abs(a_exp), abs(prep$omega),
                          abs(prep$a_exp_lower), abs(prep$a_exp_upper), 0))
    weight_tol <- tol * weight_scale
    if (any(abs(a_exp[x <= 0]) > weight_tol)) {
      stop("a_exp must be zero for unselected experiments.", call. = FALSE)
    }
    lower <- prep$a_exp_lower * x
    upper <- prep$a_exp_upper * x
    if (any(a_exp < lower - weight_tol | a_exp > upper + weight_tol)) {
      stop("a_exp violates the selection-adjusted a_exp_lower/a_exp_upper bounds.",
           call. = FALSE)
    }
    s <- .implied_gamma(a_exp, prep$omega)
    weight <- a_exp
  } else {
    if (is.null(gamma)) {
      stop("gamma must be supplied when force_gamma_unit_interval = TRUE.",
           call. = FALSE)
    }
    gamma <- as.numeric(gamma)
    if (length(gamma) != p || any(!is.finite(gamma))) {
      stop("gamma must be a finite numeric vector with length equal to length(omega).",
           call. = FALSE)
    }
    if (any(gamma < 0) || any(gamma > 1)) {
      stop("gamma must lie between zero and one when force_gamma_unit_interval = TRUE.",
           call. = FALSE)
    }
    selected <- x > 0.5
    bound_tol <- .validate_tolerance(tol)
    if (any(gamma[selected] < prep$gamma_lower[selected] - bound_tol |
            gamma[selected] > prep$gamma_upper[selected] + bound_tol)) {
      stop("gamma violates gamma_lower/gamma_upper for selected experiments.",
           call. = FALSE)
    }
    s <- x * gamma
    weight <- s
    a_exp <- prep$omega * s
  }
  a_obs <- prep$omega - a_exp
  names_out <- prep$arm_names
  allocation <- if (is.null(n_opt)) {
    .compute_allocation_with_floor(prep, weight, x = x, n_min = prep$n_min)
  } else {
    .validate_explicit_allocation(prep, n_opt, x = x, n_min = prep$n_min)
  }
  list(
    s = .make_named(s, names_out),
    gamma = .make_named(s, names_out),
    a_exp = .make_named(a_exp * prep$target_scale, names_out),
    a_obs = .make_named(a_obs * prep$target_scale, names_out),
    n_opt = .make_named(allocation, names_out),
    alpha = .compute_alpha_with_n(prep, weight, allocation) * prep$risk_rescale,
    beta = prep$compute_beta(weight) * prep$risk_rescale,
    psd_diagnostics = prep$psd_diagnostics
  )
}
