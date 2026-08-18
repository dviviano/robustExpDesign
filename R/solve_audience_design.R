.audience_cvxr_solve <- function(problem, cvxr_solver = NULL, cvxr_args = list()) {
  args <- base::c(list(object = problem), cvxr_args)
  if (!is.null(cvxr_solver)) {
    args$solver <- cvxr_solver
  }
  result <- do.call(CVXR::solve, args)
  ok <- result$status %in% c("optimal", "optimal_inaccurate", "OPTIMAL", "OPTIMAL_INACCURATE")
  if (!ok) {
    stop(sprintf("CVXR solver status: %s", result$status), call. = FALSE)
  }
  result
}

.audience_resolve_alias <- function(primary, aliases, primary_name, alias_names) {
  out <- primary
  for (i in seq_along(aliases)) {
    alias <- aliases[[i]]
    if (!is.null(alias)) {
      if (!is.null(out) && !isTRUE(all.equal(as.numeric(out), as.numeric(alias)))) {
        stop(sprintf("Use only one of %s, or supply identical values.",
                     paste(base::c(primary_name, alias_names), collapse = ", ")),
             call. = FALSE)
      }
      out <- alias
    }
  }
  out
}

.audience_psd_factor <- function(Sigma, name = "Sigma", tol = 1e-8) {
  Sigma <- as.matrix(Sigma)
  Sigma <- 0.5 * (Sigma + t(Sigma))
  eig <- eigen(Sigma, symmetric = TRUE)
  min_ev <- min(eig$values)
  if (min_ev < -sqrt(tol)) {
    stop(sprintf("%s must be positive semidefinite. Minimum eigenvalue is %.6g.",
                 name, min_ev), call. = FALSE)
  }
  vals <- pmax(eig$values, 0)
  diag(sqrt(vals), nrow = length(vals), ncol = length(vals)) %*% t(eig$vectors)
}

.audience_weighted_sum_expr <- function(x, weights) {
  weights <- as.numeric(weights)
  if (length(weights) == 0L) {
    return(0)
  }
  CVXR::sum_entries(CVXR::multiply(matrix(weights, ncol = 1L), x))
}

.audience_exp_conic_sum <- function(a_exp, n_var, v) {
  if (length(v) == 0L) {
    return(0)
  }
  terms <- vector("list", length(v))
  for (j in seq_along(v)) {
    terms[[j]] <- CVXR::quad_over_lin(v[j] * a_exp[j], n_var[j])
  }
  Reduce(`+`, terms)
}

.prepare_audience_problem <- function(Sigma_obs,
                                      v2,
                                      n_total,
                                      omega,
                                      B_grid,
                                      costs = NULL,
                                      bias_weights = NULL,
                                      Lambda_obs = NULL,
                                      Lambda_exp = NULL,
                                      h = NULL,
                                      x_max = NULL,
                                      feasible_sets = NULL,
                                      selection_constraints = NULL,
                                      min_experiments = 0L,
                                      costs_alias = NULL,
                                      bias_aliases = list(),
                                      tol = 1e-8) {
  .minimax_require("CVXR")

  costs <- .audience_resolve_alias(costs, list(costs_alias), "costs", "c")
  bias_weights <- .audience_resolve_alias(
    bias_weights, bias_aliases, "bias_weights", names(bias_aliases)
  )

  omega <- .numeric_vector(omega, "omega")
  d <- length(omega)
  if (d < 1L) {
    stop("omega must have positive length.", call. = FALSE)
  }

  Sigma_obs <- as.matrix(Sigma_obs)
  if (!is.numeric(Sigma_obs) || nrow(Sigma_obs) != ncol(Sigma_obs) || any(!is.finite(Sigma_obs))) {
    stop("Sigma_obs must be a finite numeric square matrix.", call. = FALSE)
  }
  Sigma_obs <- 0.5 * (Sigma_obs + t(Sigma_obs))
  po <- nrow(Sigma_obs)

  if (is.null(Lambda_obs)) {
    if (po != d) {
      stop("When Lambda_obs is omitted, Sigma_obs must be length(omega) by length(omega).",
           call. = FALSE)
    }
    Lambda_obs <- diag(d)
  } else {
    Lambda_obs <- as.matrix(Lambda_obs)
    if (!is.numeric(Lambda_obs) || any(!is.finite(Lambda_obs))) {
      stop("Lambda_obs must be a finite numeric matrix.", call. = FALSE)
    }
    if (nrow(Lambda_obs) != po || ncol(Lambda_obs) != d) {
      stop("Lambda_obs must have nrow(Sigma_obs) rows and length(omega) columns.",
           call. = FALSE)
    }
  }

  v2 <- .numeric_vector(v2, "v2")
  pe <- length(v2)
  if (pe < 1L) {
    stop("v2 must have positive length.", call. = FALSE)
  }
  if (any(v2 < 0)) {
    stop("v2 must be nonnegative.", call. = FALSE)
  }
  v <- sqrt(v2)

  if (is.null(Lambda_exp)) {
    if (pe != d) {
      stop("When Lambda_exp is omitted, v2 must have length equal to length(omega).",
           call. = FALSE)
    }
    Lambda_exp <- diag(d)
  } else {
    Lambda_exp <- as.matrix(Lambda_exp)
    if (!is.numeric(Lambda_exp) || any(!is.finite(Lambda_exp))) {
      stop("Lambda_exp must be a finite numeric matrix.", call. = FALSE)
    }
    if (nrow(Lambda_exp) != pe || ncol(Lambda_exp) != d) {
      stop("Lambda_exp must have length(v2) rows and length(omega) columns.",
           call. = FALSE)
    }
  }

  if (!is.numeric(n_total) || length(n_total) != 1L || !is.finite(n_total) || n_total <= 0) {
    stop("n_total must be a positive scalar budget.", call. = FALSE)
  }

  costs <- .numeric_vector(costs, "costs", pe, default = rep(1, pe))
  if (any(costs <= 0)) {
    stop("costs must be strictly positive.", call. = FALSE)
  }

  bias_weights <- .numeric_vector(bias_weights, "bias_weights", po, default = rep(1, po))
  if (any(bias_weights < 0)) {
    stop("bias_weights must be nonnegative.", call. = FALSE)
  }

  B_grid <- as.numeric(B_grid)
  if (length(B_grid) == 0L || any(!is.finite(B_grid)) || any(B_grid < 0)) {
    stop("B_grid must contain nonnegative finite values.", call. = FALSE)
  }
  B_grid <- sort(unique(B_grid))

  x_max <- .numeric_vector(x_max, "x_max", pe, default = rep(1, pe))
  if (any(!(x_max %in% c(0, 1)))) {
    stop("x_max must contain only 0 or 1.", call. = FALSE)
  }

  if (is.null(h)) {
    h <- pe
  }
  if (!.is_integerish_scalar(h) || h < 0 || h > pe) {
    stop("h must be an integer scalar between 0 and length(v2).", call. = FALSE)
  }
  h <- as.integer(round(h))

  if (!.is_integerish_scalar(min_experiments) || min_experiments < 0 || min_experiments > h) {
    stop("min_experiments must be an integer scalar between 0 and h.", call. = FALSE)
  }
  min_experiments <- as.integer(round(min_experiments))

  selection_constraints <- .validate_selection_constraints(selection_constraints, pe)
  set_info <- .process_feasible_sets(feasible_sets, pe, h, x_max, min_experiments,
                                     selection_constraints)

  experiment_names <- rownames(Lambda_exp)
  if (is.null(experiment_names) || any(experiment_names == "")) {
    experiment_names <- names(v2)
  }
  if (is.null(experiment_names) || any(experiment_names == "")) {
    experiment_names <- paste0("experiment_", seq_len(pe))
  }

  obs_names <- rownames(Lambda_obs)
  if (is.null(obs_names) || any(obs_names == "")) {
    obs_names <- rownames(Sigma_obs)
  }
  if (is.null(obs_names) || any(obs_names == "")) {
    obs_names <- paste0("obs_", seq_len(po))
  }

  param_names <- colnames(Lambda_obs)
  if (is.null(param_names) || any(param_names == "")) {
    param_names <- names(omega)
  }
  if (is.null(param_names) || any(param_names == "")) {
    param_names <- paste0("theta_", seq_len(d))
  }

  R_obs <- .audience_psd_factor(Sigma_obs, name = "Sigma_obs", tol = tol)

  c(
    list(
      p = pe,
      pe = pe,
      po = po,
      d = d,
      Sigma_obs = Sigma_obs,
      R_obs = R_obs,
      v2 = v2,
      v = v,
      n_total = n_total,
      omega = omega,
      B_grid = B_grid,
      costs = costs,
      bias_weights = bias_weights,
      Lambda_obs = Lambda_obs,
      Lambda_exp = Lambda_exp,
      h = h,
      x_max = x_max,
      min_experiments = min_experiments,
      selection_constraints = selection_constraints,
      arm_names = experiment_names,
      experiment_names = experiment_names,
      obs_names = obs_names,
      param_names = param_names
    ),
    set_info
  )
}

.audience_calibration_constraint <- function(prep, S, a_obs, a_exp = NULL) {
  lhs <- t(prep$Lambda_obs) %*% a_obs
  if (length(S) > 0L) {
    lhs <- lhs + t(prep$Lambda_exp[S, , drop = FALSE]) %*% a_exp
  }
  lhs == matrix(prep$omega, ncol = 1L)
}

.audience_full_x <- function(S, pe) {
  x <- rep(0, pe)
  if (length(S) > 0L) {
    x[S] <- 1
  }
  x
}

.audience_full_exp <- function(S, pe, a_exp) {
  out <- rep(0, pe)
  if (length(S) > 0L) {
    out[S] <- as.numeric(a_exp)
  }
  out
}

.audience_oracle_allocation <- function(prep, S, a_exp, tol = 1e-12) {
  out <- rep(0, prep$pe)
  if (length(S) == 0L) {
    return(out)
  }
  score <- abs(as.numeric(a_exp)) * prep$v[S]
  denom <- sum(score * sqrt(prep$costs[S]))
  if (denom <= tol) {
    out[S] <- prep$n_total / (length(S) * prep$costs[S])
  } else {
    out[S] <- prep$n_total * (score / sqrt(prep$costs[S])) / denom
  }
  out
}

.audience_eval <- function(prep, S, a_obs, a_exp = numeric(0), B = 0,
                           allocation = NULL, tol = 1e-12) {
  a_obs <- as.numeric(a_obs)
  a_exp <- as.numeric(a_exp)
  obs_var <- as.numeric(t(a_obs) %*% prep$Sigma_obs %*% a_obs)
  bias_radius <- sum(prep$bias_weights * abs(a_obs))
  beta <- bias_radius^2

  if (length(S) == 0L) {
    exp_var <- 0
    n_full <- rep(0, prep$pe)
  } else if (is.null(allocation)) {
    exp_radius <- sum(abs(a_exp) * prep$v[S] * sqrt(prep$costs[S]))
    exp_var <- exp_radius^2 / prep$n_total
    n_full <- .audience_oracle_allocation(prep, S, a_exp, tol = tol)
  } else {
    n_full <- as.numeric(allocation)
    exp_var <- 0
    for (j in seq_along(S)) {
      idx <- S[j]
      if (n_full[idx] <= tol) {
        if (abs(a_exp[j]) > sqrt(tol)) {
          exp_var <- Inf
        }
      } else {
        exp_var <- exp_var + (a_exp[j]^2) * prep$v2[idx] / n_full[idx]
      }
    }
  }

  risk <- obs_var + exp_var + B^2 * beta
  list(
    risk = as.numeric(risk),
    posterior_risk = as.numeric(risk),
    alpha = as.numeric(obs_var + exp_var),
    obs_variance = as.numeric(obs_var),
    exp_variance = as.numeric(exp_var),
    beta = as.numeric(beta),
    bias_radius = as.numeric(bias_radius),
    a_obs = .make_named(a_obs, prep$obs_names),
    a_exp = .make_named(.audience_full_exp(S, prep$pe, a_exp), prep$experiment_names),
    n = .make_named(n_full, prep$experiment_names),
    x = .make_named(.audience_full_x(S, prep$pe), prep$experiment_names)
  )
}

.audience_solve_pbar_fixed_set <- function(prep, S, B, cvxr_solver = NULL,
                                           cvxr_args = list(), tol = 1e-8) {
  po <- prep$po
  k <- length(S)
  a_obs <- CVXR::Variable(po)
  z_obs <- CVXR::Variable(po)

  constraints <- list(
    z_obs >= a_obs,
    z_obs >= -a_obs
  )

  exp_expr <- 0
  if (k > 0L) {
    a_exp <- CVXR::Variable(k)
    r_exp <- CVXR::Variable(k)
    constraints <- base::c(
      constraints,
      list(
        .audience_calibration_constraint(prep, S, a_obs, a_exp),
        r_exp >= a_exp,
        r_exp >= -a_exp
      )
    )
    exp_radius <- .audience_weighted_sum_expr(
      r_exp, prep$v[S] * sqrt(prep$costs[S])
    )
    exp_expr <- CVXR::square(exp_radius) / prep$n_total
  } else {
    a_exp <- NULL
    constraints <- base::c(
      constraints,
      list(.audience_calibration_constraint(prep, S, a_obs, NULL))
    )
  }

  obs_expr <- CVXR::sum_squares(prep$R_obs %*% a_obs)
  bias_expr <- CVXR::square(.audience_weighted_sum_expr(z_obs, prep$bias_weights))
  objective <- CVXR::Minimize(obs_expr + exp_expr + B^2 * bias_expr)
  result <- .audience_cvxr_solve(
    CVXR::Problem(objective, constraints), cvxr_solver = cvxr_solver, cvxr_args = cvxr_args
  )

  a_obs_val <- as.numeric(result$getValue(a_obs))
  a_exp_val <- if (k > 0L) as.numeric(result$getValue(a_exp)) else numeric(0)
  eval <- .audience_eval(prep, S, a_obs_val, a_exp_val, B = B, allocation = NULL, tol = tol)
  eval$status <- result$status
  eval$objective_value <- result$value
  eval$set <- S
  eval
}

.audience_solve_beta_fixed_set <- function(prep, S, cvxr_solver = NULL,
                                           cvxr_args = list(), tol = 1e-8) {
  po <- prep$po
  k <- length(S)
  a_obs <- CVXR::Variable(po)
  z_obs <- CVXR::Variable(po)

  constraints <- list(
    z_obs >= a_obs,
    z_obs >= -a_obs
  )

  if (k > 0L) {
    a_exp <- CVXR::Variable(k)
    constraints <- base::c(
      constraints,
      list(.audience_calibration_constraint(prep, S, a_obs, a_exp))
    )
  } else {
    a_exp <- NULL
    constraints <- base::c(
      constraints,
      list(.audience_calibration_constraint(prep, S, a_obs, NULL))
    )
  }

  radius <- .audience_weighted_sum_expr(z_obs, prep$bias_weights)
  objective <- CVXR::Minimize(radius)
  result <- .audience_cvxr_solve(
    CVXR::Problem(objective, constraints), cvxr_solver = cvxr_solver, cvxr_args = cvxr_args
  )

  a_obs_val <- as.numeric(result$getValue(a_obs))
  a_exp_val <- if (k > 0L) as.numeric(result$getValue(a_exp)) else numeric(0)
  eval <- .audience_eval(prep, S, a_obs_val, a_exp_val, B = 0, allocation = NULL, tol = tol)
  eval$status <- result$status
  eval$objective_value <- result$value
  eval$set <- S
  eval
}

.audience_solve_grid_fixed_set <- function(prep, S, P_star, beta_star,
                                           cvxr_solver = NULL,
                                           cvxr_args = list(),
                                           sample_floor = 1e-8,
                                           tol = 1e-8) {
  po <- prep$po
  G <- length(prep$B_grid)
  k <- length(S)

  r_var <- CVXR::Variable(1)
  a_obs <- CVXR::Variable(po, G)
  z_obs <- CVXR::Variable(po, G)
  a_inf_obs <- CVXR::Variable(po)
  z_inf <- CVXR::Variable(po)

  constraints <- list(r_var >= 0)

  if (k > 0L) {
    n_var <- CVXR::Variable(k)
    a_exp <- CVXR::Variable(k, G)
    a_inf_exp <- CVXR::Variable(k)
    budget_expr <- .audience_weighted_sum_expr(n_var, prep$costs[S])
    constraints <- base::c(
      constraints,
      list(n_var >= sample_floor, budget_expr == prep$n_total)
    )
  } else {
    n_var <- NULL
    a_exp <- NULL
    a_inf_exp <- NULL
  }

  for (g in seq_len(G)) {
    a_obs_g <- a_obs[, g]
    z_obs_g <- z_obs[, g]
    constraints <- base::c(
      constraints,
      list(
        z_obs_g >= a_obs_g,
        z_obs_g >= -a_obs_g
      )
    )

    if (k > 0L) {
      a_exp_g <- a_exp[, g]
      constraints <- base::c(
        constraints,
        list(.audience_calibration_constraint(prep, S, a_obs_g, a_exp_g))
      )
      exp_expr <- .audience_exp_conic_sum(a_exp_g, n_var, prep$v[S])
    } else {
      a_exp_g <- NULL
      constraints <- base::c(
        constraints,
        list(.audience_calibration_constraint(prep, S, a_obs_g, NULL))
      )
      exp_expr <- 0
    }

    obs_expr <- CVXR::sum_squares(prep$R_obs %*% a_obs_g)
    bias_expr <- CVXR::square(.audience_weighted_sum_expr(z_obs_g, prep$bias_weights))
    risk_expr <- obs_expr + exp_expr + prep$B_grid[g]^2 * bias_expr
    constraints <- base::c(constraints, list(risk_expr <= P_star[g] * r_var))
  }

  constraints <- base::c(
    constraints,
    list(
      z_inf >= a_inf_obs,
      z_inf >= -a_inf_obs
    )
  )
  if (k > 0L) {
    constraints <- base::c(
      constraints,
      list(.audience_calibration_constraint(prep, S, a_inf_obs, a_inf_exp))
    )
  } else {
    constraints <- base::c(
      constraints,
      list(.audience_calibration_constraint(prep, S, a_inf_obs, NULL))
    )
  }
  inf_radius <- .audience_weighted_sum_expr(z_inf, prep$bias_weights)
  if (beta_star <= tol) {
    constraints <- base::c(constraints, list(inf_radius <= 0))
  } else {
    constraints <- base::c(constraints, list(CVXR::square(inf_radius) <= beta_star * r_var))
  }

  result <- .audience_cvxr_solve(
    CVXR::Problem(CVXR::Minimize(r_var), constraints),
    cvxr_solver = cvxr_solver,
    cvxr_args = cvxr_args
  )

  r_val <- as.numeric(result$getValue(r_var))
  a_obs_val <- as.matrix(result$getValue(a_obs))
  if (G == 1L) {
    a_obs_val <- matrix(a_obs_val, nrow = po, ncol = 1L)
  }
  a_exp_val <- if (k > 0L) as.matrix(result$getValue(a_exp)) else matrix(0, nrow = 0L, ncol = G)
  if (k > 0L && G == 1L) {
    a_exp_val <- matrix(a_exp_val, nrow = k, ncol = 1L)
  }

  n_full <- rep(0, prep$pe)
  if (k > 0L) {
    n_full[S] <- as.numeric(result$getValue(n_var))
  }

  by_B <- vector("list", G)
  risk <- numeric(G)
  posterior_ratio <- numeric(G)
  for (g in seq_len(G)) {
    eval_g <- .audience_eval(
      prep, S,
      a_obs = a_obs_val[, g],
      a_exp = if (k > 0L) a_exp_val[, g] else numeric(0),
      B = prep$B_grid[g],
      allocation = n_full,
      tol = tol
    )
    eval_g$B <- prep$B_grid[g]
    eval_g$P_star <- P_star[g]
    eval_g$posterior_ratio <- .safe_ratio(eval_g$risk, P_star[g], tol = tol)
    by_B[[g]] <- eval_g
    risk[g] <- eval_g$risk
    posterior_ratio[g] <- eval_g$posterior_ratio
  }

  a_inf_obs_val <- as.numeric(result$getValue(a_inf_obs))
  a_inf_exp_val <- if (k > 0L) as.numeric(result$getValue(a_inf_exp)) else numeric(0)
  inf_eval <- .audience_eval(
    prep, S, a_obs = a_inf_obs_val, a_exp = a_inf_exp_val,
    B = 0, allocation = NULL, tol = tol
  )
  inf_eval$beta_ratio <- .safe_ratio(inf_eval$beta, beta_star, tol = tol)

  x <- .audience_full_x(S, prep$pe)
  list(
    set = S,
    selected = prep$experiment_names[which(x == 1L)],
    x = .make_named(x, prep$experiment_names),
    n = .make_named(n_full, prep$experiment_names),
    regret = max(base::c(posterior_ratio, inf_eval$beta_ratio), na.rm = TRUE),
    t = r_val,
    posterior_ratio = posterior_ratio,
    risk = risk,
    by_B = by_B,
    large_B = inf_eval,
    status = result$status,
    objective_value = result$value,
    result = result
  )
}

#' Solve the Section 4.1.1 audience-regret design problem on a grid of B values
#'
#' @description
#' `solve_audience_design()` implements the finite-grid profiled posterior-risk
#' program in Appendix D.1. It first computes the oracle posterior risk
#' `P_star(B_g)` at each supplied grid value and the large-B bias benchmark
#' `beta_star`. It then chooses one experiment set and one sample allocation that
#' minimize the maximum posterior-risk ratio over the grid, together with the
#' large-B bias-ratio constraint.
#'
#' By default the function uses the same direct-parameter design as
#' [solve_minimax_design()], with `Lambda_obs = I` and `Lambda_exp = I`. For a
#' more general linear reported-statistic problem, pass `Lambda_obs` and
#' `Lambda_exp`; the calibration constraint is
#' `t(Lambda_obs) %*% a_obs + t(Lambda_exp[E, ]) %*% a_exp = omega`.
#'
#' Weighted bias radii are supported through `bias_weights`, `kappa`, or `k`.
#' For the weighted l-infinity ambiguity set `|b_l| <= bias_weights[l] * B`, the
#' dual radius used in Appendix D.1 is `sum_l bias_weights[l] * |a_obs_l|`.
#'
#' @inheritParams solve_minimax_design
#' @param B_grid Numeric vector of nonnegative prior-scale values `B_g`.
#' @param Lambda_obs Optional `p_o` by `d` loading matrix for observational
#'   reported statistics. Defaults to the identity matrix.
#' @param Lambda_exp Optional `p_e` by `d` loading matrix for experimental
#'   reported statistics. Defaults to the identity matrix.
#' @param k,kappa Aliases for `bias_weights`, included for script compatibility.
#'   `kappa` matches the notation in the paper.
#' @param cvxr_solver Optional CVXR solver name, e.g. `"CLARABEL"`, `"ECOS"`,
#'   or `"SCS"` depending on what is installed.
#' @param cvxr_args Optional named list passed to `CVXR::solve()`.
#' @param sample_floor Small positive lower bound for selected-arm sample-size
#'   variables in the conic solve. This avoids numerical division by zero in
#'   quadratic-over-linear constraints.
#'
#' @return An object of class `audience_design`.
#' @export
solve_audience_design <- function(Sigma_obs,
                                  v2,
                                  n_total,
                                  omega,
                                  B_grid,
                                  costs = NULL,
                                  bias_weights = NULL,
                                  Lambda_obs = NULL,
                                  Lambda_exp = NULL,
                                  h = length(v2),
                                  x_max = NULL,
                                  feasible_sets = NULL,
                                  selection_constraints = NULL,
                                  min_experiments = 0L,
                                  cvxr_solver = NULL,
                                  cvxr_args = list(),
                                  max_sets = 200000L,
                                  sample_floor = 1e-8,
                                  tol = 1e-7,
                                  c = NULL,
                                  k = NULL,
                                  kappa = NULL) {
  prep <- .prepare_audience_problem(
    Sigma_obs = Sigma_obs,
    v2 = v2,
    n_total = n_total,
    omega = omega,
    B_grid = B_grid,
    costs = costs,
    bias_weights = bias_weights,
    Lambda_obs = Lambda_obs,
    Lambda_exp = Lambda_exp,
    h = h,
    x_max = x_max,
    feasible_sets = feasible_sets,
    selection_constraints = selection_constraints,
    min_experiments = min_experiments,
    costs_alias = c,
    bias_aliases = list(k = k, kappa = kappa),
    tol = tol
  )

  sets <- .enumerate_experiment_sets(prep, max_sets = max_sets)
  G <- length(prep$B_grid)

  oracle_by_B <- vector("list", G)
  P_star <- rep(Inf, G)
  for (g in seq_len(G)) {
    B <- prep$B_grid[g]
    sols <- vector("list", length(sets))
    vals <- rep(Inf, length(sets))
    for (i in seq_along(sets)) {
      sols[[i]] <- tryCatch(
        .audience_solve_pbar_fixed_set(
          prep, sets[[i]], B = B, cvxr_solver = cvxr_solver,
          cvxr_args = cvxr_args, tol = tol
        ),
        error = function(e) structure(list(error = conditionMessage(e)), class = "audience_error")
      )
      if (!inherits(sols[[i]], "audience_error")) {
        vals[i] <- sols[[i]]$risk
      }
    }
    best_i <- which.min(vals)
    if (!is.finite(vals[best_i])) {
      stop(sprintf("No feasible audience oracle was found for B_grid[%d] = %s.",
                   g, format(B)), call. = FALSE)
    }
    oracle_by_B[[g]] <- sols[[best_i]]
    P_star[g] <- vals[best_i]
  }

  beta_solutions <- vector("list", length(sets))
  beta_vals <- rep(Inf, length(sets))
  for (i in seq_along(sets)) {
    beta_solutions[[i]] <- tryCatch(
      .audience_solve_beta_fixed_set(
        prep, sets[[i]], cvxr_solver = cvxr_solver,
        cvxr_args = cvxr_args, tol = tol
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "audience_error")
    )
    if (!inherits(beta_solutions[[i]], "audience_error")) {
      beta_vals[i] <- beta_solutions[[i]]$beta
    }
  }
  beta_best_i <- which.min(beta_vals)
  if (!is.finite(beta_vals[beta_best_i])) {
    stop("No feasible large-B bias oracle was found.", call. = FALSE)
  }
  beta_star <- beta_vals[beta_best_i]
  beta_oracle <- beta_solutions[[beta_best_i]]

  candidates <- vector("list", length(sets))
  regret_vals <- rep(Inf, length(sets))
  errors <- rep(NA_character_, length(sets))
  for (i in seq_along(sets)) {
    candidates[[i]] <- tryCatch(
      .audience_solve_grid_fixed_set(
        prep, sets[[i]], P_star = P_star, beta_star = beta_star,
        cvxr_solver = cvxr_solver, cvxr_args = cvxr_args,
        sample_floor = sample_floor, tol = tol
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "audience_error")
    )
    if (inherits(candidates[[i]], "audience_error")) {
      errors[i] <- candidates[[i]]$error
    } else {
      regret_vals[i] <- candidates[[i]]$regret
    }
  }

  best_i <- which.min(regret_vals)
  if (!is.finite(regret_vals[best_i])) {
    stop("No feasible audience-regret solution was found.", call. = FALSE)
  }
  best <- candidates[[best_i]]

  set_label <- vapply(sets, function(S) {
    if (length(S) == 0L) "<empty>" else paste(prep$experiment_names[S], collapse = ", ")
  }, character(1))

  out <- list(
    selected = best$selected,
    x_opt = best$x,
    n_opt = best$n,
    regret = best$regret,
    t_opt = best$t,
    B_grid = prep$B_grid,
    posterior_ratio = best$posterior_ratio,
    posterior_risk = best$risk,
    P_star = P_star,
    beta_star = beta_star,
    beta_ratio_large_B = best$large_B$beta_ratio,
    a_by_B = best$by_B,
    a_large_B = best$large_B,
    oracle_by_B = oracle_by_B,
    beta_oracle = beta_oracle,
    feasible_sets_used = if (prep$use_sets) prep$sets else NULL,
    candidate_diagnostics = data.frame(
      candidate_index = seq_along(sets),
      set = set_label,
      regret = regret_vals,
      error = errors,
      stringsAsFactors = FALSE
    ),
    optimizer_result = best$result,
    inputs = list(
      Sigma_obs = prep$Sigma_obs,
      v2 = .make_named(prep$v2, prep$experiment_names),
      n_total = prep$n_total,
      omega = .make_named(prep$omega, prep$param_names),
      costs = .make_named(prep$costs, prep$experiment_names),
      bias_weights = .make_named(prep$bias_weights, prep$obs_names),
      Lambda_obs = prep$Lambda_obs,
      Lambda_exp = prep$Lambda_exp,
      h = prep$h,
      min_experiments = prep$min_experiments,
      cvxr_solver = cvxr_solver,
      sample_floor = sample_floor
    )
  )
  class(out) <- c("audience_design", "list")
  out
}

#' @rdname solve_audience_design
#' @export
solve_audience_regret_design <- solve_audience_design

#' @export
print.audience_design <- function(x, ...) {
  cat("Audience-regret design\n")
  cat("Selected experiments: ", .format_selected(x$x_opt, names(x$x_opt)), "\n", sep = "")
  cat("Regret: ", format(x$regret, digits = 6), "\n", sep = "")
  cat("Large-B bias ratio: ", format(x$beta_ratio_large_B, digits = 6), "\n", sep = "")
  invisible(x)
}
