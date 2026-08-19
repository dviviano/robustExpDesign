.audience_check_psd <- function(Sigma, tol = 1e-10,
                                action = base::c("error", "project")) {
  checked <- .validate_psd_matrix(
    Sigma, name = "Sigma_obs", tol = tol, action = match.arg(action)
  )
  out <- checked$matrix
  attr(out, "psd_diagnostics") <- checked$diagnostics
  out
}

.audience_default_lambda_grid <- function() {
  seq(0, 1, length.out = 50L)
}

.audience_validate_lambda_grid <- function(lambda_grid, tol = 1e-8,
                                           require_endpoints = TRUE,
                                           supplied = TRUE) {
  if (is.null(lambda_grid) || is.logical(lambda_grid) ||
      is.factor(lambda_grid) || !is.atomic(lambda_grid) ||
      !is.null(dim(lambda_grid)) || !is.numeric(lambda_grid)) {
    stop("lambda_grid must be a finite numeric vector.", call. = FALSE)
  }
  if (length(lambda_grid) == 0L) {
    stop("lambda_grid must contain at least two values.", call. = FALSE)
  }
  original <- as.numeric(lambda_grid)
  if (any(!is.finite(original))) {
    stop("lambda_grid must contain only finite values.", call. = FALSE)
  }
  if (any(original < -tol | original > 1 + tol)) {
    stop("lambda_grid values must lie in [0, 1].", call. = FALSE)
  }

  canonical <- original
  snap_zero <- abs(canonical) <= tol
  snap_one <- abs(canonical - 1) <= tol
  canonical[snap_zero] <- 0
  canonical[snap_one] <- 1
  if (any(canonical < 0 | canonical > 1)) {
    stop("lambda_grid values must lie in [0, 1].", call. = FALSE)
  }
  canonical <- sort(canonical)
  canonical_unique <- unique(canonical)
  duplicates_removed <- length(canonical_unique) < length(canonical)
  endpoints_snapped <- any(canonical != sort(original))
  if (duplicates_removed) {
    warning("Duplicate lambda_grid values were removed after sorting and endpoint snapping.",
            call. = FALSE)
  }
  if (length(canonical_unique) < 2L) {
    stop("lambda_grid must contain at least two distinct values.", call. = FALSE)
  }
  if (isTRUE(require_endpoints) &&
      !(identical(canonical_unique[1L], 0) &&
        identical(canonical_unique[length(canonical_unique)], 1))) {
    stop("lambda_grid must include both endpoints 0 and 1.", call. = FALSE)
  }

  gaps <- diff(canonical_unique)
  list(
    values = canonical_unique,
    diagnostics = list(
      source = if (isTRUE(supplied)) "user" else "default",
      input_points = length(original),
      canonical_points = length(canonical_unique),
      min = min(canonical_unique),
      max = max(canonical_unique),
      max_gap = if (length(gaps) > 0L) max(gaps) else NA_real_,
      duplicates_removed = duplicates_removed,
      endpoints_snapped = endpoints_snapped
    )
  )
}

.audience_scalarize <- function(alpha, beta, lambda) {
  if (lambda == 0) return(as.numeric(alpha))
  if (lambda == 1) return(as.numeric(beta))
  as.numeric((1 - lambda) * alpha + lambda * beta)
}

.audience_prepare <- function(Sigma_obs, v2, n_total, omega, lambda_grid,
                              costs = NULL, bias_weights = NULL,
                              h = NULL, x_max = NULL, feasible_sets = NULL,
                              selection_constraints = NULL,
                              min_experiments = 1L,
                              n_min = 0,
                              costs_alias = NULL,
                              bias_weights_alias = NULL,
                              lambda_grid_supplied = TRUE,
                              require_lambda_endpoints = TRUE,
                              max_sets = 200000L,
                              zero_tol = 1e-9,
                              psd_tol = 1e-10,
                              psd_action = base::c("error", "project"),
                              tol = 1e-8) {
  .minimax_require("Matrix")

  if (!is.null(costs_alias)) {
    if (!is.null(costs) &&
        !isTRUE(all.equal(as.numeric(costs), as.numeric(costs_alias)))) {
      stop("Use either costs or c, or supply identical values for both.",
           call. = FALSE)
    }
    costs <- costs_alias
  }
  if (!is.null(bias_weights_alias)) {
    if (!is.null(bias_weights) &&
        !isTRUE(all.equal(as.numeric(bias_weights),
                          as.numeric(bias_weights_alias)))) {
      stop("Use either bias_weights or k, or supply identical values for both.",
           call. = FALSE)
    }
    bias_weights <- bias_weights_alias
  }

  zero_tol <- .validate_tolerance(zero_tol, "zero_tol")
  psd_action <- match.arg(psd_action)
  n_min <- as.numeric(n_min)
  if (length(n_min) != 1L || !is.finite(n_min) || n_min < 0) {
    stop("n_min must be a nonnegative finite scalar.", call. = FALSE)
  }

  omega_names <- names(omega)
  omega_original <- .numeric_vector(omega, "omega")
  p <- length(omega_original)
  target_scale <- max(abs(omega_original))
  if (target_scale == 0) target_scale <- 1
  # Optimize in normalized target units. All risks are rescaled on output.
  omega <- omega_original / target_scale

  Sigma_obs <- as.matrix(Sigma_obs)
  if (!is.numeric(Sigma_obs) || any(dim(Sigma_obs) != base::c(p, p))) {
    stop("Sigma_obs must be a p by p numeric matrix with p = length(omega).",
         call. = FALSE)
  }
  if (any(!is.finite(Sigma_obs))) {
    stop("Sigma_obs must contain only finite values.", call. = FALSE)
  }
  checked_sigma <- .validate_psd_matrix(
    Sigma_obs, name = "Sigma_obs", tol = psd_tol, action = psd_action
  )
  Sigma_obs <- checked_sigma$matrix

  v2 <- .numeric_vector(v2, "v2", p)
  if (any(v2 < 0)) stop("v2 must be nonnegative.", call. = FALSE)
  if (!is.numeric(n_total) || length(n_total) != 1L ||
      !is.finite(n_total) || n_total <= 0) {
    stop("n_total must be a positive scalar budget.", call. = FALSE)
  }
  lambda_info <- .audience_validate_lambda_grid(
    lambda_grid, tol = tol,
    require_endpoints = require_lambda_endpoints,
    supplied = lambda_grid_supplied
  )
  lambda_grid <- lambda_info$values

  costs <- .numeric_vector(costs, "costs", p, default = rep(1, p))
  if (any(costs <= 0)) stop("costs must be strictly positive.", call. = FALSE)
  bias_weights <- .numeric_vector(
    bias_weights, "bias_weights", p, default = rep(1, p)
  )
  if (any(bias_weights < 0)) {
    stop("bias_weights must be nonnegative.", call. = FALSE)
  }
  x_max <- .numeric_vector(x_max, "x_max", p, default = rep(1, p))
  if (any(!(x_max %in% base::c(0, 1)))) {
    stop("x_max must contain only 0 or 1.", call. = FALSE)
  }
  if (is.null(h)) h <- p
  if (!.is_integerish_scalar(h) || h < 0 || h > p) {
    stop("h must be an integer scalar between 0 and length(omega).",
         call. = FALSE)
  }
  h <- as.integer(round(h))
  if (!.is_integerish_scalar(min_experiments) ||
      min_experiments < 0 || min_experiments > h) {
    stop("min_experiments must be an integer scalar between 0 and h.",
         call. = FALSE)
  }
  min_experiments <- as.integer(round(min_experiments))

  selection_constraints <- .validate_selection_constraints(selection_constraints, p)
  set_info <- .process_feasible_sets(
    feasible_sets, p, h, x_max, min_experiments, selection_constraints
  )
  prep_sets <- base::c(
    list(p = p, h = h, x_max = x_max,
         min_experiments = min_experiments,
         selection_constraints = selection_constraints),
    set_info
  )
  sets <- .enumerate_experiment_sets(prep_sets, max_sets = max_sets)
  sets <- Filter(function(S) length(S) > 0L, sets)
  budget_tol <- 1e-10 * n_total
  sets <- Filter(function(S) {
    sum(costs[S] * n_min) <= n_total + budget_tol
  }, sets)
  if (length(sets) == 0L) {
    stop("No non-empty feasible experiment set satisfies the selection restrictions and n_min budget floor.",
         call. = FALSE)
  }

  arm_names <- omega_names
  if (is.null(arm_names) || any(arm_names == "")) {
    arm_names <- paste0("arm_", seq_len(p))
  }

  alpha_scale <- as.numeric(t(omega) %*% Sigma_obs %*% omega) +
    (sum(abs(omega) * sqrt(v2) * sqrt(costs))^2) / n_total
  beta_scale <- (sum(bias_weights * abs(omega)))^2
  risk_scale <- (1 - lambda_grid) * alpha_scale + lambda_grid * beta_scale

  list(
    p = p,
    Sigma_obs = Sigma_obs,
    psd_diagnostics = checked_sigma$diagnostics,
    psd_tol = psd_tol,
    psd_action = psd_action,
    v2 = v2,
    v = sqrt(v2),
    n_total = n_total,
    n_min = n_min,
    omega = omega,
    omega_original = omega_original,
    target_scale = target_scale,
    risk_rescale = target_scale^2,
    lambda_grid = lambda_grid,
    lambda_grid_diagnostics = lambda_info$diagnostics,
    costs = costs,
    bias_weights = bias_weights,
    h = h,
    x_max = x_max,
    min_experiments = min_experiments,
    selection_constraints = selection_constraints,
    sets = sets,
    arm_names = arm_names,
    q_obs = -2 * as.numeric(Sigma_obs %*% omega),
    const_obs = as.numeric(t(omega) %*% Sigma_obs %*% omega),
    M_bias = bias_weights %o% bias_weights,
    alpha_scale = alpha_scale,
    beta_scale = beta_scale,
    risk_scale = risk_scale,
    zero_tol = zero_tol,
    loading_structure = "identity"
  )
}

.audience_abs_constraints <- function(A, rhs, sense, nvar, idx_a, idx_z,
                                      omega) {
  p <- length(omega)
  for (j in seq_len(p)) {
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- 1
    row[1L, idx_z[j]] <- -1
    A <- rbind(A, row)
    rhs <- base::c(rhs, omega[j])
    sense <- base::c(sense, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- -1
    row[1L, idx_z[j]] <- -1
    A <- rbind(A, row)
    rhs <- base::c(rhs, -omega[j])
    sense <- base::c(sense, "<")
  }
  list(A = A, rhs = rhs, sense = sense)
}

.audience_r_abs_constraints <- function(A, rhs, sense, nvar, idx_a, idx_r) {
  p <- length(idx_a)
  for (j in seq_len(p)) {
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- 1
    row[1L, idx_r[j]] <- -1
    A <- rbind(A, row)
    rhs <- base::c(rhs, 0)
    sense <- base::c(sense, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- -1
    row[1L, idx_r[j]] <- -1
    A <- rbind(A, row)
    rhs <- base::c(rhs, 0)
    sense <- base::c(sense, "<")
  }
  list(A = A, rhs = rhs, sense = sense)
}

.audience_set_bounds <- function(prep, S, nvar, idx_a = NULL, idx_n = NULL,
                                 idx_r = NULL, idx_u = NULL,
                                 idx_z = NULL, n_min = 0) {
  p <- prep$p
  selected <- rep(FALSE, p)
  selected[S] <- TRUE
  lb <- rep(-Inf, nvar)
  ub <- rep(Inf, nvar)
  if (!is.null(idx_a)) {
    lb[idx_a[!selected]] <- 0
    ub[idx_a[!selected]] <- 0
  }
  if (!is.null(idx_n)) {
    lb[idx_n] <- 0
    ub[idx_n] <- 0
    lb[idx_n[selected]] <- n_min
    ub[idx_n[selected]] <- prep$n_total / prep$costs[selected]
  }
  if (!is.null(idx_r)) {
    lb[idx_r] <- 0
    ub[idx_r[!selected]] <- 0
  }
  if (!is.null(idx_u)) {
    lb[idx_u] <- 0
    ub[idx_u[!selected]] <- 0
  }
  if (!is.null(idx_z)) {
    lb[idx_z] <- 0
  }
  list(lb = lb, ub = ub, selected = selected)
}

.audience_compute_profiled_risk <- function(prep, a_exp, lambda,
                                            n_total = NULL) {
  if (is.null(n_total)) n_total <- prep$n_total
  a_obs <- prep$omega - a_exp
  bias_term <- (sum(prep$bias_weights * abs(a_obs)))^2
  if (lambda == 1) return(as.numeric(bias_term))
  exp_term <- (sum(abs(a_exp) * prep$v * sqrt(prep$costs))^2) / n_total
  obs_term <- as.numeric(t(a_obs) %*% prep$Sigma_obs %*% a_obs)
  .audience_scalarize(obs_term + exp_term, bias_term, lambda)
}

.audience_compute_risk_with_n <- function(prep, a_exp, n_opt, lambda,
                                          tol = 1e-10) {
  a_obs <- prep$omega - a_exp
  bias_term <- (sum(prep$bias_weights * abs(a_obs)))^2
  if (lambda == 1) return(as.numeric(bias_term))
  obs_term <- as.numeric(t(a_obs) %*% prep$Sigma_obs %*% a_obs)
  exp_terms <- numeric(prep$p)
  for (j in seq_len(prep$p)) {
    if (a_exp[j] == 0 || prep$v2[j] == 0) {
      exp_terms[j] <- 0
    } else if (n_opt[j] <= 0) {
      exp_terms[j] <- Inf
    } else {
      exp_terms[j] <- prep$v2[j] * a_exp[j]^2 / n_opt[j]
    }
  }
  .audience_scalarize(obs_term + sum(exp_terms), bias_term, lambda)
}

.audience_solve_oracle_lambda_set <- function(prep, S, lambda,
                                              gurobi_params = list()) {
  .minimax_require("gurobi")
  p <- prep$p
  nvar <- 3L * p
  idx_a <- seq_len(p)
  idx_r <- p + seq_len(p)
  idx_z <- 2L * p + seq_len(p)
  bounds <- .audience_set_bounds(prep, S, nvar, idx_a = idx_a,
                                 idx_r = idx_r, idx_z = idx_z)

  model <- list(
    modelname = "audience_oracle_lambda",
    modelsense = "min",
    vtype = rep("C", nvar),
    lb = bounds$lb,
    ub = bounds$ub,
    obj = rep(0, nvar)
  )
  one_minus_lambda <- 1 - lambda
  model$obj[idx_a] <- one_minus_lambda * prep$q_obs
  M_exp <- ((prep$v * sqrt(prep$costs)) %o% (prep$v * sqrt(prep$costs))) /
    prep$n_total
  model$Q <- .block_quadratic_sum(
    list(
      list(M = one_minus_lambda * prep$Sigma_obs, idx = idx_a),
      list(M = one_minus_lambda * M_exp, idx = idx_r),
      list(M = lambda * prep$M_bias, idx = idx_z)
    ),
    nvar
  )

  A <- .empty_sparse(0L, nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .audience_abs_constraints(A, rhs, sense, nvar, idx_a, idx_z,
                                   prep$omega)
  tmp <- .audience_r_abs_constraints(tmp$A, tmp$rhs, tmp$sense, nvar,
                                     idx_a, idx_r)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  result <- .gurobi_solve(model, params = gurobi_params)
  result_status <- .gurobi_status(result)
  if (is.null(result$x) || !identical(result_status, "OPTIMAL")) {
    return(list(error = result_status, result = result))
  }
  a_exp <- as.numeric(result$x[idx_a])
  a_exp[!bounds$selected] <- 0
  risk_scale <- (1 - lambda) * prep$alpha_scale +
    lambda * prep$beta_scale
  value <- .safe_nonnegative(
    .audience_compute_profiled_risk(prep, a_exp, lambda),
    scale = risk_scale, tol = prep$zero_tol, name = "oracle risk"
  )
  list(
    value = value,
    a_exp = a_exp,
    a_obs = prep$omega - a_exp,
    x = as.integer(bounds$selected),
    result = result
  )
}

.audience_solve_oracle_lambda_set_explicit <- function(
    prep, S, lambda, gurobi_params = list()) {
  .minimax_require("gurobi")
  p <- prep$p
  selected <- rep(FALSE, p)
  selected[S] <- TRUE
  if (sum(prep$costs[selected] * prep$n_min) >
      prep$n_total * (1 + 1e-10)) {
    return(list(error = "N_MIN_EXCEEDS_BUDGET"))
  }

  idx_n <- seq_len(p)
  idx_a <- p + seq_len(p)
  idx_z <- 2L * p + seq_len(p)
  idx_u <- 3L * p + seq_len(p)
  nvar <- 4L * p

  lb <- rep(-Inf, nvar)
  ub <- rep(Inf, nvar)
  lb[idx_n] <- 0
  ub[idx_n] <- 0
  lb[idx_n[selected]] <- prep$n_min / prep$n_total
  ub[idx_n[selected]] <- 1 / prep$costs[selected]
  lb[idx_a[!selected]] <- 0
  ub[idx_a[!selected]] <- 0
  lb[idx_z] <- 0
  lb[idx_u] <- 0
  ub[idx_u[!selected]] <- 0

  A <- .empty_sparse(0L, nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .audience_add_budget_constraint(
    A, rhs, sense, nvar, idx_n, prep, selected
  )
  tmp <- .audience_abs_constraints(
    tmp$A, tmp$rhs, tmp$sense, nvar, idx_a, idx_z, prep$omega
  )

  one_minus_lambda <- 1 - lambda
  obj <- rep(0, nvar)
  obj[idx_a[selected]] <- one_minus_lambda * prep$q_obs[selected]
  obj[idx_u] <- one_minus_lambda / prep$n_total
  Q <- .block_quadratic_sum(
    list(
      list(M = one_minus_lambda *
             prep$Sigma_obs[selected, selected, drop = FALSE],
           idx = idx_a[selected]),
      list(M = lambda * prep$M_bias, idx = idx_z)
    ),
    nvar
  )
  quadcon <- list()
  if (lambda < 1) {
    quadcon <- .audience_add_rotated_constraints(
      quadcon, nvar, idx_a, idx_u, idx_n, prep, selected,
      "oracle", u_scale = 1
    )
  }

  model <- list(
    modelname = "audience_oracle_lambda_with_floor",
    modelsense = "min",
    vtype = rep("C", nvar),
    lb = lb,
    ub = ub,
    obj = obj,
    objcon = one_minus_lambda * prep$const_obs,
    Q = Q,
    A = tmp$A,
    rhs = tmp$rhs,
    sense = tmp$sense,
    quadcon = quadcon
  )

  result <- tryCatch(
    .gurobi_solve(model, params = gurobi_params),
    error = function(e) structure(
      list(error = conditionMessage(e)), class = "audience_error"
    )
  )
  if (inherits(result, "audience_error")) return(list(error = result$error))
  status <- .gurobi_status(result)
  if (is.null(result$x) || !identical(status, "OPTIMAL")) {
    return(list(error = status, result = result))
  }

  a_exp <- as.numeric(result$x[idx_a])
  a_exp[!selected] <- 0
  n_opt <- rep(0, p)
  n_opt[selected] <- prep$n_total *
    pmax(0, as.numeric(result$x[idx_n[selected]]))
  n_opt <- .repair_explicit_allocation(
    prep, n_opt, x = as.integer(selected), n_min = prep$n_min, tol = 1e-8
  )
  risk_scale <- (1 - lambda) * prep$alpha_scale +
    lambda * prep$beta_scale
  value <- .safe_nonnegative(
    .audience_compute_risk_with_n(prep, a_exp, n_opt, lambda),
    scale = risk_scale, tol = prep$zero_tol, name = "oracle risk"
  )
  list(
    value = value,
    a_exp = a_exp,
    a_obs = prep$omega - a_exp,
    n_opt = n_opt,
    x = as.integer(selected),
    result = result,
    status = status
  )
}

.audience_solve_beta_set <- function(prep, S, gurobi_params = list()) {
  .minimax_require("gurobi")
  p <- prep$p
  nvar <- 2L * p
  idx_a <- seq_len(p)
  idx_z <- p + seq_len(p)
  bounds <- .audience_set_bounds(prep, S, nvar, idx_a = idx_a, idx_z = idx_z)

  model <- list(
    modelname = "audience_beta_star",
    modelsense = "min",
    vtype = rep("C", nvar),
    lb = bounds$lb,
    ub = bounds$ub,
    obj = rep(0, nvar),
    Q = .block_quadratic(prep$M_bias, idx_z, nvar)
  )

  tmp <- .audience_abs_constraints(.empty_sparse(0L, nvar), numeric(0),
                                   character(0), nvar, idx_a, idx_z,
                                   prep$omega)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense
  result <- .gurobi_solve(model, params = gurobi_params)
  result_status <- .gurobi_status(result)
  if (is.null(result$x) || !identical(result_status, "OPTIMAL")) {
    return(list(error = result_status, result = result))
  }
  a_exp <- as.numeric(result$x[idx_a])
  a_exp[!bounds$selected] <- 0
  value <- .safe_nonnegative(
    (sum(prep$bias_weights * abs(prep$omega - a_exp)))^2,
    scale = prep$beta_scale, tol = prep$zero_tol, name = "beta oracle"
  )
  list(
    value = value,
    a_exp = a_exp,
    a_obs = prep$omega - a_exp,
    x = as.integer(bounds$selected),
    result = result
  )
}


.audience_beta_zero_feasible <- function(prep, x) {
  x <- as.numeric(x)
  if (length(x) != prep$p || any(!is.finite(x))) {
    stop("x must be a finite vector with one entry per experiment.",
         call. = FALSE)
  }
  selected <- x > 0.5
  relevant <- prep$bias_weights > 0
  all(!relevant | prep$omega == 0 | selected)
}

.audience_compute_oracles <- function(prep, gurobi_params = list()) {
  beta_best <- NULL
  for (S in prep$sets) {
    res <- .audience_solve_beta_set(prep, S, gurobi_params)
    if (!is.null(res$error)) next
    if (is.null(beta_best) || res$value < beta_best$value) {
      beta_best <- base::c(res, list(set = S))
    }
  }
  if (is.null(beta_best) || !is.finite(beta_best$value)) {
    stop("No feasible audience-regret beta oracle found.", call. = FALSE)
  }
  beta_best$value <- .safe_nonnegative(
    beta_best$value, scale = prep$beta_scale,
    tol = prep$zero_tol, name = "beta_star"
  )
  beta_zero <- .audience_beta_zero_feasible(prep, beta_best$x)
  if (beta_zero) {
    relevant <- prep$bias_weights > 0
    beta_best$a_exp[relevant] <- prep$omega[relevant]
    beta_best$a_obs <- prep$omega - beta_best$a_exp
    beta_best$value <- 0
  }

  G <- length(prep$lambda_grid)
  F_star <- rep(Inf, G)
  F_zero <- rep(FALSE, G)
  oracle_by_lambda <- vector("list", G)
  for (g in seq_len(G)) {
    lambda <- prep$lambda_grid[g]
    if (lambda == 1) {
      F_star[g] <- beta_best$value
      F_zero[g] <- beta_zero
      oracle_by_lambda[[g]] <- beta_best
      next
    }
    best <- NULL
    for (S in prep$sets) {
      res <- if (prep$n_min > 0) {
        .audience_solve_oracle_lambda_set_explicit(
          prep, S, lambda, gurobi_params
        )
      } else {
        .audience_solve_oracle_lambda_set(
          prep, S, lambda, gurobi_params
        )
      }
      if (!is.null(res$error)) next
      if (is.null(best) || res$value < best$value) {
        best <- base::c(res, list(set = S))
      }
    }
    if (is.null(best) || !is.finite(best$value)) {
      stop(sprintf(
        "No feasible audience-regret oracle found for lambda_grid[%d] = %g.",
        g, lambda
      ), call. = FALSE)
    }
    best$value <- .safe_nonnegative(
      best$value, scale = prep$risk_scale[g],
      tol = prep$zero_tol, name = "F_star"
    )
    # For lambda < 1, do not turn a genuinely small positive oracle risk into
    # an exact zero. Target normalization already removes target-unit scaling;
    # exact zero is reserved for an exactly zero validated objective.
    F_zero[g] <- identical(as.numeric(best$value), 0)
    F_star[g] <- best$value
    oracle_by_lambda[[g]] <- best
  }

  list(
    F_star = F_star,
    F_zero = F_zero,
    beta_star = beta_best$value,
    beta_zero = beta_zero,
    oracle_by_lambda = oracle_by_lambda,
    beta_oracle = beta_best
  )
}

.audience_add_budget_constraint <- function(A, rhs, sense, nvar, idx_n,
                                            prep, selected) {
  row <- .empty_sparse(1L, nvar)
  row[1L, idx_n[selected]] <- prep$costs[selected]
  A <- rbind(A, row)
  rhs <- base::c(rhs, 1)
  sense <- base::c(sense, "=")
  list(A = A, rhs = rhs, sense = sense)
}

.audience_add_rotated_constraints <- function(quadcon, nvar, idx_a, idx_u,
                                              idx_n, prep, selected, g_label,
                                              u_scale = 1) {
  for (j in which(selected)) {
    Qc <- .empty_sparse(nvar, nvar)
    if (prep$v2[j] > 0) {
      Qc[idx_a[j], idx_a[j]] <- u_scale * prep$v2[j]
    }
    Qc[idx_u[j], idx_n[j]] <- -0.5
    Qc[idx_n[j], idx_u[j]] <- -0.5
    quadcon[[length(quadcon) + 1L]] <- list(
      Qc = Qc,
      q = rep(0, nvar),
      rhs = 0,
      sense = "<",
      name = paste0("exp_var_", g_label, "_", j)
    )
  }
  quadcon
}

.audience_solve_grid_set <- function(prep, S, F_star,
                                     n_min = prep$n_min,
                                     gurobi_params = list(),
                                     tol = 1e-8) {
  .minimax_require("gurobi")
  p <- prep$p
  G <- length(prep$lambda_grid)
  selected <- rep(FALSE, p)
  selected[S] <- TRUE
  n_min <- as.numeric(n_min)
  if (!is.numeric(n_min) || length(n_min) != 1L || !is.finite(n_min) ||
      n_min < 0) {
    stop("n_min must be a nonnegative finite scalar.", call. = FALSE)
  }
  if (sum(prep$costs[selected] * n_min) >
      prep$n_total * (1 + 1e-10)) {
    return(list(error = "N_MIN_EXCEEDS_BUDGET"))
  }

  offset <- 0L
  idx_t <- 1L
  offset <- 1L
  idx_n <- offset + seq_len(p)
  offset <- offset + p
  idx_a <- idx_z <- idx_u <- vector("list", G)
  for (g in seq_len(G)) {
    idx_a[[g]] <- offset + seq_len(p)
    offset <- offset + p
    idx_z[[g]] <- offset + seq_len(p)
    offset <- offset + p
    idx_u[[g]] <- offset + seq_len(p)
    offset <- offset + p
  }
  nvar <- offset

  lb <- rep(-Inf, nvar)
  ub <- rep(Inf, nvar)
  lb[idx_t] <- 1
  lb[idx_n] <- 0
  ub[idx_n] <- 0
  lb[idx_n[selected]] <- n_min / prep$n_total
  ub[idx_n[selected]] <- 1 / prep$costs[selected]
  for (g in seq_len(G)) {
    lb[idx_a[[g]][!selected]] <- 0
    ub[idx_a[[g]][!selected]] <- 0
    lb[idx_z[[g]]] <- 0
    lb[idx_u[[g]]] <- 0
    ub[idx_u[[g]][!selected]] <- 0
  }

  A <- .empty_sparse(0L, nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .audience_add_budget_constraint(A, rhs, sense, nvar, idx_n,
                                         prep, selected)
  A <- tmp$A; rhs <- tmp$rhs; sense <- tmp$sense
  for (g in seq_len(G)) {
    tmp <- .audience_abs_constraints(A, rhs, sense, nvar, idx_a[[g]],
                                     idx_z[[g]], prep$omega)
    A <- tmp$A; rhs <- tmp$rhs; sense <- tmp$sense
  }

  quadcon <- list()
  ## Store sample shares in idx_n and n_total times each experimental-variance
  ## epigraph in idx_u. Then v2[j] * a[j]^2 <= idx_u[j] * idx_n[j],
  ## keeping both sides well scaled while preserving the original constraint.
  u_scale <- prep$n_total
  for (g in seq_len(G)) {
    lambda <- prep$lambda_grid[g]
    one_minus_lambda <- 1 - lambda
    q <- rep(0, nvar)
    q[idx_a[[g]][selected]] <-
      one_minus_lambda * prep$q_obs[selected]
    q[idx_u[[g]]] <- one_minus_lambda / u_scale
    q[idx_t] <- -F_star[g]
    rhs_g <- -one_minus_lambda * prep$const_obs
    scale_g <- .quad_scale_multi(
      list(one_minus_lambda *
             prep$Sigma_obs[selected, selected, drop = FALSE],
           lambda * prep$M_bias),
      q, rhs_g
    )
    quadcon[[length(quadcon) + 1L]] <- list(
      Qc = .block_quadratic_sum(
        list(
          list(M = one_minus_lambda *
                 prep$Sigma_obs[selected, selected, drop = FALSE],
               idx = idx_a[[g]][selected], scale = scale_g),
          list(M = lambda * prep$M_bias, idx = idx_z[[g]],
               scale = scale_g)
        ),
        nvar
      ),
      q = scale_g * q,
      rhs = scale_g * rhs_g,
      sense = "<",
      name = paste0("risk_lambda_", g)
    )
    if (lambda < 1) {
      quadcon <- .audience_add_rotated_constraints(
        quadcon, nvar, idx_a[[g]], idx_u[[g]], idx_n, prep, selected,
        paste0("lambda", g), u_scale = 1
      )
    }
  }

  model <- list(
    modelname = "audience_grid_set",
    modelsense = "min",
    vtype = rep("C", nvar),
    lb = lb,
    ub = ub,
    obj = replace(rep(0, nvar), idx_t, 1),
    A = A,
    rhs = rhs,
    sense = sense,
    quadcon = quadcon
  )

  result <- tryCatch(
    .gurobi_solve(model, params = gurobi_params),
    error = function(e) structure(list(error = conditionMessage(e)), class = "audience_error")
  )
  if (inherits(result, "audience_error")) {
    return(list(error = result$error))
  }
  result_status <- .gurobi_status(result)
  if (is.null(result$x) || !identical(result_status, "OPTIMAL")) {
    return(list(error = result_status, result = result))
  }

  sol <- result$x
  n_opt <- rep(0, p)
  n_opt[selected] <- prep$n_total *
    pmax(0, as.numeric(sol[idx_n[selected]]))
  n_opt <- .repair_explicit_allocation(
    prep, n_opt, x = as.integer(selected), n_min = n_min, tol = 1e-8
  )
  a_exp_by_lambda <- matrix(0, nrow = G, ncol = p)
  a_obs_by_lambda <- matrix(0, nrow = G, ncol = p)
  risk <- ratio <- numeric(G)
  for (g in seq_len(G)) {
    a_exp <- as.numeric(sol[idx_a[[g]]])
    a_exp[!selected] <- 0
    a_exp_by_lambda[g, ] <- a_exp
    a_obs_by_lambda[g, ] <- prep$omega - a_exp
    risk[g] <- .safe_nonnegative(
      .audience_compute_risk_with_n(
        prep, a_exp, n_opt, prep$lambda_grid[g], tol = tol
      ),
      scale = prep$risk_scale[g], tol = prep$zero_tol, name = "risk"
    )
    if (F_star[g] == 0 &&
        .relative_zero(risk[g], prep$risk_scale[g], prep$zero_tol)) {
      risk[g] <- 0
    }
    ratio[g] <- .safe_ratio(risk[g], F_star[g], tol = prep$zero_tol)
  }
  solver_r <- as.numeric(sol[idx_t])
  recomputed_r <- max(ratio)
  epigraph <- tryCatch(
    .validate_regret_epigraph(
      recomputed_r, solver_r, tol = max(1e-6, 10 * tol)
    ),
    error = function(e) structure(
      list(error = conditionMessage(e)), class = "audience_error"
    )
  )
  if (inherits(epigraph, "audience_error")) {
    return(list(error = paste0("EPIGRAPH_VALIDATION_FAILED: ",
                               epigraph$error), result = result))
  }
  r_validated <- epigraph$t
  result$validated_ratio_gap <- epigraph$validated_gap

  colnames(a_exp_by_lambda) <- colnames(a_obs_by_lambda) <- prep$arm_names
  grid_names <- paste0("lambda_", seq_len(G) - 1L)
  rownames(a_exp_by_lambda) <- rownames(a_obs_by_lambda) <- grid_names
  list(
    set = S,
    x = as.integer(selected),
    n_opt = n_opt,
    n_share = if (sum(n_opt) > 0) n_opt / sum(n_opt) else rep(0, p),
    r = r_validated,
    solver_r = solver_r,
    risk = risk,
    ratio = ratio,
    a_exp_by_lambda = a_exp_by_lambda,
    a_obs_by_lambda = a_obs_by_lambda,
    result = result,
    status = result_status
  )
}

.audience_rescale_solution <- function(solution, scale) {
  if (is.null(solution) || !is.list(solution)) return(solution)
  out <- solution
  for (name in base::c("value", "risk", "oracle")) {
    if (!is.null(out[[name]])) out[[name]] <- out[[name]] * scale^2
  }
  for (name in base::c("a_exp", "a_obs", "a_exp_by_lambda", "a_obs_by_lambda")) {
    if (!is.null(out[[name]])) out[[name]] <- out[[name]] * scale
  }
  out
}

#' Solve the audience-regret experimental design problem
#'
#' @description
#' Implements the Appendix C.2 profiled posterior/audience-regret optimization
#' on a researcher-supplied grid of scalarization weights `lambda_grid`. The
#' default is 50 equally spaced points on `[0, 1]`, including both endpoints.
#' The current implementation assumes the identity loading case, so
#' `a_obs + a_exp = omega`, enumerates feasible experiment sets, and solves the
#' fixed-set convex quadratic-conic problems with Gurobi.
#'
#' @inheritParams solve_minimax_design
#' @param lambda_grid Numeric vector of scalarization weights in `[0, 1]` that
#'   includes both endpoints. Defaults to `seq(0, 1, length.out = 50L)`.
#' @param n_min Optional nonnegative lower bound on sample size for selected
#'   experiments. The same floor is imposed on candidate designs and every
#'   oracle denominator `F_star(lambda)`.
#' @param zero_tol Relative tolerance for recognizing exact-zero oracle risks.
#' @param psd_tol Relative covariance eigenvalue tolerance.
#' @param psd_action Either `"error"` or `"project"`. Corrections beyond
#'   floating-point roundoff require explicit `"project"`; projections are
#'   reported in the fitted object.
#' @param tol Numerical tolerance for the dimensionless lambda grid and solver
#'   diagnostics.
#'
#' @return A list of class `audience_regret_design`.
#' @export
solve_audience_regret_design <- function(Sigma_obs,
                                         v2,
                                         n_total,
                                         omega,
                                         lambda_grid = seq(0, 1,
                                                           length.out = 50L),
                                         costs = NULL,
                                         bias_weights = NULL,
                                         h = length(omega),
                                         x_max = NULL,
                                         feasible_sets = NULL,
                                         selection_constraints = NULL,
                                         min_experiments = 1L,
                                         n_min = 0,
                                         gurobi_params = list(),
                                         max_sets = 200000L,
                                         zero_tol = 1e-9,
                                         psd_tol = 1e-10,
                                         psd_action = base::c("error", "project"),
                                         tol = 1e-8,
                                         c = NULL,
                                         k = NULL) {
  .minimax_require("gurobi")
  lambda_grid_supplied <- !missing(lambda_grid)
  prep <- .audience_prepare(
    Sigma_obs = Sigma_obs,
    v2 = v2,
    n_total = n_total,
    omega = omega,
    lambda_grid = lambda_grid,
    costs = costs,
    bias_weights = bias_weights,
    h = h,
    x_max = x_max,
    feasible_sets = feasible_sets,
    selection_constraints = selection_constraints,
    min_experiments = min_experiments,
    n_min = n_min,
    costs_alias = c,
    bias_weights_alias = k,
    lambda_grid_supplied = lambda_grid_supplied,
    max_sets = max_sets,
    zero_tol = zero_tol,
    psd_tol = psd_tol,
    psd_action = psd_action,
    tol = tol
  )

  oracles <- .audience_compute_oracles(prep, gurobi_params = gurobi_params)
  F_star <- oracles$F_star
  beta_star <- oracles$beta_star
  if (any(!is.finite(F_star)) || any(F_star < 0)) {
    stop("Audience-regret oracle values F_star(lambda) must be finite and nonnegative.",
         call. = FALSE)
  }

  candidate_solutions <- vector("list", length(prep$sets))
  diagnostics <- data.frame(
    candidate = seq_along(prep$sets),
    selected = vapply(prep$sets, function(S) paste(prep$arm_names[S],
                                                   collapse = ", "),
                      character(1)),
    status = NA_character_,
    objective = NA_real_,
    max_grid_ratio = NA_real_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(prep$sets)) {
    res <- .audience_solve_grid_set(
      prep, prep$sets[[i]], F_star = F_star,
      n_min = prep$n_min, gurobi_params = gurobi_params, tol = tol
    )
    candidate_solutions[[i]] <- res
    if (!is.null(res$error)) {
      diagnostics$status[i] <- res$error
    } else {
      diagnostics$status[i] <- res$status
      diagnostics$objective[i] <- res$r
      diagnostics$max_grid_ratio[i] <- max(res$ratio)
    }
  }
  feasible <- which(diagnostics$status == "OPTIMAL" &
                      is.finite(diagnostics$objective))
  if (length(feasible) == 0L) {
    stop("No feasible audience-regret grid solution was found.", call. = FALSE)
  }
  best_i <- feasible[which.min(diagnostics$objective[feasible])]
  best <- candidate_solutions[[best_i]]

  risk_by_lambda <- data.frame(
    grid_id = seq_along(prep$lambda_grid) - 1L,
    lambda = prep$lambda_grid,
    risk = best$risk,
    oracle = F_star,
    ratio = best$ratio,
    stringsAsFactors = FALSE
  )

  names_out <- prep$arm_names
  names(best$x) <- names(best$n_opt) <- names(best$n_share) <- names_out

  risk_by_lambda$risk <- risk_by_lambda$risk * prep$risk_rescale
  risk_by_lambda$oracle <- risk_by_lambda$oracle * prep$risk_rescale
  F_star_natural <- F_star * prep$risk_rescale
  beta_star_natural <- beta_star * prep$risk_rescale
  best_a_exp_natural <- best$a_exp_by_lambda * prep$target_scale
  best_a_obs_natural <- best$a_obs_by_lambda * prep$target_scale
  oracle_by_lambda_out <- lapply(
    oracles$oracle_by_lambda, .audience_rescale_solution,
    scale = prep$target_scale
  )
  beta_oracle_out <- .audience_rescale_solution(
    oracles$beta_oracle, prep$target_scale
  )
  candidate_solutions_out <- lapply(
    candidate_solutions, .audience_rescale_solution,
    scale = prep$target_scale
  )

  out <- list(
    selected = prep$arm_names[which(best$x == 1L)],
    x_opt = best$x,
    n_opt = best$n_opt,
    n_share = best$n_share,
    budget_used = .make_named(best$n_opt * prep$costs, names_out),
    budget_share = .make_named(
      if (sum(best$n_opt * prep$costs) > 0) {
        best$n_opt * prep$costs / sum(best$n_opt * prep$costs)
      } else {
        rep(0, prep$p)
      },
      names_out
    ),
    r_opt = best$r,
    solver_r_opt = best$solver_r,
    lambda_grid = prep$lambda_grid,
    lambda_grid_diagnostics = prep$lambda_grid_diagnostics,
    F_star = F_star_natural,
    F_star_zero = oracles$F_zero,
    beta_star = beta_star_natural,
    beta_zero = oracles$beta_zero,
    risk_by_lambda = risk_by_lambda,
    a_exp_by_lambda = best_a_exp_natural,
    a_obs_by_lambda = best_a_obs_natural,
    oracle_by_lambda = oracle_by_lambda_out,
    beta_oracle = beta_oracle_out,
    candidate_diagnostics = diagnostics,
    candidate_solutions = candidate_solutions_out,
    optimizer_result = best$result,
    psd_diagnostics = prep$psd_diagnostics,
    inputs = list(
      Sigma_obs = prep$Sigma_obs,
      v2 = prep$v2,
      n_total = prep$n_total,
      omega = .make_named(prep$omega_original, names_out),
      costs = .make_named(prep$costs, names_out),
      bias_weights = .make_named(prep$bias_weights, names_out),
      lambda_grid = prep$lambda_grid,
      lambda_grid_diagnostics = prep$lambda_grid_diagnostics,
      h = prep$h,
      x_max = .make_named(prep$x_max, names_out),
      min_experiments = prep$min_experiments,
      n_min = prep$n_min,
      feasible_sets = prep$sets,
      zero_tol = prep$zero_tol,
      psd_tol = prep$psd_tol,
      psd_action = prep$psd_action,
      psd_diagnostics = prep$psd_diagnostics,
      loading_structure = prep$loading_structure,
      internal_target_scale = prep$target_scale
    )
  )
  class(out) <- base::c("audience_regret_design", "list")
  out
}

#' @export
print.audience_regret_design <- function(x, ...) {
  cat("Audience-regret design\n")
  cat("Selected experiments: ", .format_selected(x$x_opt, names(x$x_opt)), "\n", sep = "")
  cat("Audience regret: ", format(x$r_opt, digits = 6), "\n", sep = "")
  cat("Max grid ratio: ", format(max(x$risk_by_lambda$ratio), digits = 6), "\n", sep = "")
  cat("Lambda grid: ", length(x$lambda_grid), " points in [",
      format(min(x$lambda_grid), digits = 4), ", ",
      format(max(x$lambda_grid), digits = 4), "] (",
      x$lambda_grid_diagnostics$source, ")\n", sep = "")
  invisible(x)
}
