.minimax_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required.", package), call. = FALSE)
  }
}

.normalize_sense <- function(sense) {
  if (is.null(sense)) {
    return(character(0))
  }
  out <- as.character(sense)
  out[out %in% c("<=", "=<")] <- "<"
  out[out %in% c(">=", "=>")] <- ">"
  if (!all(out %in% c("<", ">", "="))) {
    stop("Constraint senses must be one of '<=', '>=', '<', '>', or '='.", call. = FALSE)
  }
  out
}

.is_integerish_scalar <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && abs(x - round(x)) < 1e-9
}

.numeric_vector <- function(x, name, p = NULL, default = NULL) {
  if (is.null(x)) {
    x <- default
  }
  if (is.null(x)) {
    stop(sprintf("%s must be supplied.", name), call. = FALSE)
  }
  x <- as.numeric(x)
  if (!is.null(p) && length(x) != p) {
    stop(sprintf("%s must have length %d.", name, p), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("%s must contain only finite values.", name), call. = FALSE)
  }
  x
}

.validate_selection_constraints <- function(selection_constraints, p) {
  if (is.null(selection_constraints)) {
    return(NULL)
  }
  if (!is.list(selection_constraints)) {
    stop("selection_constraints must be a list with elements A, sense, and rhs.", call. = FALSE)
  }
  A <- selection_constraints$A
  sense <- selection_constraints$sense
  rhs <- selection_constraints$rhs
  if (is.null(A) || is.null(sense) || is.null(rhs)) {
    stop("selection_constraints must contain A, sense, and rhs.", call. = FALSE)
  }
  A <- as.matrix(A)
  if (ncol(A) != p) {
    stop("selection_constraints$A must have one column per experiment.", call. = FALSE)
  }
  sense <- .normalize_sense(sense)
  rhs <- as.numeric(rhs)
  if (nrow(A) != length(sense) || length(sense) != length(rhs)) {
    stop("selection_constraints rows, senses, and right-hand sides must have matching lengths.", call. = FALSE)
  }
  if (any(!is.finite(A)) || any(!is.finite(rhs))) {
    stop("selection_constraints must contain finite values.", call. = FALSE)
  }
  list(A = A, sense = sense, rhs = rhs)
}

.satisfies_constraint <- function(lhs, sense, rhs, tol = 1e-9) {
  if (sense == "<") {
    return(lhs <= rhs + tol)
  }
  if (sense == ">") {
    return(lhs >= rhs - tol)
  }
  abs(lhs - rhs) <= tol
}

.satisfies_x_constraints <- function(x, selection_constraints, tol = 1e-9) {
  if (is.null(selection_constraints)) {
    return(TRUE)
  }
  lhs <- as.numeric(selection_constraints$A %*% as.numeric(x))
  all(mapply(.satisfies_constraint, lhs, selection_constraints$sense,
             selection_constraints$rhs, MoreArgs = list(tol = tol)))
}

.process_feasible_sets <- function(feasible_sets, p, h, x_max, min_experiments,
                                   selection_constraints) {
  if (is.null(feasible_sets)) {
    return(list(use_sets = FALSE, sets = NULL, Aset = NULL, Msets = 0L))
  }
  if (!is.list(feasible_sets) || length(feasible_sets) == 0L) {
    stop("feasible_sets must be NULL or a non-empty list of integer vectors.", call. = FALSE)
  }

  sets <- lapply(feasible_sets, function(S) {
    S <- unique(as.integer(S))
    if (anyNA(S)) {
      stop("feasible_sets contains NA.", call. = FALSE)
    }
    if (any(S < 1L | S > p)) {
      stop("feasible_sets contains indices outside 1:length(omega).", call. = FALSE)
    }
    sort(S)
  })

  sets <- Filter(function(S) length(S) <= h, sets)
  sets <- Filter(function(S) length(S) >= min_experiments, sets)
  sets <- Filter(function(S) all(x_max[S] > 0.5), sets)
  sets <- Filter(function(S) {
    x <- rep(0, p)
    if (length(S) > 0L) {
      x[S] <- 1
    }
    .satisfies_x_constraints(x, selection_constraints)
  }, sets)

  key <- vapply(sets, function(S) {
    if (length(S) == 0L) "<empty>" else paste(S, collapse = ",")
  }, character(1))
  sets <- sets[!duplicated(key)]

  if (length(sets) == 0L) {
    stop("No feasible experiment set remains after applying the supplied restrictions.", call. = FALSE)
  }

  Aset <- matrix(0, nrow = length(sets), ncol = p)
  for (m in seq_along(sets)) {
    if (length(sets[[m]]) > 0L) {
      Aset[m, sets[[m]]] <- 1
    }
  }
  list(use_sets = TRUE, sets = sets, Aset = Aset, Msets = length(sets))
}

.empty_sparse <- function(nrow, ncol) {
  Matrix::Matrix(0, nrow = nrow, ncol = ncol, sparse = TRUE)
}

.add_row <- function(A_mat, rhs_vec, sense_vec, nvar, index, value, rhs, sense) {
  row <- .empty_sparse(1L, nvar)
  row[1L, index] <- value
  list(
    A = rbind(A_mat, row),
    rhs = c(rhs_vec, rhs),
    sense = c(sense_vec, .normalize_sense(sense))
  )
}

.add_x_constraints <- function(A_mat, rhs_vec, sense_vec, nvar, idx_x, prep) {
  p <- prep$p

  row <- .empty_sparse(1L, nvar)
  row[1L, idx_x] <- 1
  tmp <- list(A = rbind(A_mat, row), rhs = c(rhs_vec, prep$h),
              sense = c(sense_vec, "<"))
  A_mat <- tmp$A
  rhs_vec <- tmp$rhs
  sense_vec <- tmp$sense

  if (prep$min_experiments > 0) {
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_x] <- 1
    tmp <- list(A = rbind(A_mat, row), rhs = c(rhs_vec, prep$min_experiments),
                sense = c(sense_vec, ">"))
    A_mat <- tmp$A
    rhs_vec <- tmp$rhs
    sense_vec <- tmp$sense
  }

  if (!is.null(prep$selection_constraints)) {
    for (r in seq_len(nrow(prep$selection_constraints$A))) {
      row <- .empty_sparse(1L, nvar)
      row[1L, idx_x] <- prep$selection_constraints$A[r, ]
      A_mat <- rbind(A_mat, row)
      rhs_vec <- c(rhs_vec, prep$selection_constraints$rhs[r])
      sense_vec <- c(sense_vec, prep$selection_constraints$sense[r])
    }
  }

  if (prep$use_sets) {
    idx_y <- prep$current_idx_y
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_y] <- 1
    A_mat <- rbind(A_mat, row)
    rhs_vec <- c(rhs_vec, 1)
    sense_vec <- c(sense_vec, "=")

    for (j in seq_len(p)) {
      row <- .empty_sparse(1L, nvar)
      row[1L, idx_x[j]] <- 1
      row[1L, idx_y] <- -prep$Aset[, j]
      A_mat <- rbind(A_mat, row)
      rhs_vec <- c(rhs_vec, 0)
      sense_vec <- c(sense_vec, "=")
    }
  }

  list(A = A_mat, rhs = rhs_vec, sense = sense_vec)
}

.add_s_constraints <- function(A_mat, rhs_vec, sense_vec, nvar, idx_x, idx_s, prep) {
  p <- prep$p
  for (j in seq_len(p)) {
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_s[j]] <- 1
    row[1L, idx_x[j]] <- -prep$gamma_upper[j]
    A_mat <- rbind(A_mat, row)
    rhs_vec <- c(rhs_vec, 0)
    sense_vec <- c(sense_vec, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_s[j]] <- -1
    row[1L, idx_x[j]] <- prep$gamma_lower[j]
    A_mat <- rbind(A_mat, row)
    rhs_vec <- c(rhs_vec, 0)
    sense_vec <- c(sense_vec, "<")
  }
  list(A = A_mat, rhs = rhs_vec, sense = sense_vec)
}

.block_quadratic <- function(M, idx, nvar, scale = 1) {
  nz <- which(M != 0, arr.ind = TRUE)
  if (nrow(nz) == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(nvar, nvar)))
  }
  Matrix::sparseMatrix(
    i = idx[nz[, 1]],
    j = idx[nz[, 2]],
    x = as.numeric(scale * M[nz]),
    dims = c(nvar, nvar)
  )
}

.quad_scale <- function(M, q, rhs) {
  scale_base <- max(c(abs(as.numeric(M)), abs(as.numeric(q)), abs(rhs), 1e-12), na.rm = TRUE)
  min(1e8, max(1, 1 / scale_base))
}

.safe_ratio <- function(value, reference, tol = 1e-8) {
  if (abs(reference) <= tol) {
    if (abs(value) <= sqrt(tol)) {
      return(1)
    }
    return(Inf)
  }
  as.numeric(value / reference)
}

.prepare_minimax_problem <- function(Sigma_obs, v2, n_total, omega,
                                     costs = NULL, bias_weights = NULL,
                                     h = NULL, x_max = NULL,
                                     feasible_sets = NULL,
                                     selection_constraints = NULL,
                                     min_experiments = 0L,
                                     gamma_lower = NULL,
                                     gamma_upper = NULL,
                                     costs_alias = NULL,
                                     bias_weights_alias = NULL,
                                     kappa_alias = NULL) {
  .minimax_require("Matrix")

  if (!is.null(costs_alias)) {
    if (!is.null(costs) && !isTRUE(all.equal(as.numeric(costs), as.numeric(costs_alias)))) {
      stop("Use either costs or c, or supply identical values for both.", call. = FALSE)
    }
    costs <- costs_alias
  }
  bias_aliases <- list(k = bias_weights_alias, kappa = kappa_alias)
  for (alias_name in names(bias_aliases)) {
    alias_value <- bias_aliases[[alias_name]]
    if (!is.null(alias_value)) {
      if (!is.null(bias_weights) && !isTRUE(all.equal(as.numeric(bias_weights), as.numeric(alias_value)))) {
        stop("Use only one of bias_weights, k, or kappa, or supply identical values.",
             call. = FALSE)
      }
      bias_weights <- alias_value
    }
  }

  omega <- .numeric_vector(omega, "omega")
  p <- length(omega)
  if (p < 1L) {
    stop("omega must have positive length.", call. = FALSE)
  }

  Sigma_obs <- as.matrix(Sigma_obs)
  if (!is.numeric(Sigma_obs) || any(dim(Sigma_obs) != c(p, p))) {
    stop("Sigma_obs must be a p by p numeric matrix with p = length(omega).", call. = FALSE)
  }
  if (any(!is.finite(Sigma_obs))) {
    stop("Sigma_obs must contain only finite values.", call. = FALSE)
  }
  Sigma_obs <- 0.5 * (Sigma_obs + t(Sigma_obs))

  v2 <- .numeric_vector(v2, "v2", p)
  if (any(v2 < 0)) {
    stop("v2 must be nonnegative.", call. = FALSE)
  }
  if (!is.numeric(n_total) || length(n_total) != 1L || !is.finite(n_total) || n_total <= 0) {
    stop("n_total must be a positive scalar budget.", call. = FALSE)
  }
  costs <- .numeric_vector(costs, "costs", p, default = rep(1, p))
  if (any(costs <= 0)) {
    stop("costs must be strictly positive.", call. = FALSE)
  }
  bias_weights <- .numeric_vector(bias_weights, "bias_weights", p, default = rep(1, p))
  if (any(bias_weights < 0)) {
    stop("bias_weights must be nonnegative.", call. = FALSE)
  }
  x_max <- .numeric_vector(x_max, "x_max", p, default = rep(1, p))
  if (any(!(x_max %in% c(0, 1)))) {
    stop("x_max must contain only 0 or 1.", call. = FALSE)
  }

  gamma_lower <- .numeric_vector(gamma_lower, "gamma_lower", p, default = rep(0, p))
  gamma_upper <- .numeric_vector(gamma_upper, "gamma_upper", p, default = rep(1, p))
  if (any(gamma_lower < 0) || any(gamma_upper > 1) || any(gamma_lower > gamma_upper)) {
    stop("gamma_lower and gamma_upper must satisfy 0 <= lower <= upper <= 1.", call. = FALSE)
  }

  if (is.null(h)) {
    h <- p
  }
  if (!.is_integerish_scalar(h) || h < 0 || h > p) {
    stop("h must be an integer scalar between 0 and length(omega).", call. = FALSE)
  }
  h <- as.integer(round(h))
  if (!.is_integerish_scalar(min_experiments) || min_experiments < 0 || min_experiments > h) {
    stop("min_experiments must be an integer scalar between 0 and h.", call. = FALSE)
  }
  min_experiments <- as.integer(round(min_experiments))

  selection_constraints <- .validate_selection_constraints(selection_constraints, p)
  set_info <- .process_feasible_sets(feasible_sets, p, h, x_max, min_experiments,
                                     selection_constraints)

  v <- sqrt(v2)
  a <- abs(omega) * v * sqrt(costs)
  d <- bias_weights * abs(omega)
  Domega <- diag(omega, nrow = p, ncol = p)
  Qobs <- Domega %*% Sigma_obs %*% Domega
  Qobs <- 0.5 * (Qobs + t(Qobs))

  M_alpha <- (a %o% a) / n_total + Qobs
  M_alpha <- 0.5 * (M_alpha + t(M_alpha))
  q_alpha <- -2 * as.numeric(Qobs %*% rep(1, p))
  const_alpha <- as.numeric(t(rep(1, p)) %*% Qobs %*% rep(1, p))

  M_beta <- d %o% d
  q_beta <- -2 * sum(d) * d
  const_beta <- sum(d)^2

  compute_alpha <- function(s) {
    s <- as.numeric(s)
    as.numeric(t(s) %*% M_alpha %*% s + sum(q_alpha * s) + const_alpha)
  }
  compute_beta <- function(s) {
    s <- as.numeric(s)
    as.numeric((sum(d * (1 - s)))^2)
  }
  compute_allocation <- function(s) {
    s <- as.numeric(s)
    score <- abs(omega) * v * s
    denom <- sum(score * sqrt(costs))
    n_opt <- rep(0, p)
    if (denom > 0) {
      n_opt <- n_total * (score / sqrt(costs)) / denom
    }
    n_opt
  }

  arm_names <- names(omega)
  if (is.null(arm_names) || any(arm_names == "")) {
    arm_names <- paste0("arm_", seq_len(p))
  }

  c(
    list(
      p = p,
      Sigma_obs = Sigma_obs,
      v2 = v2,
      v = v,
      n_total = n_total,
      omega = omega,
      costs = costs,
      bias_weights = bias_weights,
      h = h,
      x_max = x_max,
      min_experiments = min_experiments,
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper,
      selection_constraints = selection_constraints,
      arm_names = arm_names,
      a = a,
      d = d,
      Qobs = Qobs,
      M_alpha = M_alpha,
      q_alpha = q_alpha,
      const_alpha = const_alpha,
      M_beta = M_beta,
      q_beta = q_beta,
      const_beta = const_beta,
      compute_alpha = compute_alpha,
      compute_beta = compute_beta,
      compute_allocation = compute_allocation
    ),
    set_info
  )
}

.make_named <- function(x, names) {
  names(x) <- names
  x
}

.choose_baseline_solver <- function(solver) {
  solver <- match.arg(solver, c("auto", "gurobi", "quadprog"))
  if (identical(solver, "auto")) {
    if (requireNamespace("gurobi", quietly = TRUE)) {
      return("gurobi")
    }
    return("quadprog")
  }
  solver
}

.enumerate_experiment_sets <- function(prep, max_sets = 200000L) {
  p <- prep$p
  if (prep$use_sets) {
    sets <- prep$sets
  } else {
    idx <- which(prep$x_max > 0.5)
    h_eff <- min(prep$h, length(idx))
    if (h_eff < prep$min_experiments) {
      stop("No feasible experiment set satisfies min_experiments, h, and x_max.", call. = FALSE)
    }
    sizes <- seq.int(prep$min_experiments, h_eff)
    log_terms <- vapply(sizes, function(k) lchoose(length(idx), k), numeric(1))
    approx_total <- if (length(log_terms) == 0L) 0 else {
      sum(exp(log_terms - max(log_terms))) * exp(max(log_terms))
    }
    if (!is.finite(approx_total) || approx_total > max_sets) {
      stop(sprintf(
        "The non-Gurobi solver would need to enumerate about %.0f experiment sets, above max_sets = %d. Use solver = 'gurobi', supply feasible_sets, reduce h, or increase max_sets.",
        approx_total, as.integer(max_sets)
      ), call. = FALSE)
    }
    sets <- list()
    for (k in sizes) {
      if (k == 0L) {
        sets <- c(sets, list(integer(0)))
      } else {
        sets <- c(sets, utils::combn(idx, k, simplify = FALSE))
      }
    }
    sets <- Filter(function(S) {
      x <- rep(0, p)
      if (length(S) > 0L) x[S] <- 1
      .satisfies_x_constraints(x, prep$selection_constraints)
    }, sets)
  }

  key <- vapply(sets, function(S) {
    if (length(S) == 0L) "<empty>" else paste(sort(S), collapse = ",")
  }, character(1))
  sets <- sets[!duplicated(key)]
  if (length(sets) == 0L) {
    stop("No feasible experiment set remains after applying restrictions.", call. = FALSE)
  }
  sets
}

.set_to_x <- function(S, p) {
  x <- rep(0, p)
  if (length(S) > 0L) x[S] <- 1
  x
}

.qp_min_alpha <- function(prep, lb, ub, bias_min = NULL, qp_ridge = 1e-10,
                          tol = 1e-8) {
  .minimax_require("quadprog")
  p <- prep$p
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  if (length(lb) != p || length(ub) != p) {
    stop("Internal error: lb and ub have wrong length.", call. = FALSE)
  }
  if (any(lb > ub + tol)) {
    return(NULL)
  }
  lb <- pmin(lb, ub)

  Dmat <- 2 * prep$M_alpha + diag(qp_ridge, p)
  dvec <- -prep$q_alpha

  Amat <- cbind(diag(p), -diag(p))
  bvec <- c(lb, -ub)

  if (!is.null(bias_min) && is.finite(bias_min)) {
    if (bias_min > sum(prep$d * ub) + 10 * tol) {
      return(NULL)
    }
    if (bias_min > sum(prep$d * lb) + tol) {
      Amat <- cbind(Amat, prep$d)
      bvec <- c(bvec, bias_min)
    }
  }

  out <- tryCatch(
    quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 0),
    error = function(e) NULL
  )
  if (is.null(out)) {
    return(NULL)
  }
  s <- as.numeric(out$solution)
  s <- pmax(lb, pmin(ub, s))
  if (!is.null(bias_min) && is.finite(bias_min) && sum(prep$d * s) < bias_min - 1e-6) {
    return(NULL)
  }
  s
}

.solve_alpha_oracle_enum <- function(prep, max_sets = 200000L, qp_ridge = 1e-10) {
  sets <- .enumerate_experiment_sets(prep, max_sets = max_sets)
  p <- prep$p
  best <- list(alpha = Inf, beta = Inf, s = rep(0, p), x = rep(0L, p), gamma = rep(0, p))

  for (S in sets) {
    x <- .set_to_x(S, p)
    lb <- prep$gamma_lower * x
    ub <- prep$gamma_upper * x
    s <- .qp_min_alpha(prep, lb = lb, ub = ub, qp_ridge = qp_ridge)
    if (is.null(s)) next
    alpha <- prep$compute_alpha(s)
    if (alpha < best$alpha) {
      gamma <- rep(0, p)
      gamma[x == 1L] <- s[x == 1L]
      best <- list(alpha = alpha, beta = prep$compute_beta(s), s = s,
                   x = as.integer(x), gamma = gamma)
    }
  }

  if (!is.finite(best$alpha)) {
    stop("The non-Gurobi variance oracle failed. Check feasibility and covariance inputs.",
         call. = FALSE)
  }
  best$result <- list(status = "OPTIMAL", solver = "quadprog", enumerated_sets = length(sets))
  best
}

.solve_beta_oracle_enum <- function(prep, max_sets = 200000L) {
  sets <- .enumerate_experiment_sets(prep, max_sets = max_sets)
  p <- prep$p
  best <- list(beta = Inf, alpha = Inf, s = rep(0, p), x = rep(0L, p), gamma = rep(0, p))

  for (S in sets) {
    x <- .set_to_x(S, p)
    s <- prep$gamma_upper * x
    beta <- prep$compute_beta(s)
    if (beta < best$beta) {
      gamma <- rep(0, p)
      gamma[x == 1L] <- prep$gamma_upper[x == 1L]
      best <- list(beta = beta, alpha = prep$compute_alpha(s), s = s,
                   x = as.integer(x), gamma = gamma)
    }
  }

  if (!is.finite(best$beta)) {
    stop("The non-Gurobi bias oracle failed. Check feasibility inputs.", call. = FALSE)
  }
  best$result <- list(status = "OPTIMAL", solver = "enumeration", enumerated_sets = length(sets))
  best
}

.solve_minimax_for_fixed_x_enum <- function(prep, x, alpha_star, beta_star,
                                            tol_bisect = 1e-7, bisect_iter = 80L,
                                            qp_ridge = 1e-10, tol = 1e-8) {
  p <- prep$p
  x <- as.numeric(x)
  lb <- prep$gamma_lower * x
  ub <- prep$gamma_upper * x
  total_d <- sum(prep$d)

  if (beta_star <= tol) {
    lb0 <- lb
    need_full <- prep$d > tol
    lb0[need_full] <- pmax(lb0[need_full], 1)
    if (any(ub < lb0 - tol)) {
      return(list(regret = Inf, t = Inf, alpha = Inf, beta = Inf,
                  s = rep(0, p), x = as.integer(x), gamma = rep(0, p)))
    }
    s <- .qp_min_alpha(prep, lb = lb0, ub = ub, qp_ridge = qp_ridge)
    if (is.null(s)) {
      return(list(regret = Inf, t = Inf, alpha = Inf, beta = Inf,
                  s = rep(0, p), x = as.integer(x), gamma = rep(0, p)))
    }
    alpha <- prep$compute_alpha(s)
    beta <- prep$compute_beta(s)
    gamma <- rep(0, p)
    gamma[x == 1L] <- s[x == 1L]
    alpha_ratio <- .safe_ratio(alpha, alpha_star)
    beta_ratio <- .safe_ratio(beta, beta_star)
    regret <- max(alpha_ratio, beta_ratio)
    return(list(regret = regret, t = regret, alpha = alpha, beta = beta,
                alpha_ratio = alpha_ratio, beta_ratio = beta_ratio,
                s = s, x = as.integer(x), gamma = gamma))
  }

  s_hi <- ub
  alpha_hi <- prep$compute_alpha(s_hi)
  beta_hi <- prep$compute_beta(s_hi)
  t_hi <- max(.safe_ratio(alpha_hi, alpha_star), .safe_ratio(beta_hi, beta_star), 1)
  t_lo <- 1
  best_s <- s_hi

  for (it in seq_len(bisect_iter)) {
    t_mid <- 0.5 * (t_lo + t_hi)
    bias_bound <- sqrt(max(0, t_mid * beta_star))
    bias_min <- total_d - bias_bound
    if (bias_min <= 0) {
      bias_min <- NULL
    }
    s_mid <- .qp_min_alpha(prep, lb = lb, ub = ub, bias_min = bias_min,
                           qp_ridge = qp_ridge)
    if (is.null(s_mid)) {
      t_lo <- t_mid
    } else {
      alpha_mid <- prep$compute_alpha(s_mid)
      if (alpha_mid <= t_mid * alpha_star * (1 + 1e-10)) {
        t_hi <- t_mid
        best_s <- s_mid
      } else {
        t_lo <- t_mid
      }
    }
    if ((t_hi - t_lo) <= tol_bisect) break
  }

  s <- best_s
  alpha <- prep$compute_alpha(s)
  beta <- prep$compute_beta(s)
  gamma <- rep(0, p)
  gamma[x == 1L] <- s[x == 1L]
  alpha_ratio <- .safe_ratio(alpha, alpha_star)
  beta_ratio <- .safe_ratio(beta, beta_star)
  regret <- max(alpha_ratio, beta_ratio)
  list(regret = regret, t = regret, alpha = alpha, beta = beta,
       alpha_ratio = alpha_ratio, beta_ratio = beta_ratio,
       s = s, x = as.integer(x), gamma = gamma)
}

.solve_minimax_enum <- function(prep, alpha_star, beta_star, max_sets = 200000L,
                                tol_bisect = 1e-7, bisect_iter = 80L,
                                qp_ridge = 1e-10) {
  sets <- .enumerate_experiment_sets(prep, max_sets = max_sets)
  p <- prep$p
  best <- list(regret = Inf, t = Inf, alpha = Inf, beta = Inf,
               alpha_ratio = Inf, beta_ratio = Inf, s = rep(0, p),
               x = rep(0L, p), gamma = rep(0, p))

  for (S in sets) {
    x <- .set_to_x(S, p)
    res <- .solve_minimax_for_fixed_x_enum(
      prep, x = x, alpha_star = alpha_star, beta_star = beta_star,
      tol_bisect = tol_bisect, bisect_iter = bisect_iter, qp_ridge = qp_ridge
    )
    if (is.finite(res$regret) && res$regret < best$regret) {
      best <- res
    }
  }
  if (!is.finite(best$regret)) {
    stop("The non-Gurobi minimax search failed to find a finite feasible design.",
         call. = FALSE)
  }
  best$result <- list(status = "OPTIMAL", solver = "quadprog", enumerated_sets = length(sets))
  best
}

.solve_alpha_oracle <- function(prep, gurobi_params = list(), status_label = "variance oracle") {
  .minimax_require("gurobi")
  p <- prep$p
  if (!prep$use_sets) {
    nvar <- 2L * p
    idx_x <- seq_len(p)
    idx_s <- (p + 1L):(2L * p)
    idx_y <- integer(0)
  } else {
    nvar <- 2L * p + prep$Msets
    idx_x <- seq_len(p)
    idx_s <- (p + 1L):(2L * p)
    idx_y <- (2L * p + 1L):(2L * p + prep$Msets)
  }
  prep$current_idx_y <- idx_y

  model <- list()
  model$modelname <- "alpha_star"
  model$modelsense <- "min"
  model$vtype <- c(rep("B", p), rep("C", p), if (prep$use_sets) rep("B", prep$Msets))
  model$lb <- rep(0, nvar)
  model$ub <- c(prep$x_max, prep$gamma_upper, if (prep$use_sets) rep(1, prep$Msets))
  model$Q <- .block_quadratic(prep$M_alpha, idx_s, nvar)
  model$obj <- rep(0, nvar)
  model$obj[idx_s] <- prep$q_alpha

  A <- .empty_sparse(0L, nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .add_x_constraints(A, rhs, sense, nvar, idx_x, prep)
  tmp <- .add_s_constraints(tmp$A, tmp$rhs, tmp$sense, nvar, idx_x, idx_s, prep)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  result <- gurobi::gurobi(model, params = gurobi_params)
  if (is.null(result$x)) {
    stop(sprintf("Gurobi did not return a solution for the %s. Status: %s", status_label, result$status),
         call. = FALSE)
  }
  if (!identical(result$status, "OPTIMAL")) {
    warning(sprintf("Gurobi status for the %s: %s", status_label, result$status), call. = FALSE)
  }
  sol <- result$x
  x <- as.integer(round(sol[idx_x]))
  s <- pmax(0, pmin(prep$gamma_upper, as.numeric(sol[idx_s])))
  gamma <- rep(0, p)
  gamma[x == 1L] <- s[x == 1L]
  list(
    x = x,
    s = s,
    gamma = gamma,
    alpha = prep$compute_alpha(s),
    beta = prep$compute_beta(s),
    result = result
  )
}

.solve_beta_oracle <- function(prep, gurobi_params = list()) {
  .minimax_require("gurobi")
  p <- prep$p
  if (!prep$use_sets) {
    nvar <- p
    idx_x <- seq_len(p)
    idx_y <- integer(0)
  } else {
    nvar <- p + prep$Msets
    idx_x <- seq_len(p)
    idx_y <- (p + 1L):(p + prep$Msets)
  }
  prep$current_idx_y <- idx_y

  model <- list()
  model$modelname <- "beta_star"
  model$modelsense <- "max"
  model$vtype <- c(rep("B", p), if (prep$use_sets) rep("B", prep$Msets))
  model$lb <- rep(0, nvar)
  model$ub <- c(prep$x_max, if (prep$use_sets) rep(1, prep$Msets))
  model$obj <- rep(0, nvar)
  model$obj[idx_x] <- prep$d * prep$gamma_upper

  A <- .empty_sparse(0L, nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .add_x_constraints(A, rhs, sense, nvar, idx_x, prep)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  result <- gurobi::gurobi(model, params = gurobi_params)
  if (is.null(result$x)) {
    stop(sprintf("Gurobi did not return a solution for the bias oracle. Status: %s", result$status),
         call. = FALSE)
  }
  if (!identical(result$status, "OPTIMAL")) {
    warning(sprintf("Gurobi status for the bias oracle: %s", result$status), call. = FALSE)
  }
  sol <- result$x
  x <- as.integer(round(sol[idx_x]))
  s <- prep$gamma_upper * x
  gamma <- rep(0, p)
  gamma[x == 1L] <- prep$gamma_upper[x == 1L]
  list(
    x = x,
    s = s,
    gamma = gamma,
    alpha = prep$compute_alpha(s),
    beta = prep$compute_beta(s),
    result = result
  )
}

.solve_minimax_miqcp <- function(prep, alpha_star, beta_star, gurobi_params = list()) {
  .minimax_require("gurobi")
  p <- prep$p
  if (!prep$use_sets) {
    nvar <- 2L * p + 1L
    idx_x <- seq_len(p)
    idx_s <- (p + 1L):(2L * p)
    idx_y <- integer(0)
    idx_t <- 2L * p + 1L
  } else {
    nvar <- 2L * p + prep$Msets + 1L
    idx_x <- seq_len(p)
    idx_s <- (p + 1L):(2L * p)
    idx_y <- (2L * p + 1L):(2L * p + prep$Msets)
    idx_t <- 2L * p + prep$Msets + 1L
  }
  prep$current_idx_y <- idx_y

  model <- list()
  model$modelname <- "minimax_design"
  model$modelsense <- "min"
  model$vtype <- c(rep("B", p), rep("C", p), if (prep$use_sets) rep("B", prep$Msets), "C")
  model$lb <- rep(0, nvar)
  model$ub <- c(prep$x_max, prep$gamma_upper, if (prep$use_sets) rep(1, prep$Msets), Inf)
  model$obj <- rep(0, nvar)
  model$obj[idx_t] <- 1

  A <- .empty_sparse(0L, nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .add_x_constraints(A, rhs, sense, nvar, idx_x, prep)
  tmp <- .add_s_constraints(tmp$A, tmp$rhs, tmp$sense, nvar, idx_x, idx_s, prep)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  q_alpha <- rep(0, nvar)
  q_alpha[idx_s] <- prep$q_alpha
  q_alpha[idx_t] <- -alpha_star
  rhs_alpha <- -prep$const_alpha
  scale_alpha <- .quad_scale(prep$M_alpha, q_alpha, rhs_alpha)

  q_beta <- rep(0, nvar)
  q_beta[idx_s] <- prep$q_beta
  if (beta_star > 0) {
    q_beta[idx_t] <- -beta_star
  }
  rhs_beta <- -prep$const_beta
  scale_beta <- .quad_scale(prep$M_beta, q_beta, rhs_beta)

  model$quadcon <- list(
    list(
      Qc = .block_quadratic(prep$M_alpha, idx_s, nvar, scale_alpha),
      q = scale_alpha * q_alpha,
      rhs = scale_alpha * rhs_alpha,
      sense = "<",
      name = "alpha_ratio"
    ),
    list(
      Qc = .block_quadratic(prep$M_beta, idx_s, nvar, scale_beta),
      q = scale_beta * q_beta,
      rhs = scale_beta * rhs_beta,
      sense = "<",
      name = "beta_ratio"
    )
  )

  result <- gurobi::gurobi(model, params = gurobi_params)
  if (is.null(result$x)) {
    stop(sprintf("Gurobi did not return a solution for the minimax program. Status: %s", result$status),
         call. = FALSE)
  }
  if (!identical(result$status, "OPTIMAL")) {
    warning(sprintf("Gurobi status for the minimax program: %s", result$status), call. = FALSE)
  }

  sol <- result$x
  x <- as.integer(round(sol[idx_x]))
  s <- pmax(0, pmin(prep$gamma_upper, as.numeric(sol[idx_s])))
  gamma <- rep(0, p)
  gamma[x == 1L] <- s[x == 1L]

  alpha <- prep$compute_alpha(s)
  beta <- prep$compute_beta(s)
  alpha_ratio <- .safe_ratio(alpha, alpha_star)
  beta_ratio <- .safe_ratio(beta, beta_star)
  regret <- max(alpha_ratio, beta_ratio)

  list(
    x = x,
    s = s,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    alpha_ratio = alpha_ratio,
    beta_ratio = beta_ratio,
    regret = regret,
    t = as.numeric(sol[idx_t]),
    result = result
  )
}

.format_selected <- function(x, names) {
  selected <- names[which(x == 1L)]
  if (length(selected) == 0L) {
    return("none")
  }
  paste(selected, collapse = ", ")
}
