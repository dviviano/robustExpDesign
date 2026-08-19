.minimax_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required.", package), call. = FALSE)
  }
}

.robustExpDesign_option <- function(name, default = NULL) {
  primary <- getOption(paste0("robustExpDesign", name), NULL)
  if (!is.null(primary)) return(primary)
  getOption(paste0("robustExpDesign.", name), default)
}


.gurobi_status <- function(result) {
  if (is.null(result$status) || length(result$status) == 0L) {
    "NO_STATUS"
  } else {
    as.character(result$status[[1L]])
  }
}

.require_gurobi_optimal <- function(result, label) {
  status <- .gurobi_status(result)
  if (is.null(result$x)) {
    stop(sprintf("Gurobi did not return a solution for the %s. Status: %s",
                 label, status), call. = FALSE)
  }
  if (!identical(status, "OPTIMAL")) {
    stop(sprintf(
      "Gurobi did not certify the %s as globally optimal. Status: %s",
      label, status
    ), call. = FALSE)
  }
  invisible(status)
}

.gurobi_result_scalar <- function(result, name) {
  value <- result[[name]]
  if (is.null(value) || length(value) == 0L) {
    return(NA_real_)
  }
  as.numeric(value[[1L]])
}

.record_gurobi_diagnostic <- function(model, result, objective_scale) {
  store <- .robustExpDesign_option("gurobi_diagnostics")
  if (!is.environment(store)) {
    return(invisible(NULL))
  }
  if (is.null(store$rows)) {
    store$rows <- list()
  }
  call_id <- length(store$rows) + 1L
  nvar <- if (is.null(model$obj)) {
    max(length(model$lb), length(model$ub), length(model$vtype))
  } else {
    length(model$obj)
  }
  store$rows[[call_id]] <- data.frame(
    call_id = call_id,
    model = if (is.null(model$modelname)) "" else as.character(model$modelname),
    status = if (is.null(result$status)) "NO_STATUS" else as.character(result$status),
    objective_scale = objective_scale,
    n_variables = nvar,
    n_linear_constraints = if (is.null(model$A)) 0L else nrow(model$A),
    n_quadratic_constraints = length(model$quadcon),
    has_integer_variables = !is.null(model$vtype) &&
      any(model$vtype %in% base::c("B", "I", "N")),
    runtime = .gurobi_result_scalar(result, "runtime"),
    iterations = .gurobi_result_scalar(result, "itercount"),
    barrier_iterations = .gurobi_result_scalar(result, "baritercount"),
    nodes = .gurobi_result_scalar(result, "nodecount"),
    objective_natural = .gurobi_result_scalar(result, "objval"),
    objective_bound_natural = .gurobi_result_scalar(result, "objbound"),
    mip_gap = .gurobi_result_scalar(result, "mipgap"),
    max_violation = .gurobi_result_scalar(result, "maxvio"),
    constraint_violation = .gurobi_result_scalar(result, "constrvio"),
    bound_violation = .gurobi_result_scalar(result, "boundvio"),
    integrality_violation = .gurobi_result_scalar(result, "intvio"),
    quadratic_constraint_violation = .gurobi_result_scalar(result, "qconstrvio"),
    dual_violation_scaled_objective = .gurobi_result_scalar(result, "dualvio"),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

.objective_coefficient_scale <- function(model) {
  values <- numeric(0)
  if (!is.null(model$obj)) values <- base::c(values, abs(as.numeric(model$obj)))
  if (!is.null(model$Q)) values <- base::c(values, abs(as.numeric(model$Q)))
  if (!is.null(model$objcon)) values <- base::c(values, abs(as.numeric(model$objcon)))
  values <- values[is.finite(values) & values > 0]
  if (length(values) == 0L) return(1)
  scale <- 1 / max(values)
  min(1e300, max(1e-300, scale))
}

.gurobi_solve <- function(model, params = list()) {
  objective_scale <- .robustExpDesign_option("gurobi_objective_scale", "auto")
  if (is.character(objective_scale) && length(objective_scale) == 1L &&
      identical(tolower(objective_scale), "auto")) {
    objective_scale <- .objective_coefficient_scale(model)
  } else if (!is.numeric(objective_scale) || length(objective_scale) != 1L ||
             !is.finite(objective_scale) || objective_scale <= 0) {
    stop("The Gurobi objective scale must be 'auto' or a positive finite scalar.",
         call. = FALSE)
  }
  scaled_models <- .robustExpDesign_option("gurobi_objective_scale_models")
  if (!is.null(scaled_models)) {
    model_name <- if (is.null(model$modelname)) "" else model$modelname
    if (!(model_name %in% scaled_models)) objective_scale <- 1
  }

  scaled_model <- model
  if (!is.null(scaled_model$obj)) {
    scaled_model$obj <- objective_scale * scaled_model$obj
  }
  if (!is.null(scaled_model$Q)) {
    scaled_model$Q <- objective_scale * scaled_model$Q
  }
  if (!is.null(scaled_model$objcon)) {
    scaled_model$objcon <- objective_scale * scaled_model$objcon
  }

  result <- gurobi::gurobi(scaled_model, params = params)

  # Return objective and dual quantities in the model's natural units.
  result$objective_scale <- objective_scale
  for (name in intersect(base::c("objval", "objbound", "objboundc",
                           "poolobjval", "scenobjval"),
                         names(result))) {
    result[[paste0(name, "_scaled")]] <- result[[name]]
    result[[name]] <- result[[name]] / objective_scale
  }
  for (name in intersect(base::c("pi", "rc", "qcpi"), names(result))) {
    result[[name]] <- result[[name]] / objective_scale
  }
  .record_gurobi_diagnostic(model, result, objective_scale)
  result
}

.normalize_sense <- function(sense) {
  if (is.null(sense)) {
    return(character(0))
  }
  out <- as.character(sense)
  out[out %in% base::c("<=", "=<")] <- "<"
  out[out %in% base::c(">=", "=>")] <- ">"
  if (!all(out %in% base::c("<", ">", "="))) {
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
  if (!is.numeric(x) || is.logical(x) || !is.null(dim(x))) {
    stop(sprintf("%s must be a numeric vector.", name), call. = FALSE)
  }
  x_names <- names(x)
  x <- as.numeric(x)
  if (!is.null(x_names)) names(x) <- x_names
  if (!is.null(p) && length(x) != p) {
    stop(sprintf("%s must have length %d.", name, p), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("%s must contain only finite values.", name), call. = FALSE)
  }
  x
}

.validate_tolerance <- function(tol, name = "tol") {
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop(sprintf("%s must be a nonnegative finite scalar.", name), call. = FALSE)
  }
  as.numeric(tol)
}

.relative_zero <- function(value, scale, tol = 1e-9) {
  tol <- .validate_tolerance(tol)
  value <- as.numeric(value)
  scale <- as.numeric(scale)
  if (length(value) != 1L || length(scale) != 1L ||
      !is.finite(value) || !is.finite(scale) || scale < 0) {
    stop("value and scale must be finite scalars with scale >= 0.", call. = FALSE)
  }
  if (scale == 0) return(value == 0)
  abs(value) <= tol * scale
}

.safe_nonnegative <- function(value, scale = NULL, tol = 1e-10,
                              name = "value") {
  tol <- .validate_tolerance(tol)
  value <- as.numeric(value)
  if (length(value) != 1L || !is.finite(value)) {
    stop(sprintf("%s must be a finite numeric scalar.", name), call. = FALSE)
  }
  if (is.null(scale)) scale <- abs(value)
  scale <- as.numeric(scale)
  if (length(scale) != 1L || !is.finite(scale) || scale < 0) {
    stop("scale must be a nonnegative finite scalar.", call. = FALSE)
  }
  threshold <- tol * scale
  if (value < -threshold) {
    stop(sprintf("%s is negative beyond numerical tolerance.", name),
         call. = FALSE)
  }
  if (value < 0) 0 else value
}

.clamp_relative_zero <- function(value, scale, tol = 1e-9,
                                 name = "value") {
  value <- .safe_nonnegative(value, scale = scale, tol = tol, name = name)
  if (.relative_zero(value, scale, tol = tol)) 0 else value
}

.validate_psd_matrix <- function(Sigma, name = "Sigma", tol = 1e-10,
                                 action = base::c("error", "project")) {
  tol <- .validate_tolerance(tol, "psd_tol")
  action <- match.arg(action)
  Sigma <- as.matrix(Sigma)
  if (!is.numeric(Sigma) || nrow(Sigma) != ncol(Sigma) ||
      any(!is.finite(Sigma))) {
    stop(sprintf("%s must be a finite square numeric matrix.", name),
         call. = FALSE)
  }
  Sigma <- 0.5 * (Sigma + t(Sigma))
  eig <- eigen(Sigma, symmetric = TRUE)
  values <- as.numeric(eig$values)
  spectral_scale <- if (length(values) == 0L) 0 else max(abs(values))
  min_eigenvalue <- if (length(values) == 0L) 0 else min(values)
  hard_threshold <- tol * spectral_scale
  roundoff_threshold <- 64 * .Machine$double.eps *
    max(1, nrow(Sigma)) * spectral_scale

  if (min_eigenvalue < -hard_threshold) {
    stop(sprintf(
      "%s must be positive semidefinite: minimum eigenvalue %.6g is below the relative tolerance %.6g.",
      name, min_eigenvalue, -hard_threshold
    ), call. = FALSE)
  }

  projected <- FALSE
  projection_norm <- 0
  if (min_eigenvalue < 0) {
    if (min_eigenvalue < -roundoff_threshold && identical(action, "error")) {
      stop(sprintf(
        "%s has a negative eigenvalue %.6g within psd_tol but above roundoff. Use psd_action = 'project' to opt in to projection.",
        name, min_eigenvalue
      ), call. = FALSE)
    }
    clipped <- pmax(values, 0)
    projected_matrix <- eig$vectors %*% diag(clipped, nrow = length(clipped)) %*%
      t(eig$vectors)
    projected_matrix <- 0.5 * (projected_matrix + t(projected_matrix))
    projection_norm <- sqrt(sum((projected_matrix - Sigma)^2))
    Sigma <- projected_matrix
    projected <- TRUE
  }

  list(
    matrix = Sigma,
    diagnostics = list(
      minimum_eigenvalue = min_eigenvalue,
      spectral_scale = spectral_scale,
      relative_minimum_eigenvalue = if (spectral_scale > 0) {
        min_eigenvalue / spectral_scale
      } else {
        0
      },
      psd_tol = tol,
      psd_action = action,
      projected = projected,
      projection_frobenius_norm = projection_norm
    )
  )
}

.psd_square_root <- function(Sigma, tol = 1e-10,
                             action = base::c("error", "project"),
                             name = "Sigma") {
  checked <- .validate_psd_matrix(Sigma, name = name, tol = tol,
                                  action = action)
  eig <- eigen(checked$matrix, symmetric = TRUE)
  root <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)),
                               nrow = length(eig$values))
  list(root = root, matrix = checked$matrix,
       diagnostics = checked$diagnostics)
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
  if (!is.numeric(A)) {
    stop("selection_constraints$A must be numeric.", call. = FALSE)
  }
  if (ncol(A) != p) {
    stop("selection_constraints$A must have one column per experiment.", call. = FALSE)
  }
  sense <- .normalize_sense(sense)
  if (!is.numeric(rhs) || is.logical(rhs) || !is.null(dim(rhs))) {
    stop("selection_constraints$rhs must be a numeric vector.", call. = FALSE)
  }
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
    if (!is.numeric(S) || is.logical(S) || !is.null(dim(S)) ||
        any(!is.finite(S)) || any(abs(S - round(S)) > 1e-9)) {
      stop("feasible_sets must contain finite integer index vectors.",
           call. = FALSE)
    }
    S <- unique(as.integer(round(S)))
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
    rhs = base::c(rhs_vec, rhs),
    sense = base::c(sense_vec, .normalize_sense(sense))
  )
}

.add_x_constraints <- function(A_mat, rhs_vec, sense_vec, nvar, idx_x, prep) {
  p <- prep$p

  row <- .empty_sparse(1L, nvar)
  row[1L, idx_x] <- 1
  tmp <- list(A = rbind(A_mat, row), rhs = base::c(rhs_vec, prep$h),
              sense = base::c(sense_vec, "<"))
  A_mat <- tmp$A
  rhs_vec <- tmp$rhs
  sense_vec <- tmp$sense

  if (prep$min_experiments > 0) {
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_x] <- 1
    tmp <- list(A = rbind(A_mat, row), rhs = base::c(rhs_vec, prep$min_experiments),
                sense = base::c(sense_vec, ">"))
    A_mat <- tmp$A
    rhs_vec <- tmp$rhs
    sense_vec <- tmp$sense
  }

  if (!is.null(prep$selection_constraints)) {
    for (r in seq_len(nrow(prep$selection_constraints$A))) {
      row <- .empty_sparse(1L, nvar)
      row[1L, idx_x] <- prep$selection_constraints$A[r, ]
      A_mat <- rbind(A_mat, row)
      rhs_vec <- base::c(rhs_vec, prep$selection_constraints$rhs[r])
      sense_vec <- base::c(sense_vec, prep$selection_constraints$sense[r])
    }
  }

  if (prep$use_sets) {
    idx_y <- prep$current_idx_y
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_y] <- 1
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, 1)
    sense_vec <- base::c(sense_vec, "=")

    for (j in seq_len(p)) {
      row <- .empty_sparse(1L, nvar)
      row[1L, idx_x[j]] <- 1
      row[1L, idx_y] <- -prep$Aset[, j]
      A_mat <- rbind(A_mat, row)
      rhs_vec <- base::c(rhs_vec, 0)
      sense_vec <- base::c(sense_vec, "=")
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
    rhs_vec <- base::c(rhs_vec, 0)
    sense_vec <- base::c(sense_vec, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_s[j]] <- -1
    row[1L, idx_x[j]] <- prep$gamma_lower[j]
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, 0)
    sense_vec <- base::c(sense_vec, "<")
  }
  list(A = A_mat, rhs = rhs_vec, sense = sense_vec)
}

.block_quadratic <- function(M, idx, nvar, scale = 1) {
  nz <- which(M != 0, arr.ind = TRUE)
  if (nrow(nz) == 0L) {
    return(Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = base::c(nvar, nvar)))
  }
  Matrix::sparseMatrix(
    i = idx[nz[, 1]],
    j = idx[nz[, 2]],
    x = as.numeric(scale * M[nz]),
    dims = base::c(nvar, nvar)
  )
}

.block_quadratic_sum <- function(blocks, nvar) {
  rows <- integer(0)
  cols <- integer(0)
  vals <- numeric(0)
  for (block in blocks) {
    M <- block$M
    idx <- block$idx
    scale <- if (is.null(block$scale)) 1 else block$scale
    nz <- which(M != 0, arr.ind = TRUE)
    if (nrow(nz) == 0L) {
      next
    }
    rows <- base::c(rows, idx[nz[, 1]])
    cols <- base::c(cols, idx[nz[, 2]])
    vals <- base::c(vals, as.numeric(scale * M[nz]))
  }
  Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = base::c(nvar, nvar))
}

.constraint_scale <- function(values) {
  values <- abs(as.numeric(values))
  values <- values[is.finite(values) & values > 0]
  if (length(values) == 0L) return(1)
  scale <- 1 / max(values)
  min(1e300, max(1e-300, scale))
}

.quad_scale <- function(M, q, rhs) {
  .constraint_scale(base::c(as.numeric(M), as.numeric(q), rhs))
}

.quad_scale_multi <- function(M_list, q, rhs) {
  vals <- base::c(as.numeric(q), rhs)
  for (M in M_list) vals <- base::c(vals, as.numeric(M))
  .constraint_scale(vals)
}

.safe_ratio <- function(value, reference, tol = 1e-10) {
  tol <- .validate_tolerance(tol)
  value <- as.numeric(value)
  reference <- as.numeric(reference)
  if (length(value) != 1L || length(reference) != 1L ||
      !is.finite(value) || !is.finite(reference)) {
    stop("value and reference must be finite numeric scalars.", call. = FALSE)
  }
  scale <- max(abs(value), abs(reference))
  value <- .safe_nonnegative(value, scale = scale, tol = tol, name = "value")
  reference <- .safe_nonnegative(reference, scale = scale, tol = tol,
                                 name = "reference")
  if (reference == 0) return(if (value == 0) 1 else Inf)
  as.numeric(value / reference)
}

.validate_regret_epigraph <- function(regret, solver_t, tol = 1e-6) {
  tol <- .validate_tolerance(tol)
  regret <- as.numeric(regret)
  solver_t <- as.numeric(solver_t)
  if (length(regret) != 1L || length(solver_t) != 1L ||
      !is.finite(regret) || !is.finite(solver_t)) {
    stop("Validated regret and solver epigraph value must be finite scalars.",
         call. = FALSE)
  }
  gap <- regret - solver_t
  gap_tol <- tol * max(1, abs(regret), abs(solver_t))
  if (gap > gap_tol) {
    stop(sprintf(
      "The recomputed regret exceeds the solver epigraph by %.6g, above tolerance %.6g.",
      gap, gap_tol
    ), call. = FALSE)
  }
  list(t = max(regret, solver_t), solver_t = solver_t,
       validated_gap = gap, tolerance = gap_tol)
}

.implied_gamma <- function(a_exp, omega, tol = NULL) {
  out <- rep(NA_real_, length(omega))
  ok <- as.numeric(omega) != 0
  out[ok] <- as.numeric(a_exp)[ok] / as.numeric(omega)[ok]
  out
}

.beta_zero_feasible <- function(prep, x) {
  x <- as.numeric(x)
  if (length(x) != prep$p || any(!is.finite(x))) {
    stop("x must be a finite vector with one entry per experiment.",
         call. = FALSE)
  }
  selected <- x > 0.5
  if (identical(prep$weight_mode, "a")) {
    lower <- prep$a_exp_lower * as.numeric(selected)
    upper <- prep$a_exp_upper * as.numeric(selected)
    relevant <- prep$bias_weights > 0
    return(all(!relevant |
                 (prep$omega >= lower & prep$omega <= upper)))
  }

  relevant <- prep$d > 0
  all(!relevant |
        (selected & prep$gamma_lower <= 1 & prep$gamma_upper >= 1))
}

.beta_star_is_zero <- function(prep, beta_star) {
  if (!is.null(prep$beta_zero)) return(isTRUE(prep$beta_zero))
  .relative_zero(beta_star, prep$beta_scale, tol = prep$zero_tol)
}

.add_a_constraints <- function(A_mat, rhs_vec, sense_vec, nvar, idx_x, idx_a,
                               prep) {
  p <- prep$p
  for (j in seq_len(p)) {
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- 1
    row[1L, idx_x[j]] <- -prep$a_exp_upper[j]
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, 0)
    sense_vec <- base::c(sense_vec, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- -1
    row[1L, idx_x[j]] <- prep$a_exp_lower[j]
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, 0)
    sense_vec <- base::c(sense_vec, "<")
  }
  list(A = A_mat, rhs = rhs_vec, sense = sense_vec)
}

.add_abs_a_constraints <- function(A_mat, rhs_vec, sense_vec, nvar, idx_a,
                                   idx_r, idx_z, prep) {
  p <- prep$p
  for (j in seq_len(p)) {
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- 1
    row[1L, idx_r[j]] <- -1
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, 0)
    sense_vec <- base::c(sense_vec, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- -1
    row[1L, idx_r[j]] <- -1
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, 0)
    sense_vec <- base::c(sense_vec, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- 1
    row[1L, idx_z[j]] <- -1
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, prep$omega[j])
    sense_vec <- base::c(sense_vec, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- -1
    row[1L, idx_z[j]] <- -1
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, -prep$omega[j])
    sense_vec <- base::c(sense_vec, "<")
  }
  list(A = A_mat, rhs = rhs_vec, sense = sense_vec)
}

.prepare_minimax_problem <- function(Sigma_obs, v2, n_total, omega,
                                     costs = NULL, bias_weights = NULL,
                                     h = NULL, x_max = NULL,
                                     feasible_sets = NULL,
                                     selection_constraints = NULL,
                                     min_experiments = 0L,
                                     n_min = 0,
                                     gamma_lower = NULL,
                                     gamma_upper = NULL,
                                     a_exp_lower = NULL,
                                     a_exp_upper = NULL,
                                     force_gamma_unit_interval = TRUE,
                                     costs_alias = NULL,
                                     bias_weights_alias = NULL,
                                     zero_tol = 1e-9,
                                     psd_tol = 1e-10,
                                     psd_action = base::c("error", "project")) {
  .minimax_require("Matrix")

  if (!is.null(costs_alias)) {
    if (!is.null(costs) && !isTRUE(all.equal(as.numeric(costs), as.numeric(costs_alias)))) {
      stop("Use either costs or c, or supply identical values for both.", call. = FALSE)
    }
    costs <- costs_alias
  }
  if (!is.null(bias_weights_alias)) {
    if (!is.null(bias_weights) && !isTRUE(all.equal(as.numeric(bias_weights), as.numeric(bias_weights_alias)))) {
      stop("Use either bias_weights or k, or supply identical values for both.", call. = FALSE)
    }
    bias_weights <- bias_weights_alias
  }

  omega_original <- .numeric_vector(omega, "omega")
  p <- length(omega_original)
  if (p < 1L) {
    stop("omega must have positive length.", call. = FALSE)
  }
  arm_names <- names(omega_original)
  target_scale <- max(abs(omega_original))
  if (target_scale == 0) target_scale <- 1
  # Work in normalized target units. This prevents overflow/underflow and makes
  # target-unit invariance exact up to the optimizer's dimensionless tolerances.
  omega <- as.numeric(omega_original) / target_scale
  names(omega) <- arm_names

  Sigma_obs <- as.matrix(Sigma_obs)
  if (!is.numeric(Sigma_obs) || any(dim(Sigma_obs) != base::c(p, p))) {
    stop("Sigma_obs must be a p by p numeric matrix with p = length(omega).", call. = FALSE)
  }
  if (any(!is.finite(Sigma_obs))) {
    stop("Sigma_obs must contain only finite values.", call. = FALSE)
  }
  psd_action <- match.arg(psd_action)
  zero_tol <- .validate_tolerance(zero_tol, "zero_tol")
  checked_sigma <- .validate_psd_matrix(
    Sigma_obs, name = "Sigma_obs", tol = psd_tol, action = psd_action
  )
  Sigma_obs <- checked_sigma$matrix
  psd_diagnostics <- checked_sigma$diagnostics

  v2 <- .numeric_vector(v2, "v2", p)
  if (any(v2 < 0)) {
    stop("v2 must be nonnegative.", call. = FALSE)
  }
  if (!is.numeric(n_total) || length(n_total) != 1L || !is.finite(n_total) || n_total <= 0) {
    stop("n_total must be a positive scalar budget.", call. = FALSE)
  }
  if (!is.numeric(n_min) || length(n_min) != 1L || !is.finite(n_min) || n_min < 0) {
    stop("n_min must be a nonnegative finite scalar.", call. = FALSE)
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
  if (any(!(x_max %in% base::c(0, 1)))) {
    stop("x_max must contain only 0 or 1.", call. = FALSE)
  }

  gamma_lower <- .numeric_vector(gamma_lower, "gamma_lower", p, default = rep(0, p))
  gamma_upper <- .numeric_vector(gamma_upper, "gamma_upper", p, default = rep(1, p))
  if (!is.logical(force_gamma_unit_interval) || length(force_gamma_unit_interval) != 1L ||
      is.na(force_gamma_unit_interval)) {
    stop("force_gamma_unit_interval must be TRUE or FALSE.", call. = FALSE)
  }
  weight_mode <- if (isTRUE(force_gamma_unit_interval)) "gamma" else "a"
  if (isTRUE(force_gamma_unit_interval)) {
    if (any(gamma_lower < 0) || any(gamma_upper > 1) || any(gamma_lower > gamma_upper)) {
      stop("gamma_lower and gamma_upper must satisfy 0 <= lower <= upper <= 1.", call. = FALSE)
    }
  } else {
    if (any(gamma_lower > gamma_upper)) {
      stop("When force_gamma_unit_interval = FALSE, gamma bounds must satisfy lower <= upper and be finite.", call. = FALSE)
    }
  }

  if (identical(weight_mode, "a")) {
    if (is.null(a_exp_lower) && is.null(a_exp_upper)) {
      a_from_lower <- omega * gamma_lower
      a_from_upper <- omega * gamma_upper
      a_exp_lower <- pmin(a_from_lower, a_from_upper)
      a_exp_upper <- pmax(a_from_lower, a_from_upper)
    } else if (is.null(a_exp_lower) || is.null(a_exp_upper)) {
      stop("Supply both a_exp_lower and a_exp_upper, or neither.", call. = FALSE)
    } else {
      # Direct linear-weight bounds are supplied in the user's target units.
      a_exp_lower <- .numeric_vector(a_exp_lower, "a_exp_lower", p) / target_scale
      a_exp_upper <- .numeric_vector(a_exp_upper, "a_exp_upper", p) / target_scale
    }
    a_exp_lower <- .numeric_vector(a_exp_lower, "a_exp_lower", p)
    a_exp_upper <- .numeric_vector(a_exp_upper, "a_exp_upper", p)
    if (any(a_exp_lower > a_exp_upper)) {
      stop("a_exp_lower and a_exp_upper must satisfy lower <= upper.", call. = FALSE)
    }
  } else {
    a_from_lower <- omega * gamma_lower
    a_from_upper <- omega * gamma_upper
    a_exp_lower <- pmin(a_from_lower, a_from_upper)
    a_exp_upper <- pmax(a_from_lower, a_from_upper)
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
  if (n_min > 0) {
    eligible_costs <- sort(costs[x_max > 0.5])
    if (length(eligible_costs) == 0L || h < 1L) {
      stop("n_min > 0 requires at least one eligible experiment and h >= 1.",
           call. = FALSE)
    }
    min_required_count <- max(1L, min_experiments)
    if (length(eligible_costs) < min_required_count) {
      stop("Not enough eligible experiments to satisfy min_experiments with n_min.",
           call. = FALSE)
    }
    cheapest_required <- sum(eligible_costs[seq_len(min_required_count)] * n_min)
    if (cheapest_required > n_total +
        1e-10 * max(n_total, cheapest_required)) {
      stop("n_min is infeasible: even the cheapest required selected experiments exceed n_total.",
           call. = FALSE)
    }
  }

  selection_constraints <- .validate_selection_constraints(selection_constraints, p)
  set_info <- .process_feasible_sets(feasible_sets, p, h, x_max, min_experiments,
                                     selection_constraints)
  if (n_min > 0 && isTRUE(set_info$use_sets)) {
    budget_tol <- 1e-10 * n_total
    keep <- vapply(set_info$sets, function(S) {
      sum(costs[S] * n_min) <= n_total + budget_tol
    }, logical(1))
    set_info$sets <- set_info$sets[keep]
    if (length(set_info$sets) == 0L) {
      stop("No feasible experiment set can satisfy n_min within n_total.",
           call. = FALSE)
    }
    set_info$Msets <- length(set_info$sets)
    set_info$Aset <- matrix(0, nrow = set_info$Msets, ncol = p)
    for (m in seq_along(set_info$sets)) {
      S <- set_info$sets[[m]]
      if (length(S) > 0L) set_info$Aset[m, S] <- 1
    }
  }

  v <- sqrt(v2)
  a <- abs(omega) * v * sqrt(costs)
  a_direct <- v * sqrt(costs)
  d <- bias_weights * abs(omega)
  Domega <- diag(omega, nrow = p, ncol = p)
  Qobs <- Domega %*% Sigma_obs %*% Domega
  Qobs <- 0.5 * (Qobs + t(Qobs))
  Sigma_obs <- 0.5 * (Sigma_obs + t(Sigma_obs))

  M_alpha <- (a %o% a) / n_total + Qobs
  M_alpha <- 0.5 * (M_alpha + t(M_alpha))
  q_alpha <- -2 * as.numeric(Qobs %*% rep(1, p))
  const_alpha <- as.numeric(t(rep(1, p)) %*% Qobs %*% rep(1, p))

  M_beta <- d %o% d
  q_beta <- -2 * sum(d) * d
  const_beta <- sum(d)^2

  M_alpha_a <- Sigma_obs
  M_alpha_a <- 0.5 * (M_alpha_a + t(M_alpha_a))
  M_exp_a <- (a_direct %o% a_direct) / n_total
  M_exp_a <- 0.5 * (M_exp_a + t(M_exp_a))
  q_alpha_a <- -2 * as.numeric(Sigma_obs %*% omega)
  const_alpha_a <- as.numeric(t(omega) %*% Sigma_obs %*% omega)
  M_beta_a <- bias_weights %o% bias_weights

  beta_linear_scale <- if (identical(weight_mode, "a")) {
    max_dev <- pmax(abs(omega - a_exp_lower), abs(omega - a_exp_upper),
                    abs(omega))
    sum(bias_weights * max_dev)
  } else {
    sum(d)
  }
  beta_scale <- beta_linear_scale^2
  alpha_scale <- max(base::c(abs(as.numeric(M_alpha)), abs(q_alpha),
                       abs(const_alpha), abs(as.numeric(M_alpha_a)),
                       abs(as.numeric(M_exp_a)), abs(q_alpha_a),
                       abs(const_alpha_a), 0))

  compute_alpha <- function(s) {
    s <- as.numeric(s)
    a_exp <- if (identical(weight_mode, "a")) s else omega * s
    a_obs <- omega - a_exp
    value <- as.numeric(t(a_obs) %*% Sigma_obs %*% a_obs +
                          (sum(abs(a_exp) * a_direct)^2) / n_total)
    .safe_nonnegative(value, scale = alpha_scale, tol = zero_tol,
                      name = "alpha")
  }
  compute_beta <- function(s) {
    s <- as.numeric(s)
    a_exp <- if (identical(weight_mode, "a")) s else omega * s
    a_obs <- omega - a_exp
    as.numeric((sum(bias_weights * abs(a_obs)))^2)
  }
  compute_allocation <- function(s, x = NULL) {
    s <- as.numeric(s)
    if (identical(weight_mode, "a")) {
      score <- abs(s) * v
    } else {
      score <- abs(omega) * v * abs(s)
    }
    denom <- sum(score * sqrt(costs))
    n_opt <- rep(0, p)
    if (denom > 0) {
      n_opt <- n_total * (score / sqrt(costs)) / denom
      return(n_opt)
    }

    # When the experimental contribution is identically zero, allocation does
    # not affect risk. Still return a budget-feasible allocation whenever an
    # experiment is selected, as required by the equality budget constraint.
    selected <- if (is.null(x)) {
      abs(s) > 0
    } else {
      as.numeric(x) > 0.5
    }
    if (any(selected)) {
      j <- which(selected)[which.min(costs[selected])]
      n_opt[j] <- n_total / costs[j]
    }
    n_opt
  }

  if (is.null(arm_names) || any(arm_names == "")) {
    arm_names <- paste0("arm_", seq_len(p))
  }

  base::c(
    list(
      p = p,
      Sigma_obs = Sigma_obs,
      v2 = v2,
      v = v,
      n_total = n_total,
      omega = omega,
      omega_original = omega_original,
      target_scale = target_scale,
      risk_rescale = target_scale^2,
      costs = costs,
      bias_weights = bias_weights,
      h = h,
      x_max = x_max,
      min_experiments = min_experiments,
      n_min = n_min,
      zero_tol = zero_tol,
      psd_tol = psd_tol,
      psd_action = psd_action,
      psd_diagnostics = psd_diagnostics,
      alpha_scale = alpha_scale,
      beta_linear_scale = beta_linear_scale,
      beta_scale = beta_scale,
      weight_mode = weight_mode,
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper,
      a_exp_lower = a_exp_lower,
      a_exp_upper = a_exp_upper,
      force_gamma_unit_interval = force_gamma_unit_interval,
      selection_constraints = selection_constraints,
      arm_names = arm_names,
      a = a,
      a_direct = a_direct,
      d = d,
      Qobs = Qobs,
      M_alpha = M_alpha,
      q_alpha = q_alpha,
      const_alpha = const_alpha,
      M_beta = M_beta,
      q_beta = q_beta,
      const_beta = const_beta,
      M_alpha_a = M_alpha_a,
      M_exp_a = M_exp_a,
      q_alpha_a = q_alpha_a,
      const_alpha_a = const_alpha_a,
      M_beta_a = M_beta_a,
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
  solver <- match.arg(solver, base::c("auto", "gurobi", "quadprog"))
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
        sets <- base::c(sets, list(integer(0)))
      } else {
        sets <- base::c(sets, utils::combn(idx, k, simplify = FALSE))
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

.qp_min_alpha <- function(prep, lb, ub, bias_min = NULL, bias_max = NULL,
                          bias_equal = NULL, qp_ridge = 1e-10, tol = 1e-8) {
  .minimax_require("quadprog")
  p <- prep$p
  lb <- as.numeric(lb)
  ub <- as.numeric(ub)
  if (length(lb) != p || length(ub) != p) {
    stop("Internal error: lb and ub have wrong length.", call. = FALSE)
  }
  if (any(lb > ub + tol)) return(NULL)
  lb <- pmin(lb, ub)

  # Normalize the objective before adding the ridge. This makes qp_ridge a
  # relative, dimensionless regularizer and preserves target-unit invariance.
  objective_scale <- max(base::c(abs(as.numeric(prep$M_alpha)),
                           abs(prep$q_alpha), 0))
  if (objective_scale > 0) {
    M_qp <- prep$M_alpha / objective_scale
    q_qp <- prep$q_alpha / objective_scale
  } else {
    M_qp <- prep$M_alpha
    q_qp <- prep$q_alpha
  }
  qp_ridge <- .validate_tolerance(qp_ridge, "qp_ridge")
  ridge <- if (qp_ridge > 0) qp_ridge else 64 * .Machine$double.eps
  Dmat <- 2 * M_qp + diag(ridge, p)
  dvec <- -q_qp

  Amat <- cbind(diag(p), -diag(p))
  bvec <- base::c(lb, -ub)
  meq <- 0L

  bias_scale <- prep$beta_linear_scale
  bias_vec <- if (bias_scale > 0) prep$d / bias_scale else prep$d
  normalize_bias_rhs <- function(x) {
    if (is.null(x) || !is.finite(x) || bias_scale == 0) x else x / bias_scale
  }
  bias_equal_n <- normalize_bias_rhs(bias_equal)
  bias_min_n <- normalize_bias_rhs(bias_min)
  bias_max_n <- normalize_bias_rhs(bias_max)
  constraint_tol <- max(tol, 256 * .Machine$double.eps)

  if (!is.null(bias_equal_n) && is.finite(bias_equal_n) &&
      max(abs(bias_vec)) > 0) {
    lower_sum <- sum(bias_vec * lb)
    upper_sum <- sum(bias_vec * ub)
    if (bias_equal_n < lower_sum - 10 * constraint_tol ||
        bias_equal_n > upper_sum + 10 * constraint_tol) return(NULL)
    Amat <- cbind(bias_vec, Amat)
    bvec <- base::c(bias_equal_n, bvec)
    meq <- 1L
  }

  if (!is.null(bias_min_n) && is.finite(bias_min_n)) {
    if (bias_min_n > sum(bias_vec * ub) + 10 * constraint_tol) return(NULL)
    if (bias_min_n > sum(bias_vec * lb) + constraint_tol) {
      Amat <- cbind(Amat, bias_vec)
      bvec <- base::c(bvec, bias_min_n)
    }
  }
  if (!is.null(bias_max_n) && is.finite(bias_max_n)) {
    if (bias_max_n < sum(bias_vec * lb) - 10 * constraint_tol) return(NULL)
    if (bias_max_n < sum(bias_vec * ub) - constraint_tol) {
      Amat <- cbind(Amat, -bias_vec)
      bvec <- base::c(bvec, -bias_max_n)
    }
  }

  out <- tryCatch(
    quadprog::solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat,
                       bvec = bvec, meq = meq),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  sol <- as.numeric(out$solution)
  sol <- pmax(lb, pmin(ub, sol))
  bias_value <- sum(bias_vec * sol)
  if (!is.null(bias_min_n) && is.finite(bias_min_n) &&
      bias_value < bias_min_n - 10 * constraint_tol) return(NULL)
  if (!is.null(bias_max_n) && is.finite(bias_max_n) &&
      bias_value > bias_max_n + 10 * constraint_tol) return(NULL)
  if (!is.null(bias_equal_n) && is.finite(bias_equal_n) &&
      max(abs(bias_vec)) > 0 &&
      abs(bias_value - bias_equal_n) > 10 * constraint_tol) return(NULL)
  sol
}

.beta_oracle_s_for_x <- function(prep, x, tol = 1e-12) {
  p <- prep$p
  x <- as.numeric(x)
  lb <- prep$gamma_lower * x
  ub <- prep$gamma_upper * x
  d <- prep$d
  d_scale <- max(d)
  if (d_scale == 0) return(lb)
  d <- d / d_scale
  target <- sum(d)
  lower_sum <- sum(d * lb)
  upper_sum <- sum(d * ub)
  target_proj <- min(max(target, lower_sum), upper_sum)

  sol <- lb
  remaining <- target_proj - lower_sum
  tol <- .validate_tolerance(tol)
  threshold <- tol * max(1, abs(target_proj), abs(lower_sum), abs(upper_sum))
  if (remaining <= threshold) return(sol)

  for (j in order(d, decreasing = TRUE)) {
    if (d[j] == 0) next
    add <- min(ub[j] - lb[j], remaining / d[j])
    if (add > 0) {
      sol[j] <- sol[j] + add
      remaining <- remaining - d[j] * add
    }
    if (remaining <= threshold) break
  }
  pmax(lb, pmin(ub, sol))
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
    s <- .beta_oracle_s_for_x(prep, x)
    beta <- prep$compute_beta(s)
    alpha <- prep$compute_alpha(s)
    tie_scale <- max(abs(beta), abs(best$beta))
    tied <- is.finite(tie_scale) &&
      abs(beta - best$beta) <= prep$zero_tol * tie_scale
    if (beta < best$beta || (tied && alpha < best$alpha)) {
      gamma <- rep(0, p)
      gamma[x == 1L] <- s[x == 1L]
      best <- list(beta = beta, alpha = alpha, s = s,
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

  if (.beta_star_is_zero(prep, beta_star)) {
    s <- .qp_min_alpha(prep, lb = lb, ub = ub, bias_equal = total_d,
                       qp_ridge = qp_ridge)
    if (is.null(s)) {
      return(list(regret = Inf, t = Inf, alpha = Inf, beta = Inf,
                  s = rep(0, p), x = as.integer(x), gamma = rep(0, p)))
    }
    alpha <- prep$compute_alpha(s)
    beta <- .clamp_relative_zero(
      prep$compute_beta(s), prep$beta_scale, prep$zero_tol, "beta"
    )
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
    bias_max <- total_d + bias_bound
    if (bias_min <= 0) {
      bias_min <- NULL
    }
    s_mid <- .qp_min_alpha(prep, lb = lb, ub = ub, bias_min = bias_min,
                           bias_max = bias_max,
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

.a_gurobi_indices <- function(prep, include_t = FALSE) {
  p <- prep$p
  offset <- 0L
  idx_x <- seq_len(p)
  offset <- offset + p
  idx_a <- offset + seq_len(p)
  offset <- offset + p
  idx_r <- offset + seq_len(p)
  offset <- offset + p
  idx_z <- offset + seq_len(p)
  offset <- offset + p
  if (prep$use_sets) {
    idx_y <- offset + seq_len(prep$Msets)
    offset <- offset + prep$Msets
  } else {
    idx_y <- integer(0)
  }
  if (isTRUE(include_t)) {
    idx_t <- offset + 1L
    offset <- offset + 1L
  } else {
    idx_t <- integer(0)
  }
  list(nvar = offset, idx_x = idx_x, idx_a = idx_a, idx_r = idx_r,
       idx_z = idx_z, idx_y = idx_y, idx_t = idx_t)
}

.a_gurobi_bounds <- function(prep, idx, include_t = FALSE) {
  p <- prep$p
  a_lb <- pmin(0, prep$a_exp_lower)
  a_ub <- pmax(0, prep$a_exp_upper)
  r_ub <- pmax(abs(prep$a_exp_lower), abs(prep$a_exp_upper), 0)
  z_ub <- pmax(abs(prep$omega - prep$a_exp_lower),
               abs(prep$omega - prep$a_exp_upper),
               abs(prep$omega))
  lb <- rep(0, idx$nvar)
  ub <- rep(Inf, idx$nvar)
  lb[idx$idx_x] <- 0
  ub[idx$idx_x] <- prep$x_max
  lb[idx$idx_a] <- a_lb
  ub[idx$idx_a] <- a_ub
  lb[idx$idx_r] <- 0
  ub[idx$idx_r] <- r_ub
  lb[idx$idx_z] <- 0
  ub[idx$idx_z] <- z_ub
  if (prep$use_sets) {
    lb[idx$idx_y] <- 0
    ub[idx$idx_y] <- 1
  }
  if (isTRUE(include_t)) {
    lb[idx$idx_t] <- 0
    ub[idx$idx_t] <- Inf
  }
  list(lb = lb, ub = ub)
}

.add_general_a_constraints <- function(prep, idx) {
  prep$current_idx_y <- idx$idx_y
  A <- .empty_sparse(0L, idx$nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .add_x_constraints(A, rhs, sense, idx$nvar, idx$idx_x, prep)
  tmp <- .add_a_constraints(tmp$A, tmp$rhs, tmp$sense, idx$nvar,
                            idx$idx_x, idx$idx_a, prep)
  .add_abs_a_constraints(tmp$A, tmp$rhs, tmp$sense, idx$nvar,
                         idx$idx_a, idx$idx_r, idx$idx_z, prep)
}

.solution_general_a <- function(prep, idx, result) {
  sol <- result$x
  x <- as.integer(round(sol[idx$idx_x]))
  a_exp <- pmax(prep$a_exp_lower * x, pmin(prep$a_exp_upper * x,
                                           as.numeric(sol[idx$idx_a])))
  a_obs <- prep$omega - a_exp
  gamma <- .implied_gamma(a_exp, prep$omega)
  list(
    x = x,
    s = gamma,
    gamma = gamma,
    a_exp = a_exp,
    a_obs = a_obs,
    alpha = prep$compute_alpha(a_exp),
    beta = prep$compute_beta(a_exp),
    result = result
  )
}

.compute_alpha_with_n <- function(prep, weight, n_opt, tol = 1e-10) {
  weight <- as.numeric(weight)
  n_opt <- as.numeric(n_opt)
  if (identical(prep$weight_mode, "a")) {
    a_exp <- weight
  } else {
    a_exp <- prep$omega * weight
  }
  a_obs <- prep$omega - a_exp
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
  as.numeric(t(a_obs) %*% prep$Sigma_obs %*% a_obs + sum(exp_terms))
}

.compute_allocation_with_floor <- function(prep, weight, x = NULL,
                                           n_min = prep$n_min,
                                           tol = 1e-10) {
  n_min <- as.numeric(n_min)
  if (!is.numeric(n_min) || length(n_min) != 1L || !is.finite(n_min) ||
      n_min < 0) {
    stop("n_min must be a nonnegative finite scalar.", call. = FALSE)
  }
  if (n_min == 0) {
    return(prep$compute_allocation(weight, x = x))
  }

  weight <- as.numeric(weight)
  if (identical(prep$weight_mode, "a")) {
    a_exp <- weight
  } else {
    a_exp <- prep$omega * weight
  }
  if (is.null(x)) {
    selected <- a_exp != 0
  } else {
    x <- as.numeric(x)
    if (length(x) != prep$p || any(!is.finite(x))) {
      stop("x must be a finite vector with length equal to length(omega).",
           call. = FALSE)
    }
    selected <- x > 0.5
  }

  n_opt <- rep(0, prep$p)
  if (!any(selected)) {
    return(n_opt)
  }

  floor_cost <- sum(prep$costs[selected] * n_min)
  budget_scale <- max(prep$n_total, floor_cost)
  budget_tol <- tol * budget_scale
  if (floor_cost > prep$n_total + budget_tol) {
    stop("n_min is infeasible for the selected experiments.", call. = FALSE)
  }
  n_opt[selected] <- n_min
  if (abs(floor_cost - prep$n_total) <= budget_tol) {
    return(n_opt)
  }

  score <- abs(a_exp) * prep$v
  positive <- selected & score > 0
  if (!any(positive)) {
    j <- which(selected)[which.min(prep$costs[selected])]
    n_opt[j] <- n_opt[j] + (prep$n_total - floor_cost) / prep$costs[j]
    return(n_opt)
  }

  cost_at <- function(lambda) {
    n <- n_opt
    n[positive] <- pmax(n_min, score[positive] /
                          sqrt(lambda * prep$costs[positive]))
    sum(prep$costs * n)
  }

  hi <- 1
  while (cost_at(hi) > prep$n_total && hi < .Machine$double.xmax / 2) {
    hi <- hi * 2
  }
  if (!is.finite(hi) || cost_at(hi) > prep$n_total + budget_tol) {
    stop("Could not find a feasible floor-constrained allocation.",
         call. = FALSE)
  }
  lo <- 0
  for (iter in seq_len(100L)) {
    mid <- (lo + hi) / 2
    if (cost_at(mid) > prep$n_total) {
      lo <- mid
    } else {
      hi <- mid
    }
  }

  lambda <- hi
  n_opt[positive] <- pmax(n_min, score[positive] /
                            sqrt(lambda * prep$costs[positive]))
  total_cost <- sum(prep$costs * n_opt)
  if (abs(total_cost - prep$n_total) > budget_tol) {
    j <- which(positive)[which.max(score[positive])]
    n_opt[j] <- n_opt[j] + (prep$n_total - total_cost) / prep$costs[j]
  }
  pmax(0, n_opt)
}

.validate_explicit_allocation <- function(prep, n_opt, x = NULL, n_min = 0,
                                          tol = 1e-8) {
  n_opt <- as.numeric(n_opt)
  if (length(n_opt) != prep$p || any(!is.finite(n_opt))) {
    stop("n_opt must be a finite numeric vector with length equal to length(omega).",
         call. = FALSE)
  }
  if (any(n_opt < -tol)) {
    stop("n_opt must be nonnegative.", call. = FALSE)
  }
  n_opt <- pmax(0, n_opt)
  total_cost <- sum(prep$costs * n_opt)
  budget_scale <- max(prep$n_total, total_cost)
  budget_tol <- tol * budget_scale
  if (total_cost > prep$n_total + budget_tol) {
    stop("n_opt exceeds the total sample budget.", call. = FALSE)
  }
  if (!is.null(x)) {
    x <- as.numeric(x)
    if (length(x) != prep$p || any(!is.finite(x)) ||
        any(!(x %in% base::c(0, 1)))) {
      stop("x must contain one binary entry per experiment.", call. = FALSE)
    }
    selected <- x > 0.5
    if (any(n_opt[!selected] > tol)) {
      stop("n_opt must be zero for unselected experiments.", call. = FALSE)
    }
    n_min <- as.numeric(n_min)
    if (length(n_min) != 1L || !is.finite(n_min) || n_min < 0) {
      stop("n_min must be a nonnegative finite scalar.", call. = FALSE)
    }
    if (n_min > 0 && any(n_opt[selected] < n_min - tol)) {
      stop("n_opt violates n_min for selected experiments.", call. = FALSE)
    }
    if (any(selected) && abs(total_cost - prep$n_total) > budget_tol) {
      stop("n_opt must exhaust n_total when at least one experiment is selected.",
           call. = FALSE)
    }
  }
  n_opt
}


.repair_explicit_allocation <- function(prep, n_opt, x, n_min = 0,
                                        tol = 1e-8) {
  tol <- .validate_tolerance(tol)
  n_opt <- as.numeric(n_opt)
  x <- as.numeric(x)
  if (length(n_opt) != prep$p || length(x) != prep$p ||
      any(!is.finite(n_opt)) || any(!is.finite(x))) {
    stop("Internal allocation and selection vectors have invalid dimensions.",
         call. = FALSE)
  }
  selected <- x > 0.5
  n_opt[!selected] <- 0
  n_opt[selected] <- pmax(n_min, n_opt[selected])
  if (!any(selected)) return(n_opt)

  spent <- sum(prep$costs * n_opt)
  residual <- prep$n_total - spent
  budget_tol <- tol * max(prep$n_total, abs(spent))
  if (abs(residual) > budget_tol) {
    stop(sprintf(
      "Optimizer allocation violates the budget by %.6g, beyond tolerance %.6g.",
      residual, budget_tol
    ), call. = FALSE)
  }
  if (residual != 0) {
    candidates <- which(selected)
    if (residual < 0) {
      slack <- n_opt[candidates] - n_min
      candidates <- candidates[slack > 0]
      if (length(candidates) == 0L) {
        stop("Cannot repair the optimizer allocation without violating n_min.",
             call. = FALSE)
      }
      j <- candidates[which.max(n_opt[candidates] - n_min)]
    } else {
      headroom <- prep$n_total / prep$costs[candidates] - n_opt[candidates]
      candidates <- candidates[headroom > 0]
      if (length(candidates) == 0L) {
        stop("Cannot repair the optimizer allocation within its upper bounds.",
             call. = FALSE)
      }
      j <- candidates[which.max(
        prep$n_total / prep$costs[candidates] - n_opt[candidates]
      )]
    }
    n_opt[j] <- n_opt[j] + residual / prep$costs[j]
  }
  .validate_explicit_allocation(prep, n_opt, x = as.integer(selected),
                                n_min = n_min, tol = max(tol, 1e-10))
}

.add_obs_abs_constraints <- function(A_mat, rhs_vec, sense_vec, nvar, idx_a,
                                     idx_z, omega) {
  p <- length(omega)
  for (j in seq_len(p)) {
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- 1
    row[1L, idx_z[j]] <- -1
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, omega[j])
    sense_vec <- base::c(sense_vec, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_a[j]] <- -1
    row[1L, idx_z[j]] <- -1
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, -omega[j])
    sense_vec <- base::c(sense_vec, "<")
  }
  list(A = A_mat, rhs = rhs_vec, sense = sense_vec)
}

.alloc_gurobi_indices <- function(prep, include_t = FALSE) {
  p <- prep$p
  offset <- 0L
  idx_x <- seq_len(p)
  offset <- offset + p
  idx_w <- offset + seq_len(p)
  offset <- offset + p
  if (identical(prep$weight_mode, "a")) {
    idx_z <- offset + seq_len(p)
    offset <- offset + p
  } else {
    idx_z <- integer(0)
  }
  idx_n <- offset + seq_len(p)
  offset <- offset + p
  idx_u <- offset + seq_len(p)
  offset <- offset + p
  if (prep$use_sets) {
    idx_y <- offset + seq_len(prep$Msets)
    offset <- offset + prep$Msets
  } else {
    idx_y <- integer(0)
  }
  if (isTRUE(include_t)) {
    idx_t <- offset + 1L
    offset <- offset + 1L
  } else {
    idx_t <- integer(0)
  }
  list(nvar = offset, idx_x = idx_x, idx_w = idx_w, idx_z = idx_z,
       idx_n = idx_n, idx_u = idx_u, idx_y = idx_y, idx_t = idx_t)
}

.alloc_gurobi_bounds <- function(prep, idx, include_t = FALSE) {
  lb <- rep(0, idx$nvar)
  ub <- rep(Inf, idx$nvar)
  lb[idx$idx_x] <- 0
  ub[idx$idx_x] <- prep$x_max
  if (identical(prep$weight_mode, "a")) {
    lb[idx$idx_w] <- pmin(0, prep$a_exp_lower)
    ub[idx$idx_w] <- pmax(0, prep$a_exp_upper)
    z_ub <- pmax(abs(prep$omega - prep$a_exp_lower),
                 abs(prep$omega - prep$a_exp_upper),
                 abs(prep$omega))
    lb[idx$idx_z] <- 0
    ub[idx$idx_z] <- z_ub
  } else {
    lb[idx$idx_w] <- 0
    ub[idx$idx_w] <- prep$gamma_upper
  }
  lb[idx$idx_n] <- 0
  ub[idx$idx_n] <- prep$n_total / prep$costs
  lb[idx$idx_u] <- 0
  ub[idx$idx_u] <- Inf
  if (prep$use_sets) {
    lb[idx$idx_y] <- 0
    ub[idx$idx_y] <- 1
  }
  if (isTRUE(include_t)) {
    lb[idx$idx_t] <- 0
    ub[idx$idx_t] <- Inf
  }
  list(lb = lb, ub = ub)
}

.add_allocation_constraints <- function(A_mat, rhs_vec, sense_vec, nvar,
                                        idx_x, idx_n, prep) {
  p <- prep$p
  for (j in seq_len(p)) {
    row <- .empty_sparse(1L, nvar)
    row[1L, idx_n[j]] <- 1
    row[1L, idx_x[j]] <- -prep$n_total / prep$costs[j]
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, 0)
    sense_vec <- base::c(sense_vec, "<")

    row <- .empty_sparse(1L, nvar)
    row[1L, idx_n[j]] <- -1
    row[1L, idx_x[j]] <- prep$n_min
    A_mat <- rbind(A_mat, row)
    rhs_vec <- base::c(rhs_vec, 0)
    sense_vec <- base::c(sense_vec, "<")
  }
  row <- .empty_sparse(1L, nvar)
  row[1L, idx_n] <- prep$costs
  A_mat <- rbind(A_mat, row)
  rhs_vec <- base::c(rhs_vec, prep$n_total)
  sense_vec <- base::c(sense_vec, "=")
  list(A = A_mat, rhs = rhs_vec, sense = sense_vec)
}

.add_explicit_design_constraints <- function(prep, idx) {
  prep$current_idx_y <- idx$idx_y
  A <- .empty_sparse(0L, idx$nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .add_x_constraints(A, rhs, sense, idx$nvar, idx$idx_x, prep)
  if (identical(prep$weight_mode, "a")) {
    tmp <- .add_a_constraints(tmp$A, tmp$rhs, tmp$sense, idx$nvar,
                              idx$idx_x, idx$idx_w, prep)
    tmp <- .add_obs_abs_constraints(tmp$A, tmp$rhs, tmp$sense, idx$nvar,
                                    idx$idx_w, idx$idx_z, prep$omega)
  } else {
    tmp <- .add_s_constraints(tmp$A, tmp$rhs, tmp$sense, idx$nvar,
                              idx$idx_x, idx$idx_w, prep)
  }
  .add_allocation_constraints(tmp$A, tmp$rhs, tmp$sense, idx$nvar,
                              idx$idx_x, idx$idx_n, prep)
}

.add_explicit_variance_quadcons <- function(quadcon, prep, idx) {
  for (j in seq_len(prep$p)) {
    Qc <- .empty_sparse(idx$nvar, idx$nvar)
    coeff <- prep$v2[j]
    if (identical(prep$weight_mode, "gamma")) {
      coeff <- coeff * prep$omega[j]^2
    }
    if (coeff > 0) {
      Qc[idx$idx_w[j], idx$idx_w[j]] <- coeff
    }
    Qc[idx$idx_u[j], idx$idx_n[j]] <- -0.5
    Qc[idx$idx_n[j], idx$idx_u[j]] <- -0.5
    quadcon[[length(quadcon) + 1L]] <- list(
      Qc = Qc,
      q = rep(0, idx$nvar),
      rhs = 0,
      sense = "<",
      name = paste0("exp_var_", j)
    )
  }
  quadcon
}

.solution_explicit_alloc <- function(prep, idx, result) {
  sol <- result$x
  x <- as.integer(round(sol[idx$idx_x]))
  n_opt <- rep(0, prep$p)
  n_opt[x == 1L] <- pmax(0, as.numeric(sol[idx$idx_n[x == 1L]]))
  n_opt <- .repair_explicit_allocation(
    prep, n_opt, x = x, n_min = prep$n_min, tol = 1e-8
  )
  if (identical(prep$weight_mode, "a")) {
    weight <- pmax(prep$a_exp_lower * x, pmin(prep$a_exp_upper * x,
                                              as.numeric(sol[idx$idx_w])))
    a_exp <- weight
    gamma <- .implied_gamma(a_exp, prep$omega)
    s <- gamma
  } else {
    weight <- pmax(prep$gamma_lower * x, pmin(prep$gamma_upper * x,
                                              as.numeric(sol[idx$idx_w])))
    s <- weight
    gamma <- rep(0, prep$p)
    gamma[x == 1L] <- s[x == 1L]
    a_exp <- prep$omega * s
  }
  a_exp[x == 0L] <- 0
  a_obs <- prep$omega - a_exp
  list(
    x = x,
    s = s,
    gamma = gamma,
    a_exp = a_exp,
    a_obs = a_obs,
    n = n_opt,
    u = pmax(0, as.numeric(sol[idx$idx_u])),
    alpha = .compute_alpha_with_n(prep, weight, n_opt),
    beta = prep$compute_beta(weight),
    result = result
  )
}

.solve_alpha_oracle_explicit_alloc <- function(prep, gurobi_params = list(),
                                               status_label = "variance oracle") {
  .minimax_require("gurobi")
  idx <- .alloc_gurobi_indices(prep, include_t = FALSE)
  bounds <- .alloc_gurobi_bounds(prep, idx, include_t = FALSE)

  model <- list()
  model$modelname <- "alpha_star_explicit_allocation"
  model$modelsense <- "min"
  model$vtype <- rep("C", idx$nvar)
  model$vtype[idx$idx_x] <- "B"
  if (prep$use_sets) model$vtype[idx$idx_y] <- "B"
  model$lb <- bounds$lb
  model$ub <- bounds$ub
  model$obj <- rep(0, idx$nvar)
  model$obj[idx$idx_u] <- 1
  if (identical(prep$weight_mode, "a")) {
    model$obj[idx$idx_w] <- prep$q_alpha_a
    model$Q <- .block_quadratic(prep$M_alpha_a, idx$idx_w, idx$nvar)
  } else {
    model$obj[idx$idx_w] <- prep$q_alpha
    model$Q <- .block_quadratic(prep$Qobs, idx$idx_w, idx$nvar)
  }

  tmp <- .add_explicit_design_constraints(prep, idx)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense
  model$quadcon <- .add_explicit_variance_quadcons(list(), prep, idx)

  result <- .gurobi_solve(model, params = gurobi_params)
  .require_gurobi_optimal(result, status_label)
  .solution_explicit_alloc(prep, idx, result)
}

.solve_beta_oracle_explicit_alloc <- function(prep, gurobi_params = list()) {
  .minimax_require("gurobi")
  idx <- .alloc_gurobi_indices(prep, include_t = FALSE)
  bounds <- .alloc_gurobi_bounds(prep, idx, include_t = FALSE)

  model <- list()
  model$modelname <- "beta_star_explicit_allocation"
  model$modelsense <- "min"
  model$vtype <- rep("C", idx$nvar)
  model$vtype[idx$idx_x] <- "B"
  if (prep$use_sets) model$vtype[idx$idx_y] <- "B"
  model$lb <- bounds$lb
  model$ub <- bounds$ub
  model$obj <- rep(0, idx$nvar)
  if (identical(prep$weight_mode, "a")) {
    model$Q <- .block_quadratic(prep$M_beta_a, idx$idx_z, idx$nvar)
  } else {
    model$obj[idx$idx_w] <- prep$q_beta
    model$Q <- .block_quadratic(prep$M_beta, idx$idx_w, idx$nvar)
  }

  tmp <- .add_explicit_design_constraints(prep, idx)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  result <- .gurobi_solve(model, params = gurobi_params)
  .require_gurobi_optimal(result, "bias oracle")
  .solution_explicit_alloc(prep, idx, result)
}

.solve_minimax_explicit_alloc <- function(prep, alpha_star, beta_star,
                                          gurobi_params = list()) {
  .minimax_require("gurobi")
  idx <- .alloc_gurobi_indices(prep, include_t = TRUE)
  bounds <- .alloc_gurobi_bounds(prep, idx, include_t = TRUE)

  model <- list()
  model$modelname <- "minimax_design_explicit_allocation"
  model$modelsense <- "min"
  model$vtype <- rep("C", idx$nvar)
  model$vtype[idx$idx_x] <- "B"
  if (prep$use_sets) model$vtype[idx$idx_y] <- "B"
  model$lb <- bounds$lb
  model$ub <- bounds$ub
  model$obj <- rep(0, idx$nvar)
  model$obj[idx$idx_t] <- 1

  tmp <- .add_explicit_design_constraints(prep, idx)
  if (.beta_star_is_zero(prep, beta_star) && any(prep$bias_weights > 0)) {
    row <- .empty_sparse(1L, idx$nvar)
    if (identical(prep$weight_mode, "a")) {
      row[1L, idx$idx_z] <- prep$bias_weights
      rhs_zero <- 0
    } else {
      d_normalized <- if (prep$beta_linear_scale > 0) {
        prep$d / prep$beta_linear_scale
      } else {
        prep$d
      }
      row[1L, idx$idx_w] <- d_normalized
      rhs_zero <- sum(d_normalized)
    }
    tmp$A <- rbind(tmp$A, row)
    tmp$rhs <- base::c(tmp$rhs, rhs_zero)
    tmp$sense <- base::c(tmp$sense, "=")
  }
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  q_alpha <- rep(0, idx$nvar)
  q_alpha[idx$idx_u] <- 1
  q_alpha[idx$idx_t] <- -alpha_star
  if (identical(prep$weight_mode, "a")) {
    q_alpha[idx$idx_w] <- prep$q_alpha_a
    rhs_alpha <- -prep$const_alpha_a
    M_alpha <- prep$M_alpha_a
  } else {
    q_alpha[idx$idx_w] <- prep$q_alpha
    rhs_alpha <- -prep$const_alpha
    M_alpha <- prep$Qobs
  }
  scale_alpha <- .quad_scale(M_alpha, q_alpha, rhs_alpha)
  quadcon <- list(
    list(
      Qc = .block_quadratic(M_alpha, idx$idx_w, idx$nvar, scale_alpha),
      q = scale_alpha * q_alpha,
      rhs = scale_alpha * rhs_alpha,
      sense = "<",
      name = "alpha_ratio"
    )
  )
  quadcon <- .add_explicit_variance_quadcons(quadcon, prep, idx)

  if (!.beta_star_is_zero(prep, beta_star)) {
    q_beta <- rep(0, idx$nvar)
    q_beta[idx$idx_t] <- -beta_star
    if (identical(prep$weight_mode, "a")) {
      M_beta <- prep$M_beta_a
      idx_beta <- idx$idx_z
      rhs_beta <- 0
    } else {
      M_beta <- prep$M_beta
      idx_beta <- idx$idx_w
      q_beta[idx$idx_w] <- prep$q_beta
      rhs_beta <- -prep$const_beta
    }
    scale_beta <- .quad_scale(M_beta, q_beta, rhs_beta)
    quadcon[[length(quadcon) + 1L]] <- list(
      Qc = .block_quadratic(M_beta, idx_beta, idx$nvar, scale_beta),
      q = scale_beta * q_beta,
      rhs = scale_beta * rhs_beta,
      sense = "<",
      name = "beta_ratio"
    )
  }
  model$quadcon <- quadcon

  result <- .gurobi_solve(model, params = gurobi_params)
  .require_gurobi_optimal(result, "minimax program")

  out <- .solution_explicit_alloc(prep, idx, result)
  if (.beta_star_is_zero(prep, beta_star)) {
    out$beta <- .clamp_relative_zero(
      out$beta, prep$beta_scale, prep$zero_tol, "beta"
    )
  }
  out$alpha_ratio <- .safe_ratio(out$alpha, alpha_star)
  out$beta_ratio <- .safe_ratio(out$beta, beta_star)
  out$regret <- max(out$alpha_ratio, out$beta_ratio)
  epi <- .validate_regret_epigraph(
    out$regret, as.numeric(result$x[idx$idx_t])
  )
  out$t <- epi$t
  out$solver_t <- epi$solver_t
  out$result$validated_ratio_gap <- epi$validated_gap
  out
}

.fixed_set_a_bounds <- function(prep, S, include_z = FALSE,
                                include_t = FALSE) {
  p <- prep$p
  offset <- 0L
  idx_a <- offset + seq_len(p)
  offset <- offset + p
  idx_r <- offset + seq_len(p)
  offset <- offset + p
  if (isTRUE(include_z)) {
    idx_z <- offset + seq_len(p)
    offset <- offset + p
  } else {
    idx_z <- integer(0)
  }
  if (isTRUE(include_t)) {
    idx_t <- offset + 1L
    offset <- offset + 1L
  } else {
    idx_t <- integer(0)
  }

  nvar <- offset
  lb <- ub <- rep(0, nvar)
  if (length(S) > 0L) {
    lb[idx_a[S]] <- prep$a_exp_lower[S]
    ub[idx_a[S]] <- prep$a_exp_upper[S]
    ub[idx_r[S]] <- pmax(abs(prep$a_exp_lower[S]),
                         abs(prep$a_exp_upper[S]))
  }
  if (isTRUE(include_z)) {
    z_ub <- pmax(abs(prep$omega - prep$a_exp_lower),
                 abs(prep$omega - prep$a_exp_upper),
                 abs(prep$omega))
    lb[idx_z] <- 0
    ub[idx_z] <- z_ub
  }
  if (isTRUE(include_t)) {
    lb[idx_t] <- 1
    ub[idx_t] <- Inf
  }

  list(nvar = nvar, idx_a = idx_a, idx_r = idx_r, idx_z = idx_z,
       idx_t = idx_t, lb = lb, ub = ub)
}

.fixed_set_a_abs_constraints <- function(prep, S, idx, include_z = FALSE) {
  A <- .empty_sparse(0L, idx$nvar)
  rhs <- numeric(0)
  sense <- character(0)

  for (j in S) {
    row <- .empty_sparse(1L, idx$nvar)
    row[1L, idx$idx_a[j]] <- 1
    row[1L, idx$idx_r[j]] <- -1
    A <- rbind(A, row)
    rhs <- base::c(rhs, 0)
    sense <- base::c(sense, "<")

    row <- .empty_sparse(1L, idx$nvar)
    row[1L, idx$idx_a[j]] <- -1
    row[1L, idx$idx_r[j]] <- -1
    A <- rbind(A, row)
    rhs <- base::c(rhs, 0)
    sense <- base::c(sense, "<")
  }

  if (isTRUE(include_z)) {
    for (j in seq_len(prep$p)) {
      row <- .empty_sparse(1L, idx$nvar)
      row[1L, idx$idx_a[j]] <- 1
      row[1L, idx$idx_z[j]] <- -1
      A <- rbind(A, row)
      rhs <- base::c(rhs, prep$omega[j])
      sense <- base::c(sense, "<")

      row <- .empty_sparse(1L, idx$nvar)
      row[1L, idx$idx_a[j]] <- -1
      row[1L, idx$idx_z[j]] <- -1
      A <- rbind(A, row)
      rhs <- base::c(rhs, -prep$omega[j])
      sense <- base::c(sense, "<")
    }
  }

  list(A = A, rhs = rhs, sense = sense)
}

.fixed_set_quadratic_blocks <- function(prep, S, idx) {
  if (length(S) == 0L) {
    return(list())
  }
  list(
    list(M = prep$M_alpha_a[S, S, drop = FALSE], idx = idx$idx_a[S]),
    list(M = prep$M_exp_a[S, S, drop = FALSE], idx = idx$idx_r[S])
  )
}

.fixed_set_gurobi_params <- function(gurobi_params, qcp = FALSE) {
  params <- gurobi_params
  if (is.null(params$FeasibilityTol)) params$FeasibilityTol <- 1e-9
  if (is.null(params$OptimalityTol)) params$OptimalityTol <- 1e-9
  if (is.null(params$NumericFocus)) params$NumericFocus <- 2
  if (isTRUE(qcp) && is.null(params$BarConvTol)) params$BarConvTol <- 1e-10
  if (!isTRUE(qcp)) {
    params$Method <- 0
    params$BarHomogeneous <- NULL
  }
  params
}

.solve_alpha_oracle_a_fixed_set <- function(prep, S, gurobi_params = list(),
                                            validation_tol = 1e-9) {
  idx <- .fixed_set_a_bounds(prep, S, include_z = FALSE,
                             include_t = FALSE)
  model <- list(
    modelname = "alpha_star_general_a_fixed_set",
    modelsense = "min",
    vtype = rep("C", idx$nvar),
    lb = idx$lb,
    ub = idx$ub,
    obj = rep(0, idx$nvar)
  )
  if (length(S) > 0L) {
    model$obj[idx$idx_a[S]] <- prep$q_alpha_a[S]
  }
  model$Q <- .block_quadratic_sum(
    .fixed_set_quadratic_blocks(prep, S, idx),
    idx$nvar
  )
  constraints <- .fixed_set_a_abs_constraints(
    prep, S, idx, include_z = FALSE
  )
  model$A <- constraints$A
  model$rhs <- constraints$rhs
  model$sense <- constraints$sense

  result <- .gurobi_solve(
    model,
    params = .fixed_set_gurobi_params(gurobi_params, qcp = FALSE)
  )
  if (is.null(result$x) || !identical(result$status, "OPTIMAL")) {
    return(list(error = if (is.null(result$status)) "NO_SOLUTION" else result$status,
                result = result))
  }

  a_exp <- rep(0, prep$p)
  if (length(S) > 0L) {
    a_exp[S] <- pmax(
      prep$a_exp_lower[S],
      pmin(prep$a_exp_upper[S], as.numeric(result$x[idx$idx_a[S]]))
    )
  }
  alpha <- prep$compute_alpha(a_exp)
  objective_alpha <- as.numeric(result$objval) + prep$const_alpha_a
  objective_gap <- abs(objective_alpha - alpha)
  validation_scale <- max(base::c(prep$alpha_scale, abs(alpha),
                                 abs(objective_alpha), 0))
  gap_tol <- if (validation_scale == 0) 0 else validation_tol * validation_scale
  if (!is.finite(alpha) || objective_gap > gap_tol) {
    return(list(error = "NATURAL_OBJECTIVE_VALIDATION_FAILED",
                result = result, objective_gap = objective_gap))
  }

  x <- .set_to_x(S, prep$p)
  gamma <- .implied_gamma(a_exp, prep$omega)
  result$validated_objective_gap <- objective_gap
  list(
    x = as.integer(x),
    s = gamma,
    gamma = gamma,
    a_exp = a_exp,
    a_obs = prep$omega - a_exp,
    alpha = alpha,
    beta = prep$compute_beta(a_exp),
    result = result
  )
}

.solve_alpha_oracle_a_sets <- function(prep, gurobi_params = list()) {
  sets <- prep$sets
  candidates <- vector("list", length(sets))
  diagnostics <- data.frame(
    candidate = seq_along(sets),
    set = vapply(sets, function(S) paste(S, collapse = ","), character(1)),
    status = NA_character_,
    alpha = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(sets)) {
    candidate <- .solve_alpha_oracle_a_fixed_set(
      prep, sets[[i]], gurobi_params = gurobi_params
    )
    candidates[[i]] <- candidate
    if (!is.null(candidate$error)) {
      diagnostics$status[i] <- candidate$error
    } else {
      diagnostics$status[i] <- candidate$result$status
      diagnostics$alpha[i] <- candidate$alpha
    }
  }
  failed <- which(diagnostics$status != "OPTIMAL" |
                    !is.finite(diagnostics$alpha))
  if (length(failed) > 0L) {
    stop(
      "A fixed-set variance oracle failed, so global optimality cannot be ",
      "certified. Candidate status(es): ",
      paste(diagnostics$status[failed], collapse = ", "),
      call. = FALSE
    )
  }

  best_i <- which.min(diagnostics$alpha)
  best <- candidates[[best_i]]
  best$result$solver <- "gurobi_fixed_set_enumeration"
  best$result$enumerated_sets <- length(sets)
  best$result$candidate_diagnostics <- diagnostics
  best
}

.solve_minimax_a_fixed_set <- function(prep, S, alpha_star, beta_star,
                                       gurobi_params = list(),
                                       validation_tol = 1e-6) {
  idx <- .fixed_set_a_bounds(prep, S, include_z = TRUE, include_t = TRUE)
  constraints <- .fixed_set_a_abs_constraints(
    prep, S, idx, include_z = TRUE
  )

  if (.beta_star_is_zero(prep, beta_star) && any(prep$bias_weights > 0)) {
    row <- .empty_sparse(1L, idx$nvar)
    row[1L, idx$idx_z] <- if (prep$beta_linear_scale > 0) {
      prep$bias_weights / prep$beta_linear_scale
    } else {
      prep$bias_weights
    }
    constraints$A <- rbind(constraints$A, row)
    constraints$rhs <- base::c(constraints$rhs, 0)
    constraints$sense <- base::c(constraints$sense, "=")
  }

  q_alpha <- rep(0, idx$nvar)
  if (length(S) > 0L) {
    q_alpha[idx$idx_a[S]] <- prep$q_alpha_a[S]
  }
  q_alpha[idx$idx_t] <- -alpha_star
  rhs_alpha <- -prep$const_alpha_a
  alpha_blocks <- .fixed_set_quadratic_blocks(prep, S, idx)
  alpha_matrices <- lapply(alpha_blocks, `[[`, "M")
  scale_alpha <- .quad_scale_multi(alpha_matrices, q_alpha, rhs_alpha)
  alpha_blocks <- lapply(alpha_blocks, function(block) {
    block$scale <- scale_alpha
    block
  })
  quadcon <- list(
    list(
      Qc = .block_quadratic_sum(alpha_blocks, idx$nvar),
      q = scale_alpha * q_alpha,
      rhs = scale_alpha * rhs_alpha,
      sense = "<",
      name = "alpha_ratio"
    )
  )

  if (!.beta_star_is_zero(prep, beta_star)) {
    q_beta <- rep(0, idx$nvar)
    q_beta[idx$idx_t] <- -beta_star
    scale_beta <- .quad_scale(prep$M_beta_a, q_beta, 0)
    quadcon[[length(quadcon) + 1L]] <- list(
      Qc = .block_quadratic(
        prep$M_beta_a, idx$idx_z, idx$nvar, scale_beta
      ),
      q = scale_beta * q_beta,
      rhs = 0,
      sense = "<",
      name = "beta_ratio"
    )
  }

  model <- list(
    modelname = "minimax_design_general_a_fixed_set",
    modelsense = "min",
    vtype = rep("C", idx$nvar),
    lb = idx$lb,
    ub = idx$ub,
    obj = replace(rep(0, idx$nvar), idx$idx_t, 1),
    A = constraints$A,
    rhs = constraints$rhs,
    sense = constraints$sense,
    quadcon = quadcon
  )
  result <- tryCatch(
    .gurobi_solve(
      model,
      params = .fixed_set_gurobi_params(gurobi_params, qcp = TRUE)
    ),
    error = function(e) {
      structure(list(error = conditionMessage(e)), class = "fixed_set_error")
    }
  )
  if (inherits(result, "fixed_set_error")) {
    return(list(error = result$error))
  }
  status <- if (is.null(result$status)) "NO_SOLUTION" else result$status
  if (is.null(result$x) || !identical(status, "OPTIMAL")) {
    return(list(error = status, result = result))
  }

  a_exp <- rep(0, prep$p)
  if (length(S) > 0L) {
    a_exp[S] <- pmax(
      prep$a_exp_lower[S],
      pmin(prep$a_exp_upper[S], as.numeric(result$x[idx$idx_a[S]]))
    )
  }
  alpha <- prep$compute_alpha(a_exp)
  beta <- prep$compute_beta(a_exp)
  if (.beta_star_is_zero(prep, beta_star)) {
    beta <- .clamp_relative_zero(
      beta, prep$beta_scale, prep$zero_tol, "beta"
    )
  }
  alpha_ratio <- .safe_ratio(alpha, alpha_star)
  beta_ratio <- .safe_ratio(beta, beta_star)
  regret <- max(alpha_ratio, beta_ratio)
  solver_t <- as.numeric(result$x[idx$idx_t])
  validated_gap <- regret - solver_t
  gap_tol <- validation_tol * max(1, abs(regret), abs(solver_t))
  if (!is.finite(regret) || !is.finite(solver_t) || validated_gap > gap_tol) {
    return(list(error = "NATURAL_RATIO_VALIDATION_FAILED",
                result = result, validated_gap = validated_gap,
                validation_tolerance = gap_tol))
  }

  x <- .set_to_x(S, prep$p)
  gamma <- .implied_gamma(a_exp, prep$omega)
  result$solver_t <- solver_t
  result$validated_regret <- regret
  result$validated_ratio_gap <- validated_gap
  list(
    x = as.integer(x),
    s = gamma,
    gamma = gamma,
    a_exp = a_exp,
    a_obs = prep$omega - a_exp,
    alpha = alpha,
    beta = beta,
    alpha_ratio = alpha_ratio,
    beta_ratio = beta_ratio,
    regret = regret,
    t = max(regret, solver_t),
    solver_t = solver_t,
    result = result
  )
}

.solve_minimax_a_sets <- function(prep, alpha_star, beta_star,
                                  gurobi_params = list()) {
  sets <- prep$sets
  candidates <- vector("list", length(sets))
  diagnostics <- data.frame(
    candidate = seq_along(sets),
    set = vapply(sets, function(S) paste(S, collapse = ","), character(1)),
    status = NA_character_,
    regret = NA_real_,
    solver_t = NA_real_,
    validated_gap = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_along(sets)) {
    candidate <- .solve_minimax_a_fixed_set(
      prep, sets[[i]], alpha_star = alpha_star, beta_star = beta_star,
      gurobi_params = gurobi_params
    )
    candidates[[i]] <- candidate
    if (!is.null(candidate$error)) {
      diagnostics$status[i] <- candidate$error
    } else {
      diagnostics$status[i] <- candidate$result$status
      diagnostics$regret[i] <- candidate$regret
      diagnostics$solver_t[i] <- candidate$solver_t
      diagnostics$validated_gap[i] <-
        candidate$result$validated_ratio_gap
    }
  }

  zero_bias_infeasible <- .beta_star_is_zero(prep, beta_star) &
    diagnostics$status %in% base::c("INFEASIBLE", "INF_OR_UNBD")
  failed <- which(
    !zero_bias_infeasible &
      (diagnostics$status != "OPTIMAL" | !is.finite(diagnostics$regret))
  )
  if (length(failed) > 0L) {
    stop(
      "A fixed-set minimax solve failed, so global optimality cannot be ",
      "certified. Candidate status(es): ",
      paste(diagnostics$status[failed], collapse = ", "),
      call. = FALSE
    )
  }

  feasible <- which(diagnostics$status == "OPTIMAL" &
                      is.finite(diagnostics$regret))
  if (length(feasible) == 0L) {
    stop("No fixed experiment set has finite minimax regret.", call. = FALSE)
  }
  best_i <- feasible[which.min(diagnostics$regret[feasible])]
  best <- candidates[[best_i]]
  best$result$solver <- "gurobi_fixed_set_enumeration"
  best$result$enumerated_sets <- length(sets)
  best$result$candidate_diagnostics <- diagnostics
  best
}

.solve_alpha_oracle_a <- function(prep, gurobi_params = list(),
                                  status_label = "variance oracle") {
  .minimax_require("gurobi")
  idx <- .a_gurobi_indices(prep, include_t = FALSE)
  bounds <- .a_gurobi_bounds(prep, idx, include_t = FALSE)

  model <- list()
  model$modelname <- "alpha_star_general_a"
  model$modelsense <- "min"
  model$vtype <- base::c(rep("B", prep$p), rep("C", 3L * prep$p),
                   if (prep$use_sets) rep("B", prep$Msets))
  model$lb <- bounds$lb
  model$ub <- bounds$ub
  model$Q <- .block_quadratic_sum(
    list(
      list(M = prep$M_alpha_a, idx = idx$idx_a),
      list(M = prep$M_exp_a, idx = idx$idx_r)
    ),
    idx$nvar
  )
  model$obj <- rep(0, idx$nvar)
  model$obj[idx$idx_a] <- prep$q_alpha_a

  tmp <- .add_general_a_constraints(prep, idx)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  result <- .gurobi_solve(model, params = gurobi_params)
  .require_gurobi_optimal(result, status_label)
  .solution_general_a(prep, idx, result)
}

.solve_beta_oracle_a <- function(prep, gurobi_params = list()) {
  .minimax_require("gurobi")
  idx <- .a_gurobi_indices(prep, include_t = FALSE)
  bounds <- .a_gurobi_bounds(prep, idx, include_t = FALSE)

  model <- list()
  model$modelname <- "beta_star_general_a"
  model$modelsense <- "min"
  model$vtype <- base::c(rep("B", prep$p), rep("C", 3L * prep$p),
                   if (prep$use_sets) rep("B", prep$Msets))
  model$lb <- bounds$lb
  model$ub <- bounds$ub
  model$obj <- rep(0, idx$nvar)
  model$Q <- .block_quadratic(prep$M_beta_a, idx$idx_z, idx$nvar)

  tmp <- .add_general_a_constraints(prep, idx)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  result <- .gurobi_solve(model, params = gurobi_params)
  .require_gurobi_optimal(result, "bias oracle")
  .solution_general_a(prep, idx, result)
}

.solve_minimax_a_miqcp <- function(prep, alpha_star, beta_star,
                                   gurobi_params = list()) {
  .minimax_require("gurobi")
  idx <- .a_gurobi_indices(prep, include_t = TRUE)
  bounds <- .a_gurobi_bounds(prep, idx, include_t = TRUE)

  model <- list()
  model$modelname <- "minimax_design_general_a"
  model$modelsense <- "min"
  model$vtype <- base::c(rep("B", prep$p), rep("C", 3L * prep$p),
                   if (prep$use_sets) rep("B", prep$Msets), "C")
  model$lb <- bounds$lb
  model$ub <- bounds$ub
  model$obj <- rep(0, idx$nvar)
  model$obj[idx$idx_t] <- 1

  tmp <- .add_general_a_constraints(prep, idx)
  if (.beta_star_is_zero(prep, beta_star) && any(prep$bias_weights > 0)) {
    row <- .empty_sparse(1L, idx$nvar)
    row[1L, idx$idx_z] <- if (prep$beta_linear_scale > 0) {
      prep$bias_weights / prep$beta_linear_scale
    } else {
      prep$bias_weights
    }
    tmp$A <- rbind(tmp$A, row)
    tmp$rhs <- base::c(tmp$rhs, 0)
    tmp$sense <- base::c(tmp$sense, "=")
  }
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  q_alpha <- rep(0, idx$nvar)
  q_alpha[idx$idx_a] <- prep$q_alpha_a
  q_alpha[idx$idx_t] <- -alpha_star
  rhs_alpha <- -prep$const_alpha_a
  scale_alpha <- .quad_scale_multi(list(prep$M_alpha_a, prep$M_exp_a),
                                   q_alpha, rhs_alpha)

  quadcon <- list(
    list(
      Qc = .block_quadratic_sum(
        list(
          list(M = prep$M_alpha_a, idx = idx$idx_a, scale = scale_alpha),
          list(M = prep$M_exp_a, idx = idx$idx_r, scale = scale_alpha)
        ),
        idx$nvar
      ),
      q = scale_alpha * q_alpha,
      rhs = scale_alpha * rhs_alpha,
      sense = "<",
      name = "alpha_ratio"
    )
  )

  if (!.beta_star_is_zero(prep, beta_star)) {
    q_beta <- rep(0, idx$nvar)
    q_beta[idx$idx_t] <- -beta_star
    rhs_beta <- 0
    scale_beta <- .quad_scale(prep$M_beta_a, q_beta, rhs_beta)
    quadcon[[length(quadcon) + 1L]] <- list(
      Qc = .block_quadratic(prep$M_beta_a, idx$idx_z, idx$nvar, scale_beta),
      q = scale_beta * q_beta,
      rhs = scale_beta * rhs_beta,
      sense = "<",
      name = "beta_ratio"
    )
  }
  model$quadcon <- quadcon

  result <- .gurobi_solve(model, params = gurobi_params)
  .require_gurobi_optimal(result, "minimax program")

  out <- .solution_general_a(prep, idx, result)
  if (.beta_star_is_zero(prep, beta_star)) {
    out$beta <- .clamp_relative_zero(
      out$beta, prep$beta_scale, prep$zero_tol, "beta"
    )
  }
  out$alpha_ratio <- .safe_ratio(out$alpha, alpha_star)
  out$beta_ratio <- .safe_ratio(out$beta, beta_star)
  out$regret <- max(out$alpha_ratio, out$beta_ratio)
  epi <- .validate_regret_epigraph(
    out$regret, as.numeric(result$x[idx$idx_t])
  )
  out$t <- epi$t
  out$solver_t <- epi$solver_t
  out$result$validated_ratio_gap <- epi$validated_gap
  out
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
  model$vtype <- base::c(rep("B", p), rep("C", p), if (prep$use_sets) rep("B", prep$Msets))
  model$lb <- rep(0, nvar)
  model$ub <- base::c(prep$x_max, prep$gamma_upper, if (prep$use_sets) rep(1, prep$Msets))
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

  result <- .gurobi_solve(model, params = gurobi_params)
  .require_gurobi_optimal(result, status_label)
  sol <- result$x
  x <- as.integer(round(sol[idx_x]))
  s <- pmax(prep$gamma_lower * x,
             pmin(prep$gamma_upper * x, as.numeric(sol[idx_s])))
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
  model$modelname <- "beta_star"
  model$modelsense <- "min"
  model$vtype <- base::c(rep("B", p), rep("C", p), if (prep$use_sets) rep("B", prep$Msets))
  model$lb <- rep(0, nvar)
  model$ub <- base::c(prep$x_max, prep$gamma_upper, if (prep$use_sets) rep(1, prep$Msets))
  model$obj <- rep(0, nvar)
  model$obj[idx_s] <- prep$q_beta
  model$Q <- .block_quadratic(prep$M_beta, idx_s, nvar)

  A <- .empty_sparse(0L, nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .add_x_constraints(A, rhs, sense, nvar, idx_x, prep)
  tmp <- .add_s_constraints(tmp$A, tmp$rhs, tmp$sense, nvar, idx_x, idx_s, prep)
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  result <- .gurobi_solve(model, params = gurobi_params)
  .require_gurobi_optimal(result, "bias oracle")
  sol <- result$x
  x <- as.integer(round(sol[idx_x]))
  s <- pmax(prep$gamma_lower * x,
             pmin(prep$gamma_upper * x, as.numeric(sol[idx_s])))
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
  model$vtype <- base::c(rep("B", p), rep("C", p), if (prep$use_sets) rep("B", prep$Msets), "C")
  model$lb <- rep(0, nvar)
  model$ub <- base::c(prep$x_max, prep$gamma_upper, if (prep$use_sets) rep(1, prep$Msets), Inf)
  model$obj <- rep(0, nvar)
  model$obj[idx_t] <- 1

  A <- .empty_sparse(0L, nvar)
  rhs <- numeric(0)
  sense <- character(0)
  tmp <- .add_x_constraints(A, rhs, sense, nvar, idx_x, prep)
  tmp <- .add_s_constraints(tmp$A, tmp$rhs, tmp$sense, nvar, idx_x, idx_s, prep)
  if (.beta_star_is_zero(prep, beta_star) && prep$beta_linear_scale > 0) {
    row <- .empty_sparse(1L, nvar)
    d_normalized <- if (prep$beta_linear_scale > 0) {
      prep$d / prep$beta_linear_scale
    } else {
      prep$d
    }
    row[1L, idx_s] <- d_normalized
    tmp$A <- rbind(tmp$A, row)
    tmp$rhs <- base::c(tmp$rhs, sum(d_normalized))
    tmp$sense <- base::c(tmp$sense, "=")
  }
  model$A <- tmp$A
  model$rhs <- tmp$rhs
  model$sense <- tmp$sense

  q_alpha <- rep(0, nvar)
  q_alpha[idx_s] <- prep$q_alpha
  q_alpha[idx_t] <- -alpha_star
  rhs_alpha <- -prep$const_alpha
  scale_alpha <- .quad_scale(prep$M_alpha, q_alpha, rhs_alpha)

  model$quadcon <- list(
    list(
      Qc = .block_quadratic(prep$M_alpha, idx_s, nvar, scale_alpha),
      q = scale_alpha * q_alpha,
      rhs = scale_alpha * rhs_alpha,
      sense = "<",
      name = "alpha_ratio"
    )
  )
  if (!.beta_star_is_zero(prep, beta_star)) {
    q_beta <- rep(0, nvar)
    q_beta[idx_s] <- prep$q_beta
    q_beta[idx_t] <- -beta_star
    rhs_beta <- -prep$const_beta
    scale_beta <- .quad_scale(prep$M_beta, q_beta, rhs_beta)
    model$quadcon[[length(model$quadcon) + 1L]] <- list(
      Qc = .block_quadratic(prep$M_beta, idx_s, nvar, scale_beta),
      q = scale_beta * q_beta,
      rhs = scale_beta * rhs_beta,
      sense = "<",
      name = "beta_ratio"
    )
  }

  result <- .gurobi_solve(model, params = gurobi_params)
  .require_gurobi_optimal(result, "minimax program")

  sol <- result$x
  x <- as.integer(round(sol[idx_x]))
  s <- pmax(prep$gamma_lower * x,
             pmin(prep$gamma_upper * x, as.numeric(sol[idx_s])))
  gamma <- rep(0, p)
  gamma[x == 1L] <- s[x == 1L]

  alpha <- prep$compute_alpha(s)
  beta <- prep$compute_beta(s)
  if (.beta_star_is_zero(prep, beta_star)) {
    beta <- .clamp_relative_zero(
      beta, prep$beta_scale, prep$zero_tol, "beta"
    )
  }
  alpha_ratio <- .safe_ratio(alpha, alpha_star)
  beta_ratio <- .safe_ratio(beta, beta_star)
  regret <- max(alpha_ratio, beta_ratio)
  epi <- .validate_regret_epigraph(regret, as.numeric(sol[idx_t]))
  result$validated_ratio_gap <- epi$validated_gap

  list(
    x = x,
    s = s,
    gamma = gamma,
    alpha = alpha,
    beta = beta,
    alpha_ratio = alpha_ratio,
    beta_ratio = beta_ratio,
    regret = regret,
    t = epi$t,
    solver_t = epi$solver_t,
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
