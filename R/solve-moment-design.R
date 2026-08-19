.moment_trace <- function(M) {
  sum(diag(as.matrix(M)))
}

.moment_validate_max_sign_patterns <- function(max_sign_patterns) {
  if (!.is_integerish_scalar(max_sign_patterns) || max_sign_patterns < 2) {
    stop("max_sign_patterns must be an integer of at least 2.", call. = FALSE)
  }
  as.integer(round(max_sign_patterns))
}

.moment_pinv <- function(A, tol = 1e-10) {
  tol <- .validate_tolerance(tol)
  A <- as.matrix(A)
  if (!is.numeric(A) || any(!is.finite(A))) {
    stop("A must be a finite numeric matrix.", call. = FALSE)
  }
  if (nrow(A) == 0L || ncol(A) == 0L) {
    return(matrix(0, ncol(A), nrow(A)))
  }
  s <- svd(A, nu = min(dim(A)), nv = min(dim(A)))
  if (length(s$d) == 0L || max(s$d) == 0) {
    return(matrix(0, ncol(A), nrow(A)))
  }
  cutoff <- tol * max(dim(A)) * max(s$d)
  dinv <- ifelse(s$d > cutoff, 1 / s$d, 0)
  s$v %*% (dinv * t(s$u))
}

.moment_nullspace <- function(A, tol = 1e-10) {
  tol <- .validate_tolerance(tol)
  A <- as.matrix(A)
  if (!is.numeric(A) || any(!is.finite(A))) {
    stop("A must be a finite numeric matrix.", call. = FALSE)
  }
  n <- ncol(A)
  if (n == 0L) return(matrix(0, 0, 0))
  if (nrow(A) == 0L) return(diag(n))

  # Request the full right-singular-vector basis. The default svd() truncates V
  # when nrow(A) < ncol(A), which otherwise drops valid nullspace directions.
  s <- svd(A, nu = 0L, nv = n)
  if (length(s$d) == 0L || max(s$d) == 0) return(diag(n))
  cutoff <- tol * max(dim(A)) * max(s$d)
  rank <- sum(s$d > cutoff)
  if (rank >= n) return(matrix(0, n, 0))
  s$v[, seq.int(rank + 1L, n), drop = FALSE]
}

.moment_gamma_from_W <- function(Lambda, W) {
  A <- t(Lambda) %*% W %*% Lambda
  B <- t(Lambda) %*% W
  out <- tryCatch(solve(A, B), error = function(e) NULL)
  if (is.null(out)) {
    stop("t(Lambda) %*% W %*% Lambda must be nonsingular.", call. = FALSE)
  }
  -out
}

.moment_reconstruct_W <- function(Lambda, Gamma, W_mask = NULL,
                                  tol = 1e-9) {
  tol <- .validate_tolerance(tol)
  Lambda <- as.matrix(Lambda)
  Gamma <- as.matrix(Gamma)
  m <- nrow(Lambda)
  d <- ncol(Lambda)
  if (any(dim(Gamma) != base::c(d, m))) {
    stop("Gamma must be ncol(Lambda) by nrow(Lambda).", call. = FALSE)
  }
  if (is.null(W_mask)) {
    W_mask <- matrix(TRUE, m, m)
  } else {
    if (length(W_mask) != m * m ||
        (!is.null(dim(W_mask)) && any(dim(W_mask) != base::c(m, m)))) {
      stop("W_mask must be an m by m logical matrix.", call. = FALSE)
    }
    W_mask <- matrix(as.logical(W_mask), m, m)
    if (anyNA(W_mask)) stop("W_mask contains NA.", call. = FALSE)
  }
  W <- matrix(0, m, m)
  B <- -Gamma
  for (j in seq_len(m)) {
    rows <- which(W_mask[, j])
    if (length(rows) == 0L) next
    A <- t(Lambda[rows, , drop = FALSE])
    coef <- .moment_pinv(A, tol = tol) %*% B[, j, drop = FALSE]
    W[rows, j] <- as.numeric(coef)
  }
  residual <- t(Lambda) %*% W + Gamma
  err <- if (length(residual) == 0L) 0 else max(abs(residual))
  equation_scale <- max(base::c(abs(t(Lambda) %*% W), abs(Gamma), 0))
  exact <- is.finite(err) &&
    (if (equation_scale == 0) err == 0 else err <= tol * equation_scale)
  list(W = W, max_equation_error = err, equation_scale = equation_scale,
       exact = exact)
}

.moment_bias_value <- function(M, norm = base::c("linf", "l2", "l1"),
                               max_sign_patterns = 4096L) {
  norm <- match.arg(norm)
  max_sign_patterns <- .moment_validate_max_sign_patterns(max_sign_patterns)
  M <- as.matrix(M)
  if (!is.numeric(M) || any(!is.finite(M))) {
    stop("M must be a finite numeric matrix.", call. = FALSE)
  }
  if (ncol(M) == 0L || nrow(M) == 0L) {
    return(0)
  }
  if (identical(norm, "l2")) {
    return(max(svd(M, nu = 0, nv = 0)$d)^2)
  }
  if (identical(norm, "l1")) {
    return(max(colSums(M^2)))
  }
  if (nrow(M) == 1L) {
    return(sum(abs(M))^2)
  }
  r <- ncol(M)
  if (r > floor(log2(max_sign_patterns))) {
    stop("Exact linf-dual bias for multivariate targets requires enumerating sign patterns. Increase max_sign_patterns or use norm = 'l2' or 'l1'.", call. = FALSE)
  }
  signs <- as.matrix(expand.grid(rep(list(base::c(-1, 1)), r)))
  vals <- apply(signs, 1L, function(s) sum((M %*% s)^2))
  max(vals)
}

.moment_eval_from_Gamma <- function(Gamma, Lambda, Sigma, Omega, biased_idx,
                                    bias_weights, norm = base::c("linf", "l2", "l1"),
                                    max_sign_patterns = 4096L) {
  norm <- match.arg(norm)
  Gamma <- as.matrix(Gamma)
  C <- Omega %*% Gamma
  alpha <- .moment_trace(C %*% Sigma %*% t(C))
  if (length(biased_idx) == 0L) {
    beta <- 0
  } else {
    Cb <- C[, biased_idx, drop = FALSE] %*% diag(bias_weights[biased_idx], nrow = length(biased_idx))
    beta <- .moment_bias_value(Cb, norm = norm, max_sign_patterns = max_sign_patterns)
  }
  list(Gamma = Gamma, alpha = as.numeric(alpha), beta = as.numeric(beta))
}

.moment_candidate_masks <- function(m, W_mask = NULL, W_masks = NULL,
                                    moment_mask = NULL, moment_sets = NULL) {
  out <- list()
  add_mask <- function(mask, name) {
    if (length(mask) != m * m ||
        (!is.null(dim(mask)) && any(dim(mask) != base::c(m, m)))) {
      stop("Each W mask must be an m by m logical matrix.", call. = FALSE)
    }
    mask <- matrix(as.logical(mask), m, m)
    if (anyNA(mask)) stop("W_mask contains NA.", call. = FALSE)
    out[[length(out) + 1L]] <<- list(name = name, W_mask = mask)
  }
  if (!is.null(W_masks)) {
    if (!is.list(W_masks) || length(W_masks) == 0L) {
      stop("W_masks must be NULL or a non-empty list of logical matrices.", call. = FALSE)
    }
    for (i in seq_along(W_masks)) add_mask(W_masks[[i]], paste0("W_mask_", i))
  }
  if (!is.null(W_mask)) {
    add_mask(W_mask, "W_mask")
  }
  if (!is.null(moment_sets)) {
    if (!is.list(moment_sets) || length(moment_sets) == 0L) {
      stop("moment_sets must be NULL or a non-empty list of integer vectors.", call. = FALSE)
    }
    for (i in seq_along(moment_sets)) {
      raw <- as.numeric(moment_sets[[i]])
      if (any(!is.finite(raw)) ||
          any(abs(raw - round(raw)) > 1e-9)) {
        stop("moment_sets must contain integer moment indices.", call. = FALSE)
      }
      S <- sort(unique(as.integer(round(raw))))
      if (any(S < 1L | S > m)) {
        stop("moment_sets contains invalid moment indices.", call. = FALSE)
      }
      mask <- matrix(FALSE, m, m)
      if (length(S) > 0L) mask[S, S] <- TRUE
      add_mask(mask, paste0("moment_set_", i))
    }
  }
  if (!is.null(moment_mask)) {
    mm <- as.logical(moment_mask)
    if (length(mm) != m || anyNA(mm)) {
      stop("moment_mask must be a logical vector with one entry per moment.", call. = FALSE)
    }
    mask <- matrix(FALSE, m, m)
    idx <- which(mm)
    if (length(idx) > 0L) mask[idx, idx] <- TRUE
    add_mask(mask, "moment_mask")
  }
  if (length(out) == 0L) {
    out[[1L]] <- list(name = "all_moments", W_mask = matrix(TRUE, m, m))
  }

  key <- vapply(out, function(z) paste(as.integer(z$W_mask), collapse = ""), character(1))
  out[!duplicated(key)]
}

.cvxr_wmask_constraints <- function(Gamma, Lambda, W_mask, tol = 1e-10) {
  d <- ncol(Lambda)
  m <- nrow(Lambda)
  constraints <- list()
  for (j in seq_len(m)) {
    rows <- which(W_mask[, j])
    if (length(rows) == 0L) {
      constraints <- base::c(constraints, list(Gamma[, j] == matrix(0, d, 1)))
    } else {
      A <- t(Lambda[rows, , drop = FALSE])
      N <- .moment_nullspace(t(A), tol = tol)
      if (ncol(N) > 0L) {
        constraints <- base::c(constraints, list(t(N) %*% Gamma[, j] == matrix(0, ncol(N), 1)))
      }
    }
  }
  constraints
}

.moment_zero_bias_feasible <- function(Lambda, Omega, biased_idx,
                                       bias_weights, W_mask, tol = 1e-9) {
  tol <- .validate_tolerance(tol)
  Lambda <- as.matrix(Lambda)
  Omega <- as.matrix(Omega)
  d <- ncol(Lambda)
  m <- nrow(Lambda)
  nvar <- d * m

  # vec(Gamma %*% Lambda) = (t(Lambda) kron I_d) vec(Gamma).
  Aeq <- kronecker(t(Lambda), diag(d))
  beq <- as.vector(-diag(d))

  add_equations <- function(M, rhs = NULL) {
    if (nrow(M) == 0L) return(invisible(NULL))
    Aeq <<- rbind(Aeq, M)
    beq <<- base::c(beq, if (is.null(rhs)) rep(0, nrow(M)) else rhs)
    invisible(NULL)
  }

  for (j in seq_len(m)) {
    rows <- which(W_mask[, j])
    if (length(rows) == 0L) {
      block <- matrix(0, nrow = d, ncol = nvar)
      block[, ((j - 1L) * d + 1L):(j * d)] <- diag(d)
      add_equations(block)
    } else {
      A <- t(Lambda[rows, , drop = FALSE])
      N <- .moment_nullspace(t(A), tol = tol)
      if (ncol(N) > 0L) {
        block <- matrix(0, nrow = ncol(N), ncol = nvar)
        block[, ((j - 1L) * d + 1L):(j * d)] <- t(N)
        add_equations(block)
      }
    }
  }

  relevant <- intersect(biased_idx, which(bias_weights > 0))
  if (length(relevant) > 0L) {
    for (j in relevant) {
      block <- matrix(0, nrow = nrow(Omega), ncol = nvar)
      block[, ((j - 1L) * d + 1L):(j * d)] <- Omega
      add_equations(block)
    }
  }

  solution <- .moment_pinv(Aeq, tol = tol) %*% beq
  fitted <- as.numeric(Aeq %*% solution)
  residual <- if (length(beq) == 0L) 0 else max(abs(fitted - beq))
  equation_scale <- max(base::c(abs(fitted), abs(beq), 0))
  feasible <- is.finite(residual) &&
    (if (equation_scale == 0) residual == 0 else residual <= tol * equation_scale)
  list(feasible = feasible, residual = residual,
       equation_scale = equation_scale)
}

.cvxr_bias_radius <- function(C_bias, norm = base::c("linf", "l2", "l1"),
                              max_sign_patterns = 4096L) {
  norm <- match.arg(norm)
  max_sign_patterns <- .moment_validate_max_sign_patterns(max_sign_patterns)
  dims <- dim(C_bias)
  q <- dims[1]
  r <- dims[2]
  if (r == 0L) {
    return(list(radius = 0, constraints = list()))
  }
  if (identical(norm, "l2")) {
    return(list(radius = CVXR::norm(C_bias, "2"), constraints = list()))
  }
  if (identical(norm, "l1")) {
    rho <- CVXR::Variable(1)
    constraints <- list(rho >= 0)
    for (j in seq_len(r)) {
      constraints <- base::c(constraints, list(CVXR::p_norm(C_bias[, j], 2) <= rho))
    }
    return(list(radius = rho, constraints = constraints))
  }
  if (q == 1L) {
    return(list(radius = CVXR::p_norm(CVXR::vec(C_bias), 1), constraints = list()))
  }
  if (r > floor(log2(max_sign_patterns))) {
    stop("Exact linf-dual bias for multivariate targets requires too many sign patterns. Increase max_sign_patterns or use norm = 'l2' or 'l1'.", call. = FALSE)
  }
  rho <- CVXR::Variable(1)
  constraints <- list(rho >= 0)
  signs <- as.matrix(expand.grid(rep(list(base::c(-1, 1)), r)))
  for (i in seq_len(nrow(signs))) {
    s <- matrix(as.numeric(signs[i, ]), ncol = 1)
    constraints <- base::c(constraints, list(CVXR::p_norm(C_bias %*% s, 2) <= rho))
  }
  list(radius = rho, constraints = constraints)
}

.cvxr_solve <- function(problem, cvxr_solver = NULL, cvxr_args = list()) {
  args <- base::c(list(problem), cvxr_args)

  if (!is.null(cvxr_solver)) {
    args$solver <- cvxr_solver
  }

  objective_value <- do.call(CVXR::psolve, args)
  problem_status <- as.character(CVXR::status(problem))

  ok <- problem_status %in% base::c(
    "optimal",
    "optimal_inaccurate",
    "OPTIMAL",
    "OPTIMAL_INACCURATE"
  )

  if (!ok) {
    stop(
      sprintf("CVXR solver status: %s", problem_status),
      call. = FALSE
    )
  }

  list(
    status = problem_status,
    value = objective_value,
    getValue = function(x) CVXR::value(x)
  )
}

.cvxr_optimize_gamma <- function(Lambda, Sigma, Omega, biased_idx, bias_weights,
                                 W_mask, mode = base::c("alpha", "beta", "minimax"),
                                 alpha_star = NULL, beta_star = NULL,
                                 norm = base::c("linf", "l2", "l1"),
                                 max_sign_patterns = 4096L,
                                 cvxr_solver = NULL,
                                 cvxr_args = list(),
                                 Sigma_root = NULL,
                                 beta_zero = FALSE,
                                 beta_reference = NULL,
                                 zero_tol = 1e-9,
                                 psd_tol = 1e-10,
                                 tol = 1e-8) {
  mode <- match.arg(mode)
  norm <- match.arg(norm)
  .minimax_require("CVXR")

  m <- nrow(Lambda)
  d <- ncol(Lambda)
  Gamma <- CVXR::Variable(d, m)
  if (is.null(Sigma_root)) {
    Sigma_root <- .psd_square_root(Sigma, tol = psd_tol,
                                   action = "error", name = "Sigma")$root
  }
  C <- Omega %*% Gamma
  alpha_expr <- CVXR::sum_squares(C %*% Sigma_root)

  constraints <- list(Gamma %*% Lambda == -diag(d))
  constraints <- base::c(constraints, .cvxr_wmask_constraints(Gamma, Lambda, W_mask, tol = tol))

  zero_expr <- CVXR::sum_squares(Gamma) * 0
  if (length(biased_idx) == 0L) {
    radius <- zero_expr
    radius_constraints <- list()
  } else {
    C_bias <- C[, biased_idx, drop = FALSE] %*%
      diag(bias_weights[biased_idx], nrow = length(biased_idx))
    br <- .cvxr_bias_radius(C_bias, norm = norm, max_sign_patterns = max_sign_patterns)
    radius <- br$radius
    radius_constraints <- br$constraints
  }
  constraints <- base::c(constraints, radius_constraints)

  if (identical(mode, "alpha")) {
    objective <- CVXR::Minimize(alpha_expr)
  } else if (identical(mode, "beta")) {
    objective <- CVXR::Minimize(radius)
  } else {
    if (is.null(alpha_star) || is.null(beta_star)) {
      stop("alpha_star and beta_star are required for minimax mode.", call. = FALSE)
    }
    tvar <- CVXR::Variable(1)
    constraints <- base::c(constraints, list(tvar >= 1, alpha_expr <= alpha_star * tvar))
    if (isTRUE(beta_zero)) {
      constraints <- base::c(constraints, list(radius <= 0))
    } else {
      constraints <- base::c(constraints, list(radius <= sqrt(beta_star * tvar)))
    }
    objective <- CVXR::Minimize(tvar)
  }

  problem <- CVXR::Problem(objective, constraints)
  result <- .cvxr_solve(problem, cvxr_solver = cvxr_solver, cvxr_args = cvxr_args)
  Gamma_val <- as.matrix(result$getValue(Gamma))
  eval <- .moment_eval_from_Gamma(
    Gamma_val, Lambda = Lambda, Sigma = Sigma, Omega = Omega,
    biased_idx = biased_idx, bias_weights = bias_weights, norm = norm,
    max_sign_patterns = max_sign_patterns
  )
  eval$status <- result$status
  eval$objective_value <- result$value
  if (identical(mode, "minimax")) {
    eval$t <- as.numeric(result$getValue(tvar))
    if (isTRUE(beta_zero)) {
      C_val <- Omega %*% Gamma_val
      if (length(biased_idx) == 0L) {
        bias_residual <- 0
        bias_scale <- 0
      } else {
        relevant <- intersect(biased_idx, which(bias_weights > 0))
        if (length(relevant) == 0L) {
          bias_residual <- 0
          bias_scale <- 0
        } else {
          C_bias_val <- C_val[, relevant, drop = FALSE] %*%
            diag(bias_weights[relevant], nrow = length(relevant))
          bias_residual <- max(abs(C_bias_val))
          bias_scale <- max(abs(Omega)) * max(abs(Gamma_val)) *
            max(1, ncol(Omega)) * max(bias_weights[relevant])
        }
      }
      bias_tol <- max(tol, sqrt(.Machine$double.eps)) *
        max(bias_scale, .Machine$double.xmin)
      if (bias_residual > bias_tol) {
        stop(sprintf(
          "CVXR zero-bias solution violates the zero-bias equations by %.6g (tolerance %.6g).",
          bias_residual, bias_tol
        ), call. = FALSE)
      }
      eval$beta <- 0
      eval$zero_bias_residual <- bias_residual
      eval$zero_bias_tolerance <- bias_tol
    }
    eval$alpha_ratio <- .safe_ratio(eval$alpha, alpha_star, tol = zero_tol)
    eval$beta_ratio <- .safe_ratio(eval$beta, beta_star, tol = zero_tol)
    eval$regret <- max(eval$alpha_ratio, eval$beta_ratio)
  }
  eval
}

#' Evaluate or optimize a general-loading moment design
#'
#' @description
#' Implements the general-loading, fixed-covariance GMM/moment-selection regret formulas. With a
#' user-supplied `W` and `optimize_W = FALSE`, the function evaluates
#' `alpha_Omega(W, Sigma)`, `beta_l,Omega(W)`, and the corresponding GMM linear
#' expansion. With `optimize_W = TRUE`, it solves the convex equivalent problem
#' over the linear GMM influence matrix `Gamma`, subject to the supplied `W_mask`,
#' `W_masks`, `moment_mask`, or `moment_sets`.
#'
#' @param Lambda Numeric `m` by `d` Jacobian of moments with respect to parameters.
#' @param Sigma Numeric `m` by `m` covariance matrix for the moment vector.
#' @param Omega Numeric `q` by `d` target-sensitivity matrix. For a scalar target,
#'   you may pass `omega` instead.
#' @param omega Optional length-`d` vector used as `Omega = matrix(omega, 1)`.
#' @param biased_moments Integer or logical index for moments that may be biased
#'   (`I^c` in the paper). Defaults to all moments.
#' @param bias_weights Optional nonnegative length-`m` vector. These rescale the
#'   bias radius by multiplying the biased columns of the matrix product of
#'   `Omega` and `Gamma`.
#' @param norm Bias ambiguity norm: `"linf"`, `"l2"`, or `"l1"`.
#' @param W Optional user-supplied weighting matrix. Required when
#'   `optimize_W = FALSE`.
#' @param W_mask Optional logical `m` by `m` mask of allowed entries in `W`.
#' @param W_masks Optional list of `W_mask` matrices. When supplied and
#'   `optimize_W = TRUE`, the function solves each mask and selects the best
#'   minimax-regret solution.
#' @param moment_mask Optional logical length-`m` vector of moments that may enter.
#' @param moment_sets Optional list of integer moment sets to compare.
#' @param optimize_W Logical. If `FALSE`, evaluate the supplied `W`; if `TRUE`,
#'   optimize over the convex equivalent linear GMM estimator under the mask(s).
#' @param cvxr_solver Optional CVXR solver name, e.g. `"ECOS"`, `"CLARABEL"`, or
#'   `"SCS"` depending on what is installed.
#' @param cvxr_args Optional named list passed to `CVXR::psolve()`.
#' @param max_sign_patterns Maximum sign patterns for exact multivariate
#'   `norm = "linf"` bias calculations.
#' @param reconstruct_W Logical. If `TRUE`, reconstructs one weighting matrix that
#'   implements the optimized `Gamma` and respects the selected mask when possible.
#' @param zero_tol Relative tolerance for homogeneous zero-risk diagnostics.
#' @param psd_tol Relative covariance eigenvalue tolerance.
#' @param psd_action Either `"error"` or `"project"`. Corrections beyond
#'   floating-point roundoff require explicit `"project"`; projections are
#'   reported in the fitted object.
#' @param tol Numerical tolerance for dimensionless linear-algebra checks.
#'
#' @return A list of class `moment_design`.
#' @export
solve_moment_design <- function(Lambda,
                                Sigma,
                                Omega = NULL,
                                omega = NULL,
                                biased_moments = NULL,
                                bias_weights = NULL,
                                norm = base::c("linf", "l2", "l1"),
                                W = NULL,
                                W_mask = NULL,
                                W_masks = NULL,
                                moment_mask = NULL,
                                moment_sets = NULL,
                                optimize_W = is.null(W),
                                cvxr_solver = NULL,
                                cvxr_args = list(),
                                max_sign_patterns = 4096L,
                                reconstruct_W = TRUE,
                                zero_tol = 1e-9,
                                psd_tol = 1e-10,
                                psd_action = base::c("error", "project"),
                                tol = 1e-8) {
  norm <- match.arg(norm)
  max_sign_patterns <- .moment_validate_max_sign_patterns(max_sign_patterns)
  Lambda <- as.matrix(Lambda)
  Sigma <- as.matrix(Sigma)
  if (!is.numeric(Lambda) || any(!is.finite(Lambda))) {
    stop("Lambda must be a finite numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(Sigma) || any(!is.finite(Sigma))) {
    stop("Sigma must be a finite numeric matrix.", call. = FALSE)
  }
  m <- nrow(Lambda)
  d <- ncol(Lambda)
  if (m < 1L || d < 1L) {
    stop("Lambda must have at least one row and one column.", call. = FALSE)
  }
  if (any(dim(Sigma) != base::c(m, m))) {
    stop("Sigma must be an m by m matrix with m = nrow(Lambda).", call. = FALSE)
  }
  zero_tol <- .validate_tolerance(zero_tol, "zero_tol")
  psd_action <- match.arg(psd_action)
  checked_sigma <- .validate_psd_matrix(
    Sigma, name = "Sigma", tol = psd_tol, action = psd_action
  )
  Sigma <- checked_sigma$matrix
  Sigma_root <- .psd_square_root(
    Sigma, tol = psd_tol, action = "error", name = "Sigma"
  )$root
  if (is.null(Omega)) {
    if (is.null(omega)) {
      stop("Supply either Omega or omega.", call. = FALSE)
    }
    omega <- as.numeric(omega)
    if (length(omega) != d) {
      stop("omega must have length ncol(Lambda).", call. = FALSE)
    }
    Omega <- matrix(omega, nrow = 1L)
  } else {
    Omega <- as.matrix(Omega)
  }
  if (!is.numeric(Omega) || ncol(Omega) != d || any(!is.finite(Omega))) {
    stop("Omega must be a finite numeric matrix with ncol(Omega) = ncol(Lambda).", call. = FALSE)
  }
  Omega_original <- Omega
  target_scale <- max(abs(Omega_original))
  if (target_scale == 0) target_scale <- 1
  Omega_internal <- Omega_original / target_scale

  if (is.null(biased_moments)) {
    biased_idx <- seq_len(m)
  } else if (is.logical(biased_moments)) {
    if (length(biased_moments) != m || anyNA(biased_moments)) {
      stop("logical biased_moments must have length nrow(Lambda).", call. = FALSE)
    }
    biased_idx <- which(biased_moments)
  } else {
    biased_raw <- as.numeric(biased_moments)
    if (any(!is.finite(biased_raw)) ||
        any(abs(biased_raw - round(biased_raw)) > 1e-9)) {
      stop("biased_moments must contain integer indices.", call. = FALSE)
    }
    biased_idx <- sort(unique(as.integer(round(biased_raw))))
    if (any(biased_idx < 1L | biased_idx > m)) {
      stop("biased_moments contains invalid indices.", call. = FALSE)
    }
  }
  bias_weights <- .numeric_vector(bias_weights, "bias_weights", m, default = rep(1, m))
  if (any(bias_weights < 0)) {
    stop("bias_weights must be nonnegative.", call. = FALSE)
  }

  if (!isTRUE(optimize_W)) {
    if (is.null(W)) {
      stop("W must be supplied when optimize_W = FALSE.", call. = FALSE)
    }
    W <- as.matrix(W)
    if (!is.numeric(W) || any(dim(W) != base::c(m, m)) || any(!is.finite(W))) {
      stop("W must be a finite numeric m by m matrix.", call. = FALSE)
    }
    if (!is.null(W_mask)) {
      if (length(W_mask) != m * m ||
          (!is.null(dim(W_mask)) && any(dim(W_mask) != base::c(m, m)))) {
        stop("W_mask must be an m by m logical matrix.", call. = FALSE)
      }
      W_mask <- matrix(as.logical(W_mask), m, m)
      if (anyNA(W_mask)) stop("W_mask contains NA.", call. = FALSE)
      mask_scale <- max(base::c(abs(W), 0))
      outside_error <- if (any(!W_mask)) max(abs(W[!W_mask])) else 0
      if (!.relative_zero(outside_error, mask_scale, tol = tol)) {
        stop("W has nonzero entries outside W_mask.", call. = FALSE)
      }
    }
    Gamma <- .moment_gamma_from_W(Lambda, W)
    eval <- .moment_eval_from_Gamma(
      Gamma, Lambda = Lambda, Sigma = Sigma, Omega = Omega_internal,
      biased_idx = biased_idx, bias_weights = bias_weights, norm = norm,
      max_sign_patterns = max_sign_patterns
    )
    eval$alpha <- eval$alpha * target_scale^2
    eval$beta <- eval$beta * target_scale^2
    out <- base::c(eval, list(
      W = W,
      W_mask = W_mask,
      biased_moments = biased_idx,
      norm = norm,
      optimize_W = FALSE,
      psd_diagnostics = checked_sigma$diagnostics,
      inputs = list(Lambda = Lambda, Sigma = Sigma, Omega = Omega_original,
                    bias_weights = bias_weights, psd_tol = psd_tol,
                    psd_action = psd_action)
    ))
    class(out) <- base::c("moment_design", "list")
    return(out)
  }

  .minimax_require("CVXR")
  Omega <- Omega_internal
  candidates <- .moment_candidate_masks(
    m = m, W_mask = W_mask, W_masks = W_masks,
    moment_mask = moment_mask, moment_sets = moment_sets
  )
  zero_bias_checks <- lapply(candidates, function(candidate) {
    .moment_zero_bias_feasible(
      Lambda = Lambda, Omega = Omega, biased_idx = biased_idx,
      bias_weights = bias_weights, W_mask = candidate$W_mask, tol = tol
    )
  })
  zero_bias_feasible <- vapply(
    zero_bias_checks, `[[`, logical(1), "feasible"
  )

  alpha_solutions <- vector("list", length(candidates))
  beta_solutions <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    mask <- candidates[[i]]$W_mask
    alpha_solutions[[i]] <- tryCatch(
      .cvxr_optimize_gamma(
        Lambda, Sigma, Omega, biased_idx, bias_weights, mask,
        mode = "alpha", norm = norm, max_sign_patterns = max_sign_patterns,
        cvxr_solver = cvxr_solver, cvxr_args = cvxr_args,
        Sigma_root = Sigma_root, zero_tol = zero_tol,
        psd_tol = psd_tol, tol = tol
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "moment_error")
    )
    beta_solutions[[i]] <- tryCatch(
      .cvxr_optimize_gamma(
        Lambda, Sigma, Omega, biased_idx, bias_weights, mask,
        mode = "beta", norm = norm, max_sign_patterns = max_sign_patterns,
        cvxr_solver = cvxr_solver, cvxr_args = cvxr_args,
        Sigma_root = Sigma_root, zero_tol = zero_tol,
        psd_tol = psd_tol, tol = tol
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "moment_error")
    )
  }

  alpha_vals <- vapply(alpha_solutions, function(z) if (inherits(z, "moment_error")) Inf else z$alpha, numeric(1))
  beta_vals <- vapply(beta_solutions, function(z) if (inherits(z, "moment_error")) Inf else z$beta, numeric(1))
  successful_zero <- zero_bias_feasible & is.finite(beta_vals)
  if (any(successful_zero)) {
    for (i in which(successful_zero)) {
      beta_solutions[[i]]$beta <- 0
      beta_vals[i] <- 0
    }
  }
  beta_zero <- any(successful_zero)
  alpha_star <- min(alpha_vals)
  beta_star <- if (beta_zero) 0 else min(beta_vals)
  beta_reference <- max(base::c(beta_vals[is.finite(beta_vals)], 0))
  if (!is.finite(alpha_star) || alpha_star <= 0) {
    stop("No feasible moment-design variance oracle was found.", call. = FALSE)
  }
  if (!is.finite(beta_star)) {
    stop("No feasible moment-design bias oracle was found.", call. = FALSE)
  }

  minimax_solutions <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    mask <- candidates[[i]]$W_mask
    minimax_solutions[[i]] <- tryCatch(
      .cvxr_optimize_gamma(
        Lambda, Sigma, Omega, biased_idx, bias_weights, mask,
        mode = "minimax", alpha_star = alpha_star, beta_star = beta_star,
        norm = norm, max_sign_patterns = max_sign_patterns,
        cvxr_solver = cvxr_solver, cvxr_args = cvxr_args,
        Sigma_root = Sigma_root, beta_zero = beta_zero,
        beta_reference = beta_reference, zero_tol = zero_tol,
        psd_tol = psd_tol, tol = tol
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "moment_error")
    )
  }
  regret_vals <- vapply(minimax_solutions, function(z) if (inherits(z, "moment_error")) Inf else z$regret, numeric(1))
  best_i <- which.min(regret_vals)
  if (!is.finite(regret_vals[best_i])) {
    stop("No feasible moment-design minimax solution was found.", call. = FALSE)
  }
  best <- minimax_solutions[[best_i]]
  best_mask <- candidates[[best_i]]$W_mask

  W_rec <- NULL
  W_rec_error <- NA_real_
  W_rec_exact <- NA
  if (isTRUE(reconstruct_W)) {
    rec <- .moment_reconstruct_W(Lambda, best$Gamma, best_mask, tol = tol)
    W_rec <- rec$W
    W_rec_error <- rec$max_equation_error
    W_rec_exact <- rec$exact
  }

  rescale_solution <- function(z) {
    if (inherits(z, "moment_error") || !is.list(z)) return(z)
    for (nm in base::c("alpha", "beta")) {
      if (!is.null(z[[nm]])) z[[nm]] <- z[[nm]] * target_scale^2
    }
    z
  }
  best <- rescale_solution(best)
  alpha_solutions_out <- lapply(alpha_solutions, rescale_solution)
  beta_solutions_out <- lapply(beta_solutions, rescale_solution)
  minimax_solutions_out <- lapply(minimax_solutions, rescale_solution)
  alpha_star_natural <- alpha_star * target_scale^2
  beta_star_natural <- beta_star * target_scale^2

  out <- list(
    selected_candidate = candidates[[best_i]]$name,
    candidate_index = best_i,
    W_mask = best_mask,
    Gamma_opt = best$Gamma,
    W_reconstructed = W_rec,
    W_reconstruction_error = W_rec_error,
    W_reconstruction_exact = W_rec_exact,
    alpha_opt = best$alpha,
    beta_opt = best$beta,
    alpha_star = alpha_star_natural,
    beta_star = beta_star_natural,
    beta_zero = beta_zero,
    alpha_ratio = best$alpha_ratio,
    beta_ratio = best$beta_ratio,
    regret = best$regret,
    t_opt = best$t,
    norm = norm,
    biased_moments = biased_idx,
    optimize_W = TRUE,
    alpha_oracle = alpha_solutions_out[[which.min(alpha_vals)]],
    beta_oracle = beta_solutions_out[[which.min(beta_vals)]],
    candidates = candidates,
    minimax_solutions = minimax_solutions_out,
    candidate_diagnostics = data.frame(
      candidate = vapply(candidates, `[[`, character(1), "name"),
      alpha_oracle = alpha_vals * target_scale^2,
      beta_oracle = beta_vals * target_scale^2,
      zero_bias_feasible = zero_bias_feasible,
      zero_bias_residual = vapply(
        zero_bias_checks, `[[`, numeric(1), "residual"
      ),
      minimax_regret = regret_vals,
      stringsAsFactors = FALSE
    ),
    psd_diagnostics = checked_sigma$diagnostics,
    inputs = list(Lambda = Lambda, Sigma = Sigma, Omega = Omega_original,
                  bias_weights = bias_weights, cvxr_solver = cvxr_solver,
                  zero_tol = zero_tol, psd_tol = psd_tol,
                  psd_action = psd_action,
                  internal_target_scale = target_scale)
  )
  class(out) <- base::c("moment_design", "list")
  out
}

#' @export
print.moment_design <- function(x, ...) {
  cat("Moment-selection design\n")
  if (isTRUE(x$optimize_W)) {
    cat("Selected candidate: ", x$selected_candidate, "\n", sep = "")
    cat("Regret: ", format(x$regret, digits = 6), "\n", sep = "")
    cat("Variance ratio: ", format(x$alpha_ratio, digits = 6), "\n", sep = "")
    cat("Bias ratio: ", format(x$beta_ratio, digits = 6), "\n", sep = "")
  } else {
    cat("User-supplied W evaluation\n")
    cat("alpha: ", format(x$alpha, digits = 6), "\n", sep = "")
    cat("beta: ", format(x$beta, digits = 6), "\n", sep = "")
  }
  invisible(x)
}
