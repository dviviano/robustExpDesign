.moment_trace <- function(M) {
  sum(diag(as.matrix(M)))
}

.moment_pinv <- function(A, tol = 1e-10) {
  A <- as.matrix(A)
  s <- svd(A)
  if (length(s$d) == 0L) {
    return(matrix(0, ncol(A), nrow(A)))
  }
  cutoff <- tol * max(dim(A)) * max(s$d)
  dinv <- ifelse(s$d > cutoff, 1 / s$d, 0)
  s$v %*% (dinv * t(s$u))
}

.moment_nullspace <- function(A, tol = 1e-10) {
  A <- as.matrix(A)
  n <- ncol(A)
  if (n == 0L) {
    return(matrix(0, 0, 0))
  }
  s <- svd(A)
  if (length(s$d) == 0L) {
    return(diag(n))
  }
  cutoff <- tol * max(dim(A)) * max(s$d)
  rank <- sum(s$d > cutoff)
  if (rank >= n) {
    return(matrix(0, n, 0))
  }
  s$v[, (rank + 1L):n, drop = FALSE]
}

.moment_gamma_from_W <- function(Lambda, W) {
  A <- t(Lambda) %*% W %*% Lambda
  B <- t(Lambda) %*% W
  -solve(A, B)
}

.moment_reconstruct_W <- function(Lambda, Gamma, W_mask = NULL, tol = 1e-9) {
  Lambda <- as.matrix(Lambda)
  Gamma <- as.matrix(Gamma)
  m <- nrow(Lambda)
  d <- ncol(Lambda)
  if (is.null(W_mask)) {
    W_mask <- matrix(TRUE, m, m)
  } else {
    W_mask <- matrix(as.logical(W_mask), m, m)
  }
  W <- matrix(0, m, m)
  B <- -Gamma
  for (j in seq_len(m)) {
    rows <- which(W_mask[, j])
    if (length(rows) == 0L) {
      next
    }
    A <- t(Lambda[rows, , drop = FALSE])
    coef <- .moment_pinv(A, tol = tol) %*% B[, j, drop = FALSE]
    W[rows, j] <- as.numeric(coef)
  }
  err <- max(abs(t(Lambda) %*% W + Gamma))
  list(W = W, max_equation_error = err, exact = is.finite(err) && err <= sqrt(tol))
}

.moment_resolve_bias_weights <- function(bias_weights, k = NULL, kappa = NULL) {
  out <- bias_weights
  aliases <- list(k = k, kappa = kappa)
  for (nm in names(aliases)) {
    alias <- aliases[[nm]]
    if (!is.null(alias)) {
      if (!is.null(out) && !isTRUE(all.equal(as.numeric(out), as.numeric(alias)))) {
        stop("Use only one of bias_weights, k, or kappa, or supply identical values.",
             call. = FALSE)
      }
      out <- alias
    }
  }
  out
}

.moment_active_from_mask <- function(W_mask) {
  W_mask <- as.matrix(W_mask)
  which(rowSums(W_mask != 0) > 0 | colSums(W_mask != 0) > 0)
}

.moment_bias_value <- function(M, norm = c("linf", "l2", "l1"),
                               max_sign_patterns = 4096L) {
  norm <- match.arg(norm)
  M <- as.matrix(M)
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
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), r)))
  vals <- apply(signs, 1L, function(s) sum((M %*% s)^2))
  max(vals)
}

.moment_eval_from_Gamma <- function(Gamma, Lambda, Sigma, Omega, biased_idx,
                                    bias_weights, norm = c("linf", "l2", "l1"),
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
      S <- sort(unique(as.integer(moment_sets[[i]])))
      if (anyNA(S) || any(S < 1L | S > m)) {
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
      constraints <- c(constraints, list(Gamma[, j] == matrix(0, d, 1)))
    } else {
      A <- t(Lambda[rows, , drop = FALSE])
      N <- .moment_nullspace(t(A), tol = tol)
      if (ncol(N) > 0L) {
        constraints <- c(constraints, list(t(N) %*% Gamma[, j] == matrix(0, ncol(N), 1)))
      }
    }
  }
  constraints
}

.cvxr_bias_radius <- function(C_bias, norm = c("linf", "l2", "l1"),
                              max_sign_patterns = 4096L) {
  norm <- match.arg(norm)
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
      constraints <- c(constraints, list(CVXR::p_norm(C_bias[, j], 2) <= rho))
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
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), r)))
  for (i in seq_len(nrow(signs))) {
    s <- matrix(as.numeric(signs[i, ]), ncol = 1)
    constraints <- c(constraints, list(CVXR::p_norm(C_bias %*% s, 2) <= rho))
  }
  list(radius = rho, constraints = constraints)
}

.cvxr_solve <- function(problem, cvxr_solver = NULL, cvxr_args = list()) {
  args <- c(list(object = problem), cvxr_args)
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

.cvxr_optimize_gamma <- function(Lambda, Sigma, Omega, biased_idx, bias_weights,
                                 W_mask, mode = c("alpha", "beta", "minimax"),
                                 alpha_star = NULL, beta_star = NULL,
                                 norm = c("linf", "l2", "l1"),
                                 max_sign_patterns = 4096L,
                                 cvxr_solver = NULL,
                                 cvxr_args = list(),
                                 tol = 1e-8) {
  mode <- match.arg(mode)
  norm <- match.arg(norm)
  .minimax_require("CVXR")

  m <- nrow(Lambda)
  d <- ncol(Lambda)
  Gamma <- CVXR::Variable(d, m)
  R <- chol(Sigma)
  C <- Omega %*% Gamma
  alpha_expr <- CVXR::sum_squares(C %*% t(R))

  constraints <- list(Gamma %*% Lambda == -diag(d))
  constraints <- c(constraints, .cvxr_wmask_constraints(Gamma, Lambda, W_mask, tol = tol))

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
  constraints <- c(constraints, radius_constraints)

  if (identical(mode, "alpha")) {
    objective <- CVXR::Minimize(alpha_expr)
  } else if (identical(mode, "beta")) {
    objective <- CVXR::Minimize(radius)
  } else {
    if (is.null(alpha_star) || is.null(beta_star)) {
      stop("alpha_star and beta_star are required for minimax mode.", call. = FALSE)
    }
    tvar <- CVXR::Variable(1)
    constraints <- c(constraints, list(tvar >= 1, alpha_expr <= alpha_star * tvar))
    if (beta_star <= tol) {
      constraints <- c(constraints, list(radius <= 0))
    } else {
      constraints <- c(constraints, list(radius <= CVXR::sqrt(beta_star * tvar)))
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
    eval$alpha_ratio <- .safe_ratio(eval$alpha, alpha_star)
    eval$beta_ratio <- .safe_ratio(eval$beta, beta_star)
    eval$regret <- max(eval$alpha_ratio, eval$beta_ratio)
  }
  eval
}

#' Evaluate or optimize a Section 4 moment-selection design
#'
#' @description
#' Implements the Section 4 GMM/moment-selection regret formulas. With a
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
#'   bias radius by multiplying the biased columns of `Omega %*% Gamma`.
#' @param k,kappa Optional aliases for `bias_weights`. `kappa` matches the
#'   weighted-bias notation `|b_j| <= kappa[j] * B`.
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
#' @param cvxr_args Optional named list passed to `CVXR::solve()`.
#' @param max_sign_patterns Maximum sign patterns for exact multivariate
#'   `norm = "linf"` bias calculations.
#' @param reconstruct_W Logical. If `TRUE`, reconstructs one weighting matrix that
#'   implements the optimized `Gamma` and respects the selected mask when possible.
#' @param tol Numerical tolerance.
#'
#' @return A list of class `moment_design`.
#' @export
solve_moment_design <- function(Lambda,
                                Sigma,
                                Omega = NULL,
                                omega = NULL,
                                biased_moments = NULL,
                                bias_weights = NULL,
                                k = NULL,
                                kappa = NULL,
                                norm = c("linf", "l2", "l1"),
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
                                tol = 1e-8) {
  norm <- match.arg(norm)
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
  if (any(dim(Sigma) != c(m, m))) {
    stop("Sigma must be an m by m matrix with m = nrow(Lambda).", call. = FALSE)
  }
  Sigma <- 0.5 * (Sigma + t(Sigma))
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

  if (is.null(biased_moments)) {
    biased_idx <- seq_len(m)
  } else if (is.logical(biased_moments)) {
    if (length(biased_moments) != m || anyNA(biased_moments)) {
      stop("logical biased_moments must have length nrow(Lambda).", call. = FALSE)
    }
    biased_idx <- which(biased_moments)
  } else {
    biased_idx <- sort(unique(as.integer(biased_moments)))
    if (anyNA(biased_idx) || any(biased_idx < 1L | biased_idx > m)) {
      stop("biased_moments contains invalid indices.", call. = FALSE)
    }
  }
  bias_weights <- .moment_resolve_bias_weights(bias_weights, k = k, kappa = kappa)
  bias_weights <- .numeric_vector(bias_weights, "bias_weights", m, default = rep(1, m))
  if (any(bias_weights < 0)) {
    stop("bias_weights must be nonnegative.", call. = FALSE)
  }

  if (!isTRUE(optimize_W)) {
    if (is.null(W)) {
      stop("W must be supplied when optimize_W = FALSE.", call. = FALSE)
    }
    W <- as.matrix(W)
    if (!is.numeric(W) || any(dim(W) != c(m, m)) || any(!is.finite(W))) {
      stop("W must be a finite numeric m by m matrix.", call. = FALSE)
    }
    if (!is.null(W_mask)) {
      W_mask <- matrix(as.logical(W_mask), m, m)
      if (any(abs(W[!W_mask]) > tol)) {
        stop("W has nonzero entries outside W_mask.", call. = FALSE)
      }
    }
    Gamma <- .moment_gamma_from_W(Lambda, W)
    eval <- .moment_eval_from_Gamma(
      Gamma, Lambda = Lambda, Sigma = Sigma, Omega = Omega,
      biased_idx = biased_idx, bias_weights = bias_weights, norm = norm,
      max_sign_patterns = max_sign_patterns
    )
    out <- c(eval, list(
      W = W,
      W_mask = W_mask,
      biased_moments = biased_idx,
      active_moments = if (is.null(W_mask)) seq_len(m) else .moment_active_from_mask(W_mask),
      norm = norm,
      optimize_W = FALSE,
      inputs = list(Lambda = Lambda, Sigma = Sigma, Omega = Omega,
                    bias_weights = bias_weights)
    ))
    class(out) <- c("moment_design", "list")
    return(out)
  }

  .minimax_require("CVXR")
  candidates <- .moment_candidate_masks(
    m = m, W_mask = W_mask, W_masks = W_masks,
    moment_mask = moment_mask, moment_sets = moment_sets
  )

  alpha_solutions <- vector("list", length(candidates))
  beta_solutions <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    mask <- candidates[[i]]$W_mask
    alpha_solutions[[i]] <- tryCatch(
      .cvxr_optimize_gamma(
        Lambda, Sigma, Omega, biased_idx, bias_weights, mask,
        mode = "alpha", norm = norm, max_sign_patterns = max_sign_patterns,
        cvxr_solver = cvxr_solver, cvxr_args = cvxr_args, tol = tol
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "moment_error")
    )
    beta_solutions[[i]] <- tryCatch(
      .cvxr_optimize_gamma(
        Lambda, Sigma, Omega, biased_idx, bias_weights, mask,
        mode = "beta", norm = norm, max_sign_patterns = max_sign_patterns,
        cvxr_solver = cvxr_solver, cvxr_args = cvxr_args, tol = tol
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "moment_error")
    )
  }

  alpha_vals <- vapply(alpha_solutions, function(z) if (inherits(z, "moment_error")) Inf else z$alpha, numeric(1))
  beta_vals <- vapply(beta_solutions, function(z) if (inherits(z, "moment_error")) Inf else z$beta, numeric(1))
  alpha_star <- min(alpha_vals)
  beta_star <- min(beta_vals)
  if (!is.finite(alpha_star) || alpha_star <= 0) {
    stop("No feasible Section 4 variance oracle was found.", call. = FALSE)
  }
  if (!is.finite(beta_star)) {
    stop("No feasible Section 4 bias oracle was found.", call. = FALSE)
  }

  minimax_solutions <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    mask <- candidates[[i]]$W_mask
    minimax_solutions[[i]] <- tryCatch(
      .cvxr_optimize_gamma(
        Lambda, Sigma, Omega, biased_idx, bias_weights, mask,
        mode = "minimax", alpha_star = alpha_star, beta_star = beta_star,
        norm = norm, max_sign_patterns = max_sign_patterns,
        cvxr_solver = cvxr_solver, cvxr_args = cvxr_args, tol = tol
      ),
      error = function(e) structure(list(error = conditionMessage(e)), class = "moment_error")
    )
  }
  regret_vals <- vapply(minimax_solutions, function(z) if (inherits(z, "moment_error")) Inf else z$regret, numeric(1))
  best_i <- which.min(regret_vals)
  if (!is.finite(regret_vals[best_i])) {
    stop("No feasible Section 4 minimax solution was found.", call. = FALSE)
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

  out <- list(
    selected_candidate = candidates[[best_i]]$name,
    candidate_index = best_i,
    W_mask = best_mask,
    active_moments = .moment_active_from_mask(best_mask),
    Gamma_opt = best$Gamma,
    W_reconstructed = W_rec,
    W_reconstruction_error = W_rec_error,
    W_reconstruction_exact = W_rec_exact,
    alpha_opt = best$alpha,
    beta_opt = best$beta,
    alpha_star = alpha_star,
    beta_star = beta_star,
    alpha_ratio = best$alpha_ratio,
    beta_ratio = best$beta_ratio,
    regret = best$regret,
    t_opt = best$t,
    norm = norm,
    biased_moments = biased_idx,
    optimize_W = TRUE,
    alpha_oracle = alpha_solutions[[which.min(alpha_vals)]],
    beta_oracle = beta_solutions[[which.min(beta_vals)]],
    candidates = candidates,
    candidate_diagnostics = data.frame(
      candidate = vapply(candidates, `[[`, character(1), "name"),
      alpha_oracle = alpha_vals,
      beta_oracle = beta_vals,
      minimax_regret = regret_vals,
      stringsAsFactors = FALSE
    ),
    inputs = list(Lambda = Lambda, Sigma = Sigma, Omega = Omega,
                  bias_weights = bias_weights, cvxr_solver = cvxr_solver)
  )
  class(out) <- c("moment_design", "list")
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
