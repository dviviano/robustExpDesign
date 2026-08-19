#' Sweep total budget values and collect design diagnostics
#'
#' @description
#' Runs [solve_minimax_design()] over a grid of total budgets and experiment-count
#' limits. Optionally also computes the variance-optimal benchmark. The returned
#' data frames are suitable for allocation, shrinkage, and regret diagnostics.
#'
#' @inheritParams solve_minimax_design
#' @param n_grid Numeric vector of positive budget values.
#' @param h_values Integer vector of experiment-count upper bounds.
#' @param experiment_names Optional character vector of coordinate labels.
#' @param include_variance_design Logical. If `TRUE`, also runs
#'   [solve_variance_design()].
#' @param make_plots Logical. If `TRUE`, returns `ggplot2` plots when `ggplot2`
#'   is installed.
#' @param plot_total_n Logical. If `TRUE`, allocation plots use `n_opt`; otherwise
#'   they use the allocation share among selected experiments.
#' @param n_cores Integer number of cores for non-Windows parallel execution.
#' @param zero_tol Relative tolerance for homogeneous zero-risk diagnostics.
#' @param psd_tol Relative covariance eigenvalue tolerance.
#' @param psd_action Either `"error"` or `"project"`. Corrections beyond
#'   floating-point roundoff require explicit `"project"`; projections are
#'   reported in fitted objects.
#'
#' @return A list with `data`, `plots`, and `fits`.
#' @export
#' @importFrom rlang .data
sweep_n_total <- function(Sigma_obs,
                          v2,
                          omega,
                          n_grid,
                          h_values = 1L,
                          costs = NULL,
                          bias_weights = NULL,
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
                          experiment_names = NULL,
                          include_variance_design = TRUE,
                          make_plots = TRUE,
                          plot_total_n = FALSE,
                          n_cores = 1L,
                          solver = base::c("auto", "gurobi", "quadprog"),
                          gurobi_params = list(),
                          max_sets = 200000L,
                          tol_bisect = 1e-7,
                          bisect_iter = 80L,
                          qp_ridge = 1e-10,
                          zero_tol = 1e-9,
                          psd_tol = 1e-10,
                          psd_action = base::c("error", "project"),
                          c = NULL,
                          k = NULL) {
  n_grid <- as.numeric(n_grid)
  if (length(n_grid) == 0L || any(!is.finite(n_grid)) || any(n_grid <= 0)) {
    stop("n_grid must contain positive finite values.", call. = FALSE)
  }
  n_min <- as.numeric(n_min)
  if (!(length(n_min) %in% base::c(1L, length(n_grid))) ||
      any(!is.finite(n_min)) || any(n_min < 0)) {
    stop("n_min must be a nonnegative finite scalar or have length equal to n_grid.",
         call. = FALSE)
  }
  n_min_grid <- if (length(n_min) == 1L) rep(n_min, length(n_grid)) else n_min
  h_values_raw <- as.numeric(h_values)
  if (length(h_values_raw) == 0L || any(!is.finite(h_values_raw)) ||
      any(abs(h_values_raw - round(h_values_raw)) > 1e-9) ||
      any(h_values_raw < 0)) {
    stop("h_values must contain nonnegative integers.", call. = FALSE)
  }
  h_values <- as.integer(round(h_values_raw))

  solver <- match.arg(solver)
  psd_action <- match.arg(psd_action)
  p <- length(omega)
  if (any(h_values > p)) {
    stop("h_values cannot exceed length(omega).", call. = FALSE)
  }
  if (is.null(experiment_names)) {
    experiment_names <- names(omega)
  }
  if (is.null(experiment_names) || any(experiment_names == "")) {
    experiment_names <- paste0("arm_", seq_len(p))
  }
  if (length(experiment_names) != p) {
    stop("experiment_names must have length equal to length(omega).", call. = FALSE)
  }

  methods <- "minimax"
  if (isTRUE(include_variance_design)) {
    methods <- base::c(methods, "variance")
  }
  # Keep an explicit index so repeated n_grid values can carry distinct n_min
  # values instead of being collapsed by match().
  tasks <- expand.grid(
    n_index = seq_along(n_grid),
    h = h_values,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  tasks$n_total <- n_grid[tasks$n_index]
  tasks$n_min <- n_min_grid[tasks$n_index]
  tasks$task_id <- seq_len(nrow(tasks))

  run_one <- function(i) {
    task <- tasks[i, ]
    common <- list(
      Sigma_obs = Sigma_obs,
      v2 = v2,
      n_total = task$n_total,
      omega = omega,
      costs = costs,
      bias_weights = bias_weights,
      h = task$h,
      x_max = x_max,
      feasible_sets = feasible_sets,
      selection_constraints = selection_constraints,
      min_experiments = min_experiments,
      n_min = task$n_min,
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper,
      a_exp_lower = a_exp_lower,
      a_exp_upper = a_exp_upper,
      force_gamma_unit_interval = force_gamma_unit_interval,
      solver = solver,
      gurobi_params = gurobi_params,
      max_sets = max_sets,
      qp_ridge = qp_ridge,
      zero_tol = zero_tol,
      psd_tol = psd_tol,
      psd_action = psd_action,
      c = c,
      k = k
    )
    fit <- if (identical(task$method, "minimax")) {
      common$tol_bisect <- tol_bisect
      common$bisect_iter <- bisect_iter
      do.call(solve_minimax_design, common)
    } else {
      do.call(solve_variance_design, common)
    }
    total_alloc <- sum(fit$n_opt)
    total_budget_used <- sum(as.numeric(fit$n_opt) *
                             as.numeric(fit$inputs$costs))
    arms <- data.frame(
      task_id = task$task_id,
      n_index = task$n_index,
      n_total = task$n_total,
      n_min = task$n_min,
      h = task$h,
      method = task$method,
      solver = fit$solver,
      arm = experiment_names,
      x = as.numeric(fit$x_opt),
      gamma = as.numeric(fit$gamma_opt),
      s = as.numeric(fit$s_opt),
      a_exp = if (!is.null(fit$a_exp_opt)) as.numeric(fit$a_exp_opt) else rep(NA_real_, p),
      a_obs = if (!is.null(fit$a_obs_opt)) as.numeric(fit$a_obs_opt) else rep(NA_real_, p),
      n_opt = as.numeric(fit$n_opt),
      n_share = if (total_alloc > 0) {
        as.numeric(fit$n_opt) / total_alloc
      } else {
        rep(0, p)
      },
      budget_used = as.numeric(fit$n_opt) * as.numeric(fit$inputs$costs),
      budget_share = if (total_budget_used > 0) {
        as.numeric(fit$n_opt) * as.numeric(fit$inputs$costs) /
          total_budget_used
      } else {
        rep(0, p)
      },
      stringsAsFactors = FALSE
    )
    regret <- data.frame(
      task_id = task$task_id,
      n_index = task$n_index,
      n_total = task$n_total,
      n_min = task$n_min,
      h = task$h,
      method = task$method,
      solver = fit$solver,
      alpha = fit$alpha_opt,
      beta = fit$beta_opt,
      alpha_star = fit$alpha_star,
      beta_star = fit$beta_star,
      alpha_ratio = fit$alpha_ratio,
      beta_ratio = fit$beta_ratio,
      regret = fit$regret,
      stringsAsFactors = FALSE
    )
    list(arms = arms, regret = regret, fit = fit)
  }

  n_cores <- as.integer(n_cores)
  if (is.na(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }
  if (n_cores > 1L && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(seq_len(nrow(tasks)), run_one, mc.cores = n_cores)
  } else {
    results <- lapply(seq_len(nrow(tasks)), run_one)
  }

  data <- list(
    arms = do.call(rbind, lapply(results, `[[`, "arms")),
    regret = do.call(rbind, lapply(results, `[[`, "regret"))
  )
  rownames(data$arms) <- NULL
  rownames(data$regret) <- NULL

  plots <- list()
  if (isTRUE(make_plots)) {
    plots <- .make_sweep_plots(data, plot_total_n = plot_total_n)
  }
  list(data = data, plots = plots, fits = lapply(results, `[[`, "fit"))
}

.make_sweep_plots <- function(data, plot_total_n = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is not installed; returning data without plots.", call. = FALSE)
    return(list())
  }
  arms <- data$arms
  regret <- data$regret

  y_alloc <- if (isTRUE(plot_total_n)) "n_opt" else "n_share"
  y_alloc_label <- if (isTRUE(plot_total_n)) "Allocated sample / budget" else "Allocation share"

  p_allocation <- ggplot2::ggplot(
    arms,
    ggplot2::aes(
      x = .data$n_total,
      y = .data[[y_alloc]],
      color = .data$arm,
      linetype = .data$method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(ggplot2::vars(h = .data$h), labeller = ggplot2::label_both) +
    ggplot2::labs(x = "n_total", y = y_alloc_label, color = "Arm", linetype = "Design") +
    ggplot2::theme_minimal()

  p_gamma <- ggplot2::ggplot(
    arms,
    ggplot2::aes(
      x = .data$n_total,
      y = .data$gamma,
      color = .data$arm,
      linetype = .data$method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(ggplot2::vars(h = .data$h), labeller = ggplot2::label_both) +
    ggplot2::labs(x = "n_total", y = "gamma", color = "Arm", linetype = "Design") +
    ggplot2::theme_minimal()

  regret_long <- rbind(
    data.frame(regret[, base::c("n_total", "h", "method")], component = "alpha_ratio",
               value = regret$alpha_ratio, stringsAsFactors = FALSE),
    data.frame(regret[, base::c("n_total", "h", "method")], component = "beta_ratio",
               value = regret$beta_ratio, stringsAsFactors = FALSE),
    data.frame(regret[, base::c("n_total", "h", "method")], component = "regret",
               value = regret$regret, stringsAsFactors = FALSE)
  )

  p_regret <- ggplot2::ggplot(
    regret_long,
    ggplot2::aes(
      x = .data$n_total,
      y = .data$value,
      color = .data$component,
      linetype = .data$method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::facet_wrap(ggplot2::vars(h = .data$h), labeller = ggplot2::label_both) +
    ggplot2::labs(x = "n_total", y = "Ratio", color = "Component", linetype = "Design") +
    ggplot2::theme_minimal()

  list(allocation = p_allocation, gamma = p_gamma, regret = p_regret)
}

#' Sweep total sample size and create design plots
#'
#' @description
#' Compatibility wrapper around [sweep_n_total()]. The old function name is kept
#' for scripts that used the original Gurobi-only API, but the `solver` argument
#' now allows non-Gurobi execution.
#'
#' @inheritParams solve_minimax_design
#' @param n_grid Numeric vector of positive budget values.
#' @param h_multi Upper bound used when `BO_experiments = TRUE`.
#' @param h_single Upper bound for the single-count run.
#' @param BO_experiments Logical. If `TRUE`, run both `h_single` and `h_multi`.
#' @param ncor Number of cores.
#' @param experiment_names Optional labels.
#' @param plot_total_n Logical. If `TRUE`, plot allocation levels rather than
#'   allocation shares.
#'
#' @return Output from [sweep_n_total()].
#' @export
sweep_n_tot_and_plot_gurobi <- function(Sigma_obs,
                                        v2,
                                        omega,
                                        n_grid,
                                        costs = NULL,
                                        bias_weights = NULL,
                                        x_max = NULL,
                                        h_multi,
                                        h_single = 1L,
                                        BO_experiments = TRUE,
                                        ncor = 1L,
                                        solver = base::c("auto", "gurobi", "quadprog"),
                                        gurobi_params = list(),
                                        max_sets = 200000L,
                                        experiment_names = NULL,
                                        plot_total_n = FALSE,
                                        feasible_sets = NULL,
                                        selection_constraints = NULL,
                                        min_experiments = 0L,
                                        n_min = 0,
                                        gamma_lower = NULL,
                                        gamma_upper = NULL,
                                        a_exp_lower = NULL,
                                        a_exp_upper = NULL,
                                        force_gamma_unit_interval = TRUE,
                                        tol_bisect = 1e-7,
                                        bisect_iter = 80L,
                                        qp_ridge = 1e-10,
                                        zero_tol = 1e-9,
                                        psd_tol = 1e-10,
                                        psd_action = base::c("error", "project"),
                                        c = NULL,
                                        k = NULL) {
  h_values <- h_single
  if (isTRUE(BO_experiments)) {
    h_values <- unique(base::c(h_values, h_multi))
  }
  sweep_n_total(
    Sigma_obs = Sigma_obs,
    v2 = v2,
    omega = omega,
    n_grid = n_grid,
    h_values = h_values,
    costs = costs,
    bias_weights = bias_weights,
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
    experiment_names = experiment_names,
    include_variance_design = FALSE,
    make_plots = TRUE,
    plot_total_n = plot_total_n,
    n_cores = ncor,
    solver = solver,
    gurobi_params = gurobi_params,
    max_sets = max_sets,
    tol_bisect = tol_bisect,
    bisect_iter = bisect_iter,
    qp_ridge = qp_ridge,
    zero_tol = zero_tol,
    psd_tol = psd_tol,
    psd_action = psd_action,
    c = c,
    k = k
  )
}
