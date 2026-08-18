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
                          gamma_lower = NULL,
                          gamma_upper = NULL,
                          experiment_names = NULL,
                          include_variance_design = TRUE,
                          make_plots = TRUE,
                          plot_total_n = FALSE,
                          n_cores = 1L,
                          solver = c("auto", "gurobi", "quadprog"),
                          gurobi_params = list(),
                          max_sets = 200000L,
                          tol_bisect = 1e-7,
                          bisect_iter = 80L,
                          qp_ridge = 1e-10,
                          c = NULL,
                          k = NULL,
                          kappa = NULL) {
  n_grid <- as.numeric(n_grid)
  if (length(n_grid) == 0L || any(!is.finite(n_grid)) || any(n_grid <= 0)) {
    stop("n_grid must contain positive finite values.", call. = FALSE)
  }
  h_values <- as.integer(h_values)
  if (length(h_values) == 0L || any(is.na(h_values)) || any(h_values < 0)) {
    stop("h_values must contain nonnegative integers.", call. = FALSE)
  }

  solver <- match.arg(solver)
  p <- length(omega)
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
  tasks <- expand.grid(
    n_total = n_grid,
    h = h_values,
    method = methods,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

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
      gamma_lower = gamma_lower,
      gamma_upper = gamma_upper,
      solver = solver,
      gurobi_params = gurobi_params,
      max_sets = max_sets,
      qp_ridge = qp_ridge,
      c = c,
      k = k,
      kappa = kappa
    )
    fit <- if (identical(task$method, "minimax")) {
      common$tol_bisect <- tol_bisect
      common$bisect_iter <- bisect_iter
      do.call(solve_minimax_design, common)
    } else {
      do.call(solve_variance_design, common)
    }
    total_alloc <- sum(fit$n_opt)
    arms <- data.frame(
      n_total = task$n_total,
      h = task$h,
      method = task$method,
      solver = fit$solver,
      arm = experiment_names,
      x = as.numeric(fit$x_opt),
      gamma = as.numeric(fit$gamma_opt),
      s = as.numeric(fit$s_opt),
      n_opt = as.numeric(fit$n_opt),
      n_share = if (total_alloc > 0) as.numeric(fit$n_opt) / total_alloc else rep(0, p),
      stringsAsFactors = FALSE
    )
    regret <- data.frame(
      n_total = task$n_total,
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
    data.frame(regret[, c("n_total", "h", "method")], component = "alpha_ratio",
               value = regret$alpha_ratio, stringsAsFactors = FALSE),
    data.frame(regret[, c("n_total", "h", "method")], component = "beta_ratio",
               value = regret$beta_ratio, stringsAsFactors = FALSE),
    data.frame(regret[, c("n_total", "h", "method")], component = "regret",
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
#' Compatibility wrapper around [sweep_n_total()]. This function name is retained
#' for scripts using the Gurobi-only API while the `solver` argument
#' also allows non-Gurobi execution.
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
#' @param k Alias for `bias_weights`, included for script compatibility.
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
                                        solver = c("auto", "gurobi", "quadprog"),
                                        gurobi_params = list(),
                                        max_sets = 200000L,
                                        experiment_names = NULL,
                                        plot_total_n = FALSE,
                                        feasible_sets = NULL,
                                        selection_constraints = NULL,
                                        min_experiments = 0L,
                                        gamma_lower = NULL,
                                        gamma_upper = NULL,
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
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    experiment_names = experiment_names,
    include_variance_design = FALSE,
    make_plots = TRUE,
    plot_total_n = plot_total_n,
    n_cores = ncor,
    solver = solver,
    gurobi_params = gurobi_params,
    max_sets = max_sets,
    c = c,
    k = k
  )
}
