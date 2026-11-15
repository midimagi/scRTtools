# Track-length velocity models.
# summarize_track_lengths(), fit_track_velocity_linear() and
# fit_track_length_models() merge the inline code of notebook
# 2.3.model_track_length_velocity (src and DS trees, identical apart from
# paths) with fit_one_state() / fit_track_length_models() of
# jupyter_downsample_frag_quantile/retreat_figure6_helpers.R (the newer,
# refactored copy of the same code). fit_track_velocity_ode() comes from
# notebook 2.3.model_track_length_velocity_ode; predict_rep_percent() from the
# inverse_predict() helpers of both 2.3 notebooks. Plots are not included.

#' Per-cell median track length
#'
#' Summarise a track table per cell (and optional extra grouping columns):
#' `median_length = median(width_bp, na.rm = TRUE)`, `log_length =
#' log(median_length)` and `rep_percent` = first `rep_percent` of the cell
#' times 100.
#'
#' @param track_lengths A `data.frame` with columns `cell_id`, `width_bp` and
#'   `rep_percent` (fraction), e.g. from [calc_rep_track_lengths()].
#' @param by Extra grouping columns (e.g. `"state"`).
#'
#' @return A `data.frame` with the grouping columns, `median_length`,
#'   `log_length` and `rep_percent` (in percent), ordered by the grouping
#'   columns (C locale).
#'
#' @export
summarize_track_lengths <- function(track_lengths, by = character()) {
  keys <- c("cell_id", by)
  need <- c(keys, "width_bp", "rep_percent")
  if (!is.data.frame(track_lengths) || !all(need %in% names(track_lengths))) {
    stop("'track_lengths' must be a data.frame with columns: ",
         paste(need, collapse = ", "))
  }
  if (nrow(track_lengths) == 0L) {
    out <- track_lengths[0, keys, drop = FALSE]
    out$median_length <- numeric(0)
    out$log_length <- numeric(0)
    out$rep_percent <- numeric(0)
    return(out)
  }
  key_df <- track_lengths[keys]
  ord <- do.call(order, c(unname(as.list(key_df)), list(method = "radix")))
  tl <- track_lengths[ord, , drop = FALSE]
  grp <- do.call(paste, c(unname(as.list(tl[keys])), sep = "\r"))
  first <- !duplicated(grp)
  idx <- split(seq_len(nrow(tl)), factor(grp, levels = unique(grp)))

  out <- tl[first, keys, drop = FALSE]
  out$median_length <- vapply(idx, function(i) stats::median(tl$width_bp[i], na.rm = TRUE),
                              numeric(1), USE.NAMES = FALSE)
  out$log_length <- log(out$median_length)
  out$rep_percent <- tl$rep_percent[first] * 100
  rownames(out) <- NULL
  out
}

#' Robust quadratic model of track length versus replication percentage
#'
#' Fit `log_length ~ rep_percent + I(rep_percent^2)` with
#' `robustbase::lmrob()` (plus `lm()` with a linear term and, optionally, a
#' median regression with `quantreg::rq()`), predict on an even grid with
#' confidence and prediction intervals, derive the length velocity
#' `exp(fit) * (b1 + 2 b2 x)` and validate on a random hold-out split.
#'
#' Cells with non-finite `rep_percent` or `median_length`, or
#' `median_length <= 0`, are dropped and (with `sort_by_rep = TRUE`) the data
#' are sorted by `rep_percent`. Intervals use `qt((1 + level) / 2, df.residual)` with the
#' lmrob standard errors and residual scale.
#'
#' `lmrob()` draws random subsamples for its initial estimate; with
#' `fit_seed = NULL` the current random-number state is used, otherwise the fits run under `withr::with_seed(fit_seed)`.
#'
#' @param data Output of [summarize_track_lengths()] (`rep_percent`,
#'   `median_length`, `log_length`).
#' @param n_grid Number of prediction points (e.g. `100` or `200`).
#' @param level Confidence level of the intervals.
#' @param quantile Logical; also fit the median regression (needs quantreg).
#' @param holdout_prop Fraction of cells used for training in the hold-out
#'   validation (`floor(holdout_prop * n)`).
#' @param seed Seed of the hold-out split.
#' @param fit_seed `NULL` or a seed for the full-data fits (see Details).
#' @param positive_velocity_only Logical; drop prediction rows with
#'   `length_velocity <= 0`.
#' @param min_cells,min_distinct Minimum number of cells and of distinct
#'   `rep_percent` values (e.g. `20` and `8`).
#' @param sort_by_rep Logical; sort the cells by `rep_percent` before
#'   fitting (default); `FALSE` keeps the input order. The order changes the hold-out split and the row order of
#'   `data`, and moves the `lmrob()` fits within its convergence tolerance.
#'
#' @return An object of class `scrt_track_velocity` (a list) with `model =
#'   "linear"`, the fits `fit_linear`, `fit_lmrob`, `fit_rq`, `data` (with
#'   `is_outlier` = zero robustness weight), `prediction`, `linear_residuals`,
#'   `robust_residuals`, `validation`, `validation_metrics` (`rmse`,
#'   `pseudo_r2`, `correlation`), `model_comparison` (long table of the
#'   lmrob / rq predicted lengths), `minimal_predicted_length`
#'   (`exp(b0 - b1^2 / (4 b2))`) and `rep_range`.
#'
#' @seealso [fit_track_length_models()], [fit_track_velocity_ode()],
#'   [predict_rep_percent()]
#' @export
fit_track_velocity_linear <- function(data,
                                      n_grid = 100L,
                                      level = 0.95,
                                      quantile = TRUE,
                                      holdout_prop = 0.5,
                                      seed = 2026L,
                                      fit_seed = NULL,
                                      positive_velocity_only = TRUE,
                                      min_cells = 20L,
                                      min_distinct = 8L,
                                      sort_by_rep = TRUE) {
  .fit_track_lm_one(data, group_name = NULL, group_col = NULL, seed = seed,
                    fit_seed = fit_seed, n_grid = n_grid, level = level,
                    quantile = quantile, holdout_prop = holdout_prop,
                    positive_velocity_only = positive_velocity_only,
                    min_cells = min_cells, min_distinct = min_distinct,
                    sort_by_rep = sort_by_rep)
}

#' Fit track-length velocity models per group
#'
#' Run [fit_track_velocity_linear()] separately for each level of
#' `group_col` and stack the result tables. Group `i` uses
#' hold-out seed `seed + i`.
#'
#' @param data Output of [summarize_track_lengths()] with a grouping column.
#' @param group_col Grouping column.
#' @param levels Group levels to fit (in this order); default: sorted unique
#'   values.
#' @param seed Base seed of the hold-out splits.
#' @param n_grid Number of prediction points per group (e.g. `200`).
#' @param ... Further arguments of [fit_track_velocity_linear()] (`level`,
#'   `quantile`, `holdout_prop`, `fit_seed`, `positive_velocity_only`,
#'   `min_cells`, `min_distinct`, `sort_by_rep`).
#'
#' @return A list with `state_models` (per-group results), `has_quantreg`, and
#'   the stacked tables `model_data`, `model_prediction`, `linear_residuals`,
#'   `robust_residuals`, `validation`, `validation_metrics` and
#'   `model_comparison`, each with `group_col` as a factor.
#'
#' @export
fit_track_length_models <- function(data,
                                    group_col = "state",
                                    levels = NULL,
                                    seed = 2026L,
                                    n_grid = 200L,
                                    ...) {
  if (!is.data.frame(data) || !group_col %in% names(data)) {
    stop("'data' must be a data.frame with column '", group_col, "'.")
  }
  if (is.null(levels)) {
    levels <- sort(unique(as.character(data[[group_col]])))
  }
  dots <- list(...)
  has_quantreg <- if (is.null(dots$quantile)) {
    requireNamespace("quantreg", quietly = TRUE)
  } else {
    isTRUE(dots$quantile)
  }
  dots$quantile <- has_quantreg

  models <- stats::setNames(lapply(seq_along(levels), function(i) {
    lev <- levels[[i]]
    do.call(.fit_track_lm_one, c(
      list(data = data[!is.na(data[[group_col]]) & data[[group_col]] == lev, , drop = FALSE],
           group_name = lev, group_col = group_col, seed = seed + i, n_grid = n_grid),
      dots
    ))
  }), levels)

  tables <- c("data", "prediction", "linear_residuals", "robust_residuals",
              "validation", "validation_metrics", "model_comparison")
  result <- lapply(tables, function(tab) {
    obj <- do.call(rbind, lapply(models, `[[`, tab))
    rownames(obj) <- NULL
    obj[[group_col]] <- factor(obj[[group_col]], levels = levels)
    obj
  })
  names(result) <- sub("^(data|prediction)$", "model_\\1", tables)
  c(list(state_models = models, has_quantreg = has_quantreg), result)
}

#' One linear (quadratic robust) track-length fit
#' @keywords internal
#' @noRd
.fit_track_lm_one <- function(data, group_name, group_col, seed, fit_seed = NULL,
                              n_grid = 100L, level = 0.95, quantile = TRUE,
                              holdout_prop = 0.5, positive_velocity_only = TRUE,
                              min_cells = 20L, min_distinct = 8L,
                              sort_by_rep = TRUE) {
  rlang::check_installed("robustbase", reason = "to fit robust track-length models.")
  if (isTRUE(quantile)) {
    rlang::check_installed("quantreg", reason = "for median regression (or set quantile = FALSE).")
  }
  .validate_track_model_data(data)
  if (!is.numeric(level) || length(level) != 1L || !(level > 0 && level < 1)) {
    stop("'level' must be a single number in (0, 1).")
  }
  label <- if (is.null(group_name)) "Data" else paste0("State '", group_name, "'")

  data <- as.data.frame(data)
  keep <- is.finite(data$rep_percent) & is.finite(data$median_length) & data$median_length > 0
  data <- data[keep, , drop = FALSE]
  if (isTRUE(sort_by_rep)) data <- data[order(data$rep_percent), , drop = FALSE]
  rownames(data) <- NULL
  if (nrow(data) < min_cells || length(unique(data$rep_percent)) < min_distinct) {
    stop(label, " has too few cells/replication values for modeling.")
  }

  form_q <- log_length ~ rep_percent + I(rep_percent^2)
  run_fits <- function() {
    list(
      fit_linear = stats::lm(log_length ~ rep_percent, data = data),
      fit_lmrob = robustbase::lmrob(form_q, data = data),
      fit_rq = if (isTRUE(quantile)) quantreg::rq(form_q, tau = 0.5, data = data) else NULL
    )
  }
  fits <- if (is.null(fit_seed)) run_fits() else withr::with_seed(fit_seed, run_fits())
  fit_linear <- fits$fit_linear
  fit_lmrob <- fits$fit_lmrob
  fit_rq <- fits$fit_rq

  prediction <- data.frame(
    rep_percent = seq(min(data$rep_percent, na.rm = TRUE),
                      max(data$rep_percent, na.rm = TRUE),
                      length.out = n_grid)
  )
  lmrob_pred <- stats::predict(fit_lmrob, newdata = prediction, se.fit = TRUE)
  rq_pred <- if (isTRUE(quantile)) {
    as.numeric(stats::predict(fit_rq, newdata = prediction))
  } else {
    rep(NA_real_, nrow(prediction))
  }
  critical_value <- stats::qt((1 + level) / 2, df = fit_lmrob$df.residual)
  residual_scale <- if (!is.null(lmrob_pred$residual.scale)) {
    lmrob_pred$residual.scale
  } else {
    fit_lmrob$scale
  }
  half_ci <- critical_value * lmrob_pred$se.fit
  half_pi <- critical_value * sqrt(lmrob_pred$se.fit^2 + residual_scale^2)

  if (!is.null(group_name)) prediction[[group_col]] <- group_name
  lmrob_log <- as.numeric(lmrob_pred$fit)
  prediction$lmrob_log <- lmrob_log
  prediction$lmrob_log_lower_ci <- lmrob_log - half_ci
  prediction$lmrob_log_upper_ci <- lmrob_log + half_ci
  prediction$lmrob_log_lower_pi <- lmrob_log - half_pi
  prediction$lmrob_log_upper_pi <- lmrob_log + half_pi
  prediction$lmrob_length <- exp(lmrob_log)
  prediction$lmrob_length_lower_ci <- exp(prediction$lmrob_log_lower_ci)
  prediction$lmrob_length_upper_ci <- exp(prediction$lmrob_log_upper_ci)
  prediction$lmrob_length_lower_pi <- exp(prediction$lmrob_log_lower_pi)
  prediction$lmrob_length_upper_pi <- exp(prediction$lmrob_log_upper_pi)
  prediction$rq_length <- exp(rq_pred)

  coefs <- stats::coef(fit_lmrob)
  b0 <- coefs["(Intercept)"]
  b1 <- coefs["rep_percent"]
  b2 <- coefs["I(rep_percent^2)"]
  prediction$length_velocity <- prediction$lmrob_length * (b1 + 2 * b2 * prediction$rep_percent)
  if (isTRUE(positive_velocity_only)) {
    prediction <- prediction[prediction$length_velocity > 0, , drop = FALSE]
    rownames(prediction) <- NULL
  }

  data$is_outlier <- fit_lmrob$rweights == 0

  add_group <- function(df) {
    if (is.null(group_name)) return(df)
    cbind(stats::setNames(data.frame(group_name, stringsAsFactors = FALSE), group_col), df)
  }
  linear_residuals <- add_group(data.frame(fitted = stats::fitted(fit_linear),
                                           resid = stats::residuals(fit_linear)))
  robust_residuals <- add_group(data.frame(fitted = stats::fitted(fit_lmrob),
                                           resid = stats::residuals(fit_lmrob)))

  holdout <- withr::with_seed(seed, {
    train_indices <- sample(seq_len(nrow(data)), size = floor(holdout_prop * nrow(data)))
    train_data <- data[train_indices, , drop = FALSE]
    test_data <- data[-train_indices, , drop = FALSE]
    fit_train <- robustbase::lmrob(form_q, data = train_data)
    list(actual = test_data$log_length,
         predicted = as.numeric(stats::predict(fit_train, newdata = test_data)))
  })
  validation <- add_group(data.frame(actual = holdout$actual, predicted = holdout$predicted))
  validation_metrics <- add_group(.holdout_metrics(holdout$actual, holdout$predicted))

  model_comparison <- .long_model_comparison(prediction, group_col, group_name)

  structure(list(
    model = "linear",
    fit_linear = fit_linear,
    fit_lmrob = fit_lmrob,
    fit_rq = fit_rq,
    data = data,
    prediction = prediction,
    linear_residuals = linear_residuals,
    robust_residuals = robust_residuals,
    validation = validation,
    validation_metrics = validation_metrics,
    model_comparison = model_comparison,
    minimal_predicted_length = unname(exp(b0 - b1^2 / 4 / b2)),
    rep_range = range(data$rep_percent)
  ), class = "scrt_track_velocity")
}

#' @keywords internal
#' @noRd
.validate_track_model_data <- function(data) {
  need <- c("rep_percent", "median_length", "log_length")
  if (!is.data.frame(data) || !all(need %in% names(data))) {
    stop("'data' must be a data.frame with columns: ", paste(need, collapse = ", "),
         " (see summarize_track_lengths()).")
  }
  invisible(TRUE)
}

#' RMSE, pseudo-R2 and correlation of a hold-out prediction
#' @keywords internal
#' @noRd
.holdout_metrics <- function(actual, predicted) {
  sse <- sum((actual - predicted)^2)
  sst <- sum((actual - mean(actual))^2)
  data.frame(
    rmse = sqrt(mean((actual - predicted)^2)),
    pseudo_r2 = 1 - sse / sst,
    correlation = stats::cor(actual, predicted)
  )
}

#' Long table of lmrob / rq predicted lengths (tidyr::pivot_longer order)
#' @keywords internal
#' @noRd
.long_model_comparison <- function(prediction, group_col, group_name,
                                   cols = c("lmrob_length", "rq_length")) {
  n <- nrow(prediction)
  out <- data.frame(rep_percent = rep(prediction$rep_percent, each = length(cols)))
  if (!is.null(group_name)) {
    out <- cbind(stats::setNames(data.frame(rep(group_name, nrow(out)),
                                            stringsAsFactors = FALSE), group_col), out)
  }
  out$model <- rep(cols, times = n)
  out$predicted_length <- as.vector(t(as.matrix(prediction[cols])))
  out <- out[is.finite(out$predicted_length), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Robust ODE model of track length versus replication percentage
#'
#' Fit the closed-form solution of the ODE
#' `d log L / dx = alpha + beta * exp(gamma * (x - rep_center))`,
#' `log L(x) = log_L_ref + alpha * d + beta / gamma * expm1(gamma * d)` with
#' `d = x - rep_center`, by robust nonlinear regression
#' (`robustbase::nlrob(method = "M")`, `beta` and `gamma` on the log scale)
#' and optionally by median regression (`quantreg::nlrq()`, started at the
#' robust estimates). Starting values come from a robust quadratic fit
#' (`gamma = 0.03`, `beta = max(2 * b2 / gamma, 1e-4)`, `alpha = b1 - beta`).
#'
#' Predictions on an even grid get delta-method confidence intervals and
#' prediction intervals with `qnorm((1 + level) / 2)`; the length velocity is
#' `L(x) * (alpha + beta * exp(gamma * d))`. A random hold-out split refits
#' the model (with `rep_center` of the full data) and reports RMSE,
#' pseudo-R2 and correlation on the test cells.
#'
#' @inheritParams fit_track_velocity_linear
#' @param rep_center Centre of the parameterisation; default: median of
#'   `data$rep_percent`.
#' @param start `NULL` (computed) or a named list of starting values
#'   (`log_L_ref`, `alpha`, `log_beta`, `log_gamma`).
#' @param maxiter `maxiter` of `stats::nls.control()` for `nlrob()`.
#'
#' @return An object of class `scrt_track_velocity` (a list) with `model =
#'   "ode"`, `rep_center`, `start`, `fit_ode`, `fit_nlrq`, `parameters`
#'   (`rep_center`, `log_L_ref`, `alpha`, `beta`, `gamma`), `data` (with
#'   `is_outlier`), `residuals`, `prediction` (all grid rows),
#'   `velocity` (rows with positive velocity), `validation`,
#'   `validation_metrics`, `rse`, `pseudo_r2`, `turning_rep_percent`,
#'   `minimal_predicted_length` and `rep_range`.
#'
#' @export
fit_track_velocity_ode <- function(data,
                                   rep_center = NULL,
                                   start = NULL,
                                   n_grid = 100L,
                                   level = 0.95,
                                   quantile = TRUE,
                                   holdout_prop = 0.5,
                                   seed = 2026L,
                                   fit_seed = NULL,
                                   maxiter = 500L) {
  rlang::check_installed("robustbase", reason = "to fit the robust ODE model.")
  if (isTRUE(quantile)) {
    rlang::check_installed("quantreg", reason = "for median regression (or set quantile = FALSE).")
  }
  .validate_track_model_data(data)
  data <- as.data.frame(data)
  if (is.null(rep_center)) rep_center <- stats::median(data$rep_percent, na.rm = TRUE)
  if (!is.numeric(rep_center) || length(rep_center) != 1L || !is.finite(rep_center)) {
    stop("'rep_center' must be a single finite number.")
  }

  ode_log_length <- .ode_log_length_fun(rep_center)
  # nlrob() / nlrq() evaluate the model outside the formula environment, so the
  # notebook's ode_log_length() call is expanded into the formula (same
  # operations in the same order; rep_center inlined as a constant).
  ode_formula <- .ode_formula(rep_center)
  fit_ode_model <- function(d, st = NULL) {
    if (is.null(st)) st <- .ode_start(d, rep_center)
    robustbase::nlrob(
      ode_formula,
      data = d,
      start = st,
      method = "M",
      doCov = TRUE,
      model = TRUE,
      control = stats::nls.control(maxiter = maxiter),
      trace = FALSE
    )
  }

  run_fits <- function() {
    st <- if (is.null(start)) .ode_start(data, rep_center) else start
    fit_ode <- fit_ode_model(data, st)
    fit_nlrq <- if (isTRUE(quantile)) {
      quantreg::nlrq(
        ode_formula,
        tau = 0.5,
        data = data,
        start = as.list(stats::coef(fit_ode)),
        trace = FALSE
      )
    } else {
      NULL
    }
    list(start = st, fit_ode = fit_ode, fit_nlrq = fit_nlrq)
  }
  fits <- if (is.null(fit_seed)) run_fits() else withr::with_seed(fit_seed, run_fits())
  fit_ode <- fits$fit_ode
  fit_nlrq <- fits$fit_nlrq
  ode_coef <- stats::coef(fit_ode)

  parameters <- data.frame(
    parameter = c("rep_center", "log_L_ref", "alpha", "beta", "gamma"),
    estimate = c(rep_center, ode_coef["log_L_ref"], ode_coef["alpha"],
                 exp(ode_coef["log_beta"]), exp(ode_coef["log_gamma"])),
    row.names = NULL
  )
  residuals <- data.frame(fitted = stats::fitted(fit_ode), resid = stats::residuals(fit_ode))

  min_x <- min(data$rep_percent, na.rm = TRUE)
  max_x <- max(data$rep_percent, na.rm = TRUE)
  prediction <- data.frame(rep_percent = seq(min_x, max_x, length.out = n_grid))

  ode_pred_log <- .ode_predict_from_coef(prediction$rep_percent, ode_coef, ode_log_length)
  ode_grad <- .ode_gradient(prediction$rep_percent, ode_coef, rep_center)
  ode_pred_se <- sqrt(pmax(rowSums((ode_grad %*% stats::vcov(fit_ode)) * ode_grad), 0))
  crit <- stats::qnorm((1 + level) / 2)
  half_ci <- crit * ode_pred_se
  half_pi <- crit * sqrt(ode_pred_se^2 + stats::sigma(fit_ode)^2)
  nlrq_pred_log <- if (isTRUE(quantile)) {
    .ode_predict_from_coef(prediction$rep_percent, stats::coef(fit_nlrq), ode_log_length)
  } else {
    rep(NA_real_, nrow(prediction))
  }

  prediction$ode_log <- ode_pred_log
  prediction$ode_log_lower_ci <- ode_pred_log - half_ci
  prediction$ode_log_upper_ci <- ode_pred_log + half_ci
  prediction$ode_log_lower_pi <- ode_pred_log - half_pi
  prediction$ode_log_upper_pi <- ode_pred_log + half_pi
  prediction$ode_length <- exp(ode_pred_log)
  prediction$ode_length_lower_ci <- exp(prediction$ode_log_lower_ci)
  prediction$ode_length_upper_ci <- exp(prediction$ode_log_upper_ci)
  prediction$ode_length_lower_pi <- exp(prediction$ode_log_lower_pi)
  prediction$ode_length_upper_pi <- exp(prediction$ode_log_upper_pi)
  prediction$nlrq_log <- nlrq_pred_log
  prediction$nlrq_length <- exp(nlrq_pred_log)

  alpha <- ode_coef["alpha"]
  beta <- exp(ode_coef["log_beta"])
  gamma <- exp(ode_coef["log_gamma"])
  prediction$diff_term <- alpha + beta * exp(gamma * (prediction$rep_percent - rep_center))
  prediction$length_velocity <- prediction$ode_length * prediction$diff_term
  velocity <- prediction[prediction$length_velocity > 0, , drop = FALSE]
  rownames(velocity) <- NULL

  holdout <- withr::with_seed(seed, {
    train_indices <- sample(seq_len(nrow(data)), size = floor(holdout_prop * nrow(data)))
    train_data <- data[train_indices, , drop = FALSE]
    test_data <- data[-train_indices, , drop = FALSE]
    fit_train <- fit_ode_model(train_data)
    list(actual = test_data$log_length,
         predicted = .ode_predict_from_coef(test_data$rep_percent, stats::coef(fit_train),
                                            ode_log_length))
  })
  validation <- data.frame(actual = holdout$actual, predicted = holdout$predicted)

  data$is_outlier <- stats::weights(fit_ode, type = "robustness") == 0
  rse <- stats::sigma(fit_ode)
  pseudo_r2 <- 1 - sum(stats::residuals(fit_ode)^2) /
    sum((data$log_length - mean(data$log_length))^2)

  turning <- .ode_turning_point(alpha, beta, gamma, rep_center, min_x)
  turning <- min(max(turning, min_x), max_x)
  minimal_length <- exp(ode_log_length(turning, ode_coef["log_L_ref"], alpha,
                                       log(beta), log(gamma)))

  structure(list(
    model = "ode",
    rep_center = rep_center,
    start = fits$start,
    fit_ode = fit_ode,
    fit_nlrq = fit_nlrq,
    parameters = parameters,
    data = data,
    residuals = residuals,
    prediction = prediction,
    velocity = velocity,
    validation = validation,
    validation_metrics = .holdout_metrics(holdout$actual, holdout$predicted),
    rse = unname(rse),
    pseudo_r2 = unname(pseudo_r2),
    turning_rep_percent = unname(turning),
    minimal_predicted_length = unname(minimal_length),
    rep_range = c(min_x, max_x)
  ), class = "scrt_track_velocity")
}

#' ODE log-length function with a fixed rep_center
#' @keywords internal
#' @noRd
.ode_log_length_fun <- function(rep_center) {
  force(rep_center)
  function(rep_percent, log_L_ref, alpha, log_beta, log_gamma) {
    beta <- exp(log_beta)
    gamma <- exp(log_gamma)
    rep_delta <- rep_percent - rep_center
    log_L_ref +
      alpha * rep_delta +
      beta / gamma * expm1(gamma * rep_delta)
  }
}

#' Model formula with the ODE solution written out
#' @keywords internal
#' @noRd
.ode_formula <- function(rep_center) {
  stats::as.formula(bquote(
    log_length ~ log_L_ref + alpha * (rep_percent - .(rep_center)) +
      exp(log_beta) / exp(log_gamma) * expm1(exp(log_gamma) * (rep_percent - .(rep_center)))
  ), env = baseenv())
}

#' Starting values from a robust quadratic fit
#' @keywords internal
#' @noRd
.ode_start <- function(data, rep_center) {
  start_data <- data
  start_data$rep_delta <- start_data$rep_percent - rep_center
  fit_quadratic_start <- robustbase::lmrob(
    log_length ~ rep_delta + I(rep_delta^2),
    data = start_data
  )
  start_gamma <- 0.03
  quadratic_coef <- stats::coef(fit_quadratic_start)
  start_beta <- max(2 * quadratic_coef["I(rep_delta^2)"] / start_gamma, 1e-04)
  start_alpha <- quadratic_coef["rep_delta"] - start_beta
  list(
    log_L_ref = unname(quadratic_coef["(Intercept)"]),
    alpha = unname(start_alpha),
    log_beta = log(unname(start_beta)),
    log_gamma = log(start_gamma)
  )
}

#' @keywords internal
#' @noRd
.ode_predict_from_coef <- function(rep_percent, model_coef, ode_log_length) {
  ode_log_length(
    rep_percent = rep_percent,
    log_L_ref = model_coef["log_L_ref"],
    alpha = model_coef["alpha"],
    log_beta = model_coef["log_beta"],
    log_gamma = model_coef["log_gamma"]
  )
}

#' Gradient of the ODE log length w.r.t. the parameters (delta method)
#' @keywords internal
#' @noRd
.ode_gradient <- function(rep_percent, model_coef, rep_center) {
  beta <- exp(model_coef["log_beta"])
  gamma <- exp(model_coef["log_gamma"])
  rep_delta <- rep_percent - rep_center
  exp_term <- exp(gamma * rep_delta)
  gradient <- cbind(
    log_L_ref = rep(1, length(rep_percent)),
    alpha = rep_delta,
    log_beta = beta / gamma * expm1(gamma * rep_delta),
    log_gamma = beta * (rep_delta * exp_term - expm1(gamma * rep_delta) / gamma)
  )
  gradient[, names(model_coef), drop = FALSE]
}

#' rep_percent of minimal ODE length (before clamping)
#' @keywords internal
#' @noRd
.ode_turning_point <- function(alpha, beta, gamma, rep_center, lower) {
  if (alpha < 0 && -alpha / beta > 0) {
    rep_center + log(-alpha / beta) / gamma
  } else {
    lower
  }
}

#' Replication percentage needed to reach given track lengths
#'
#' Invert a fitted track-length model: for each length, the replication
#' percentage at which the predicted median track length equals it, and the
#' average length velocity between consecutive lengths,
#' `diff(track_length) / diff(rep_percent)`.
#'
#' For the linear model the right-hand root of the quadratic
#' `b2 x^2 + b1 x + b0 = log(length)` is returned (`NA` if there is none).
#' For the ODE model the root is searched with `stats::uniroot(tol = 1e-8)`
#' between the turning point (clamped to the data range) and the largest
#' observed `rep_percent` (`NA` if the length is out of range).
#'
#' @param fit Output of [fit_track_velocity_linear()] or
#'   [fit_track_velocity_ode()].
#' @param track_length Positive numeric vector of track lengths (bp).
#'
#' @return A `data.frame` with `track_length`, `rep_percent` and
#'   `length_velocity` (`NA` in the first row).
#'
#' @export
predict_rep_percent <- function(fit, track_length) {
  if (!inherits(fit, "scrt_track_velocity")) {
    stop("'fit' must come from fit_track_velocity_linear() or fit_track_velocity_ode().")
  }
  if (!is.numeric(track_length) || length(track_length) == 0L) {
    stop("'track_length' must be a non-empty numeric vector.")
  }
  rep_percent <- if (identical(fit$model, "linear")) {
    co <- stats::coef(fit$fit_lmrob)
    vapply(track_length, .inverse_predict_linear, numeric(1),
           b0 = co["(Intercept)"], b1 = co["rep_percent"], b2 = co["I(rep_percent^2)"])
  } else {
    co <- stats::coef(fit$fit_ode)
    vapply(track_length, .inverse_predict_ode, numeric(1),
           log_L_ref = co["log_L_ref"], alpha = co["alpha"],
           beta = exp(co["log_beta"]), gamma = exp(co["log_gamma"]),
           rep_center = fit$rep_center,
           lower = fit$rep_range[1], upper = fit$rep_range[2])
  }
  out <- data.frame(track_length = track_length, rep_percent = unname(rep_percent))
  out$length_velocity <- NA_real_
  if (nrow(out) >= 2L) {
    for (j in 2:nrow(out)) {
      out$length_velocity[j] <- (out$track_length[j] - out$track_length[j - 1]) /
        (out$rep_percent[j] - out$rep_percent[j - 1])
    }
  }
  out
}

#' @keywords internal
#' @noRd
.inverse_predict_linear <- function(length, b0, b1, b2) {
  log_length <- log(length)
  a <- b2
  b <- b1
  c <- b0 - log_length
  delta <- b^2 - 4 * a * c
  if (delta < 0) {
    return(NA_real_)
  }
  unname((-b + sqrt(delta)) / (2 * a))
}

#' @keywords internal
#' @noRd
.inverse_predict_ode <- function(length, log_L_ref, alpha, beta, gamma, rep_center,
                                 lower, upper) {
  if (!is.finite(length) || length <= 0) {
    return(NA_real_)
  }
  ode_log_length <- .ode_log_length_fun(rep_center)
  turning_point <- .ode_turning_point(alpha, beta, gamma, rep_center, lower)
  lower_right <- min(max(turning_point, lower), upper)
  objective <- function(rep_percent) {
    ode_log_length(rep_percent, log_L_ref, alpha, log(beta), log(gamma)) - log(length)
  }
  objective_lower <- objective(lower_right)
  objective_upper <- objective(upper)
  if (!is.finite(objective_lower) || !is.finite(objective_upper)) {
    return(NA_real_)
  }
  if (abs(objective_lower) < .Machine$double.eps^0.5) {
    return(unname(lower_right))
  }
  if (objective_lower * objective_upper > 0) {
    return(NA_real_)
  }
  stats::uniroot(objective, interval = c(lower_right, upper), tol = 1e-08)$root
}
