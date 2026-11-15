test_that("summarize_track_lengths matches the notebook dplyr code", {
  skip_if_not_installed("dplyr")
  `%>%` <- dplyr::`%>%`
  tt <- synth_track_table(n_cells = 30L)
  tt <- tt[sample(nrow(tt)), ]
  old <- tt %>% dplyr::group_by(cell_id, state) %>%
    dplyr::summarise(median_length = stats::median(width_bp, na.rm = TRUE),
                     log_length = log(median_length),
                     rep_percent = dplyr::first(rep_percent) * 100, .groups = "drop")
  expect_identical(summarize_track_lengths(tt, by = "state"), as.data.frame(old))
  old1 <- tt %>% dplyr::group_by(cell_id) %>%
    dplyr::summarise(median_length = stats::median(width_bp), log_length = log(median_length),
                     rep_percent = rep_percent[1] * 100, .groups = "drop")
  expect_identical(summarize_track_lengths(tt), as.data.frame(old1))
  expect_identical(nrow(summarize_track_lengths(tt[0, ])), 0L)
  expect_error(summarize_track_lengths(tt[, names(tt) != "cell_id"]), "cell_id")
})

test_that("fit_track_velocity_linear fits, validates and inverts", {
  skip_if_not_installed("robustbase")
  skip_if_not_installed("quantreg")
  d <- summarize_track_lengths(synth_track_table(states = "ctr"))
  fit <- fit_track_velocity_linear(d, fit_seed = 1)
  expect_s3_class(fit, "scrt_track_velocity")
  expect_identical(fit$model, "linear")
  expect_identical(fit, fit_track_velocity_linear(d, fit_seed = 1))
  p <- fit$prediction
  expect_true(all(p$length_velocity > 0))
  expect_true(all(p$lmrob_length_lower_ci <= p$lmrob_length & p$lmrob_length <= p$lmrob_length_upper_ci))
  expect_true(all(p$lmrob_length_lower_pi <= p$lmrob_length_lower_ci))
  expect_true(any(fit$data$is_outlier))
  expect_gt(fit$validation_metrics$pseudo_r2, 0.5)
  expect_equal(nrow(fit$validation), nrow(d) - floor(nrow(d) / 2))
  expect_identical(nrow(fit$model_comparison), 2L * nrow(p))
  inv <- predict_rep_percent(fit, c(2e5, 4e5, 8e5))
  co <- stats::coef(fit$fit_lmrob)
  expect_equal(unname(co[1] + co[2] * inv$rep_percent + co[3] * inv$rep_percent^2),
               log(c(2e5, 4e5, 8e5)))
  expect_true(is.na(inv$length_velocity[1]))
  nq <- fit_track_velocity_linear(d, quantile = FALSE, positive_velocity_only = FALSE, n_grid = 50)
  expect_null(nq$fit_rq)
  expect_identical(nrow(nq$prediction), 50L)
  expect_error(fit_track_velocity_linear(d[1:10, ]), "too few")
})

test_that("fit_track_velocity_linear(sort_by_rep = FALSE) matches notebook 2.3", {
  skip_if_not_installed("robustbase")
  skip_if_not_installed("quantreg")
  d <- summarize_track_lengths(synth_track_table(states = "ctr"))
  expect_true(is.unsorted(d$rep_percent))   # input order differs from rep_percent order
  # notebook 2.3 inline code
  set.seed(5)
  fit_linear <- lm(log_length ~ rep_percent, data = d)
  fit_lmrob <- robustbase::lmrob(log_length ~ rep_percent + I(rep_percent^2), data = d)
  fit_rq <- quantreg::rq(log_length ~ rep_percent + I(rep_percent^2), tau = 0.5, data = d)
  df_pred <- data.frame(rep_percent = seq(min(d$rep_percent), max(d$rep_percent), length.out = 100))
  lmrob_pred <- predict(fit_lmrob, newdata = df_pred, se.fit = TRUE)
  rq_pred <- predict(fit_rq, newdata = df_pred, se.fit = TRUE)
  crit <- qt(0.975, df = fit_lmrob$df.residual)
  half_CI <- crit * lmrob_pred$se.fit
  half_PI <- crit * sqrt(lmrob_pred$se.fit^2 + lmrob_pred$residual.scale^2)
  b1 <- coef(fit_lmrob)["rep_percent"]; b2 <- coef(fit_lmrob)["I(rep_percent^2)"]
  velocity <- exp(lmrob_pred$fit) * (b1 + 2 * b2 * df_pred$rep_percent)
  set.seed(2026)
  train_indices <- sample(1:nrow(d), size = 0.5 * nrow(d))
  fit_train <- robustbase::lmrob(log_length ~ rep_percent + I(rep_percent^2), data = d[train_indices, ])
  pred_test <- predict(fit_train, newdata = d[-train_indices, ])

  set.seed(5)
  new <- fit_track_velocity_linear(d, sort_by_rep = FALSE, positive_velocity_only = FALSE)
  p <- new$prediction
  expect_identical(coef(new$fit_lmrob), coef(fit_lmrob))
  expect_identical(coef(new$fit_linear), coef(fit_linear))
  expect_identical(coef(new$fit_rq), coef(fit_rq))
  expect_identical(p$lmrob_log, unname(lmrob_pred$fit))
  expect_identical(p$lmrob_log_lower_ci, unname(lmrob_pred$fit - half_CI))
  expect_identical(p$lmrob_length_upper_pi, unname(exp(lmrob_pred$fit + half_PI)))
  expect_identical(p$rq_length, unname(exp(rq_pred)))
  expect_identical(p$length_velocity, unname(velocity))
  expect_identical(new$data$cell_id, d$cell_id)
  expect_identical(new$data$is_outlier, unname(fit_lmrob$rweights == 0))
  expect_identical(new$validation$predicted, unname(pred_test))
  expect_identical(new$validation$actual, d$log_length[-train_indices])
  expect_identical(new$minimal_predicted_length,
                   unname(exp(coef(fit_lmrob)[1] - b1^2 / 4 / b2)))
  # default (sorted) still runs and differs only in order-dependent parts
  set.seed(5)
  srt <- fit_track_velocity_linear(d, positive_velocity_only = FALSE)
  expect_false(is.unsorted(srt$data$rep_percent))
  expect_equal(unname(coef(srt$fit_lmrob)), unname(coef(fit_lmrob)), tolerance = 1e-4)
})

test_that("fit_track_length_models matches the fig. 6 helpers", {
  skip_on_cran()
  skip_if_not_installed("robustbase")
  skip_if_not_installed("quantreg")
  skip_if_not_installed("tidyr")
  helper <- file.path(local_dir("SCRTTOOLS_NOTEBOOKS_DS"), "retreat_figure6_helpers.R")
  skip_if_not(file.exists(helper), "retreat_figure6_helpers.R not available")
  env <- new.env(parent = asNamespace("dplyr"))
  env$lmrob <- robustbase::lmrob
  env$pivot_longer <- tidyr::pivot_longer
  sys.source(helper, envir = env, toplevel.env = env)
  d <- summarize_track_lengths(synth_track_table(), by = "state")
  strip <- function(df) {
    df <- as.data.frame(df)
    df[] <- lapply(df, function(v) if (is.factor(v)) v else unname(v))
    rownames(df) <- NULL
    df
  }

  set.seed(11)
  old <- env$fit_track_length_models(d, c("ctr", "ama"))
  set.seed(11)
  new <- fit_track_length_models(d, levels = c("ctr", "ama"))
  expect_identical(names(new), names(old))
  # later groups: lmrob() subsampling sees a different random-number state
  # (the original calls set.seed() without restoring it), so the fits agree
  # only up to lmrob's convergence tolerance
  for (tab in setdiff(names(old), c("state_models", "has_quantreg"))) {
    expect_equal(strip(new[[tab]]), strip(old[[tab]]), tolerance = 1e-4, info = tab)
  }
  # the first group sees the same random-number stream: identical
  o1 <- old$state_models$ctr
  n1 <- new$state_models$ctr
  for (tab in c("prediction", "validation_metrics", "linear_residuals", "robust_residuals")) {
    expect_identical(strip(n1[[tab]]), strip(o1[[tab]]), info = tab)
  }
  expect_identical(stats::coef(n1$fit_lmrob), stats::coef(o1$fit_lmrob))
})

test_that("fit_track_velocity_ode matches the notebook code", {
  skip_on_cran()
  skip_if_not_installed("robustbase")
  skip_if_not_installed("quantreg")
  d <- summarize_track_lengths(synth_track_table(states = "ctr", ode = TRUE))
  nb <- notebook_ode_env(d)
  # the notebook defines these in the global environment, where nlrob() finds them
  assign("ode_log_length", nb$ode_log_length, envir = globalenv())
  assign("rep_center", nb$rep_center, envir = globalenv())
  withr::defer(rm(list = c("ode_log_length", "rep_center"), envir = globalenv()))
  set.seed(5)
  fit_ode <- nb$fit_ode_model(d, start = nb$get_ode_start(d))
  fit_nlrq <- quantreg::nlrq(
    log_length ~ ode_log_length(rep_percent, log_L_ref, alpha, log_beta, log_gamma),
    tau = 0.5, data = d, start = as.list(coef(fit_ode)), trace = FALSE
  )
  x_seq <- seq(min(d$rep_percent), max(d$rep_percent), length.out = 100)
  ode_pred_log <- nb$ode_predict_from_coef(x_seq, coef(fit_ode))
  g <- nb$ode_gradient(x_seq, coef(fit_ode))
  se <- sqrt(pmax(rowSums((g %*% vcov(fit_ode)) * g), 0))
  half_pi <- qnorm(0.975) * sqrt(se^2 + sigma(fit_ode)^2)
  set.seed(2026)
  tr <- sample(1:nrow(d), size = 0.5 * nrow(d))
  fit_train <- nb$fit_ode_model(d[tr, ])
  pred_test <- nb$ode_predict_from_coef(d$rep_percent[-tr], coef(fit_train))

  set.seed(5)
  new <- fit_track_velocity_ode(d)
  expect_identical(unname(coef(new$fit_ode)), unname(coef(fit_ode)))
  expect_identical(coef(new$fit_nlrq), coef(fit_nlrq))
  expect_identical(new$prediction$ode_log, unname(ode_pred_log))
  expect_identical(new$prediction$ode_log_upper_pi, unname(ode_pred_log + half_pi))
  expect_identical(new$validation$predicted, unname(pred_test))
  expect_identical(new$data$is_outlier, unname(weights(fit_ode, type = "robustness") == 0))
  expect_true(all(new$velocity$length_velocity > 0))
  expect_identical(new$parameters$parameter, c("rep_center", "log_L_ref", "alpha", "beta", "gamma"))

  # inverse prediction: notebook inverse_predict
  co <- coef(fit_ode)
  inv_nb <- function(length, lower = min(d$rep_percent), upper = max(d$rep_percent)) {
    alpha <- co["alpha"]; beta <- exp(co["log_beta"]); gamma <- exp(co["log_gamma"])
    turning_point <- if (alpha < 0 && -alpha / beta > 0) {
      nb$rep_center + log(-alpha / beta) / gamma
    } else lower
    lower_right <- min(max(turning_point, lower), upper)
    objective <- function(r) nb$ode_log_length(r, co["log_L_ref"], alpha, log(beta), log(gamma)) - log(length)
    lo <- objective(lower_right); up <- objective(upper)
    if (!is.finite(lo) || !is.finite(up)) return(NA_real_)
    if (abs(lo) < .Machine$double.eps^0.5) return(lower_right)
    if (lo * up > 0) return(NA_real_)
    uniroot(objective, interval = c(lower_right, upper), tol = 1e-08)$root
  }
  L <- c(240000, 420000, 740000, 1420000, 2430000, 4470000)
  inv <- predict_rep_percent(new, L)
  expect_identical(inv$rep_percent, unname(vapply(L, inv_nb, numeric(1))))
  expect_error(predict_rep_percent(list(), 1), "fit_track_velocity")
})

test_that("count_track_centers counts filtered centres", {
  x <- synth_states()
  x$gr$is_mappable <- TRUE
  rp <- colSums(x$st == 4, na.rm = TRUE) / colSums(!is.na(x$st))
  tr <- define_replication_tracks(x$st, x$gr, rp, verbose = FALSE)
  bins <- GenomicRanges::tile(GenomicRanges::reduce(x$gr), width = 200000)
  bins <- BiocGenerics::unlist(bins)
  res <- count_track_centers(tr, bins)
  expect_identical(sum(res$score), sum(lengths(tr)))
  f <- count_track_centers(tr, bins, track_filter = (flag == 0) & (num_bins > 3))
  expect_identical(sum(f$score), sum(vapply(tr, function(g) sum(g$flag == 0 & g$num_bins > 3), integer(1))))
  thr <- 3
  f2 <- count_track_centers(tr, bins, track_filter = num_bins > thr)   # caller variable
  expect_identical(sum(f2$score), sum(vapply(tr, function(g) sum(g$num_bins > 3), integer(1))))
  expect_error(count_track_centers(tr, bins, track_filter = nope > 1), "track_filter")
  orig <- source_orig_file("analysis/count_track_centers.R")
  expect_identical(res, quiet_deprecated(orig$count_track_centers(tr, bins)))
  expect_identical(f, quiet_deprecated(
    orig$count_track_centers(tr, bins, track_filter = (flag == 0) & (num_bins > 3))))
})
