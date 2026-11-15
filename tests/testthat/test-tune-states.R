tune_inputs <- function() {
  x <- synth_g1(n_cells = 20L)
  gb <- synth_genome_bins(x$gr)
  md <- data.frame(row.names = rev(colnames(x$m)), cell_id = rev(colnames(x$m)))
  list(x = x, gb = gb, md = md)
}

test_that("evaluate_state_parameters summarises each grid row", {
  d <- tune_inputs()
  grid <- data.frame(init_ratio_lower = c(1.2, 1.2, 3), final_ratio_lower = c(1.2, 1.4, 1.4))
  ia <- list(hmm_aggr_width = 15L)
  qa <- list(min_rep_percent = 0.01, max_sigma = 0.3, neighbor_window = 5,
             outlier_threshold = 0.2)
  res <- evaluate_state_parameters(d$x$m, d$x$gr, d$gb, d$md, grid, infer_args = ia,
                                   qc_args = qa, keep_results = TRUE, verbose = FALSE)
  expect_s3_class(res, "data.frame")
  expect_identical(nrow(res), 3L)
  expect_identical(res$status, c("ok", "ok", "error"))
  expect_match(res$error[3], "init_ratio_upper")
  expect_identical(res$n_total, c(20L, 20L, 20L))
  # row 2 equals running the pipeline by hand
  a <- assign_replication_states(d$x$m, d$x$gr, d$gb,
                                 infer_args = list(hmm_aggr_width = 15L, init_ratio_lower = 1.2,
                                                   final_ratio_lower = 1.4),
                                 verbose = FALSE)
  r <- do.call(compute_rep_percent, c(list(a$genobins_rep_states, d$gb,
                                           d$md[colnames(d$x$m), , drop = FALSE],
                                           verbose = FALSE), qa))
  expect_identical(attr(res, "results")[[2]], r)
  expect_null(attr(res, "results")[[3]])
  rp <- r$meta_data$rep_percent
  expect_identical(res$n_target[2], sum(!is.na(rp) & rp >= 0.001 & rp <= 0.2))
  expect_identical(res$n_pass_qc[2], sum(r$meta_data$pass_qc))
  # pilot subset: reproducible, keeps matrix order
  p1 <- evaluate_state_parameters(d$x$m, d$x$gr, d$gb, d$md, grid[1, ], infer_args = ia,
                                  n_cells = 8, seed = 3, verbose = FALSE)
  p2 <- evaluate_state_parameters(d$x$m, d$x$gr, d$gb, d$md, grid[1, ], infer_args = ia,
                                  n_cells = 8, seed = 3, verbose = FALSE)
  expect_identical(p1$runtime_min > -1, TRUE)
  expect_identical(attr(p1, "cells"), attr(p2, "cells"))
  expect_identical(attr(p1, "cells"), colnames(d$x$m)[colnames(d$x$m) %in% attr(p1, "cells")])
  expect_identical(p1$n_total, 8L)
  # validation
  expect_error(evaluate_state_parameters(d$x$m, d$x$gr, d$gb, d$md, data.frame(foo = 1)),
               "foo")
  expect_error(evaluate_state_parameters(d$x$m, d$x$gr, d$gb, d$md[1:5, , drop = FALSE], grid),
               "missing from meta_data")
  expect_error(evaluate_state_parameters(d$x$m, d$x$gr, d$gb, d$md, grid,
                                         qc_args = list(bpparam = 1)), "qc_args")
})

test_that("summary and ranking match the notebook code", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tibble")
  `%>%` <- dplyr::`%>%`
  d <- tune_inputs()
  st <- synth_states()
  rep_result <- compute_rep_percent(st$st, st$gr, st$md, 0.01, 0.3, verbose = FALSE)
  rep_result$meta_data$rep_percent[c(2, 7)] <- NA
  # notebook summarize_parameter_result, verbatim apart from namespacing
  args <- list(target_lower = 0.001, target_upper = 0.2)
  summarize_parameter_result <- function(rep_result, init_ratio_lower, final_ratio_lower,
                                         stage, runtime_min) {
    result_meta <- rep_result$meta_data
    rep_percent <- result_meta$rep_percent
    pass_qc <- !is.na(result_meta$pass_qc) & result_meta$pass_qc
    in_target <- !is.na(rep_percent) & rep_percent >= args$target_lower &
      rep_percent <= args$target_upper
    n_total <- length(rep_percent)
    n_non_missing <- sum(!is.na(rep_percent))
    n_target <- sum(in_target)
    n_pass_qc <- sum(pass_qc)
    tibble::tibble(
      stage = stage, init_ratio_lower = init_ratio_lower,
      final_ratio_lower = final_ratio_lower, status = "ok", error = NA_character_,
      n_total = n_total, n_non_missing = n_non_missing, n_target = n_target,
      n_pass_qc = n_pass_qc, target_prop_all_cells = n_target / n_total,
      target_prop_non_missing = ifelse(n_non_missing > 0, n_target / n_non_missing, NA_real_),
      target_prop_pass_qc = ifelse(n_pass_qc > 0, sum(in_target & pass_qc) / n_pass_qc, NA_real_),
      non_missing_prop = n_non_missing / n_total, pass_qc_prop = n_pass_qc / n_total,
      runtime_min = runtime_min
    )
  }
  old <- as.data.frame(summarize_parameter_result(rep_result, 1.2, 1.4, "pilot", 0.5))
  new <- scRTtools:::.summarize_state_result(rep_result$meta_data, c(0.001, 0.2), 0.5)
  expect_identical(new, old[, names(new)])

  # notebook pilot_ranking, on a grid with ties and NAs
  withr::local_seed(4)
  res <- data.frame(
    init_ratio_lower = rep(c(1.0, 1.1, 1.2, 1.3), 3),
    final_ratio_lower = rep(c(1.2, 1.3, 1.4), each = 4),
    status = c(rep("ok", 10), "error", "ok"),
    target_prop_all_cells = sample(c(0.5, 0.6, NA), 12, TRUE),
    target_prop_non_missing = sample(c(0.7, 0.8), 12, TRUE),
    pass_qc_prop = sample(c(0.2, 0.3), 12, TRUE)
  )
  old <- res %>%
    dplyr::filter(status == "ok") %>%
    dplyr::mutate(distance_to_current = abs(init_ratio_lower - 1.2) + abs(final_ratio_lower - 1.4)) %>%
    dplyr::arrange(dplyr::desc(target_prop_all_cells), dplyr::desc(target_prop_non_missing),
                   dplyr::desc(pass_qc_prop), distance_to_current, init_ratio_lower,
                   final_ratio_lower)
  new <- rank_state_parameters(res, c(init_ratio_lower = 1.2, final_ratio_lower = 1.4))
  expect_identical(new, old)
  expect_error(rank_state_parameters(res[11, ], c(init_ratio_lower = 1.2)), "failed")
  expect_error(rank_state_parameters(res, c(nope = 1)), "reference")
})

test_that("pilot sampling matches the notebook", {
  d <- tune_inputs()
  set.seed(20260811L)
  sampled_cells <- sample(colnames(d$x$m), size = min(8L, ncol(d$x$m)), replace = FALSE)
  pilot_cells <- colnames(d$x$m)[colnames(d$x$m) %in% sampled_cells]
  res <- evaluate_state_parameters(d$x$m, d$x$gr, d$gb, d$md,
                                   data.frame(init_ratio_lower = 1.2), n_cells = 8L,
                                   seed = 20260811L, verbose = FALSE)
  expect_identical(attr(res, "cells"), pilot_cells)
})
