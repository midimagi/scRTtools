# Parameter tuning for assign_replication_states().
# Extracted from the notebooks
# jupyter_downsample_frag_quantile/1.4.tune_{init,final}_ratio_assign_replication_states_scRT6_rep1
# (summarize_parameter_result, evaluate_parameters and the pilot ranking cell).
# Plotting and checkpoint files of the notebooks are not part of the package.

#' Evaluate replication-state parameters on a grid
#'
#' For every row of `param_grid`, run [assign_replication_states()] followed by
#' [compute_rep_percent()] and summarise how many cells end up with a
#' replicated fraction inside `target_range`. Optionally the evaluation is run
#' on a random pilot subset of `n_cells` cells. A failing combination is
#' recorded with `status = "error"` instead of stopping the whole grid. Use
#' [rank_state_parameters()] to order the result.
#'
#' Summary columns, with `in_target` = `rep_percent` within `target_range`
#' (inclusive) and `pass` = `pass_qc` (`NA` counts as `FALSE`):
#' `n_total` (cells), `n_non_missing` (non-`NA` `rep_percent`), `n_target`,
#' `n_pass_qc`, `target_prop_all_cells` = `n_target / n_total`,
#' `target_prop_non_missing` = `n_target / n_non_missing`,
#' `target_prop_pass_qc` = `sum(in_target & pass) / n_pass_qc`,
#' `non_missing_prop`, `pass_qc_prop` and `runtime_min`.
#'
#' @param g1bins_matrix,g1_granges,genome_granges As in
#'   [assign_replication_states()].
#' @param meta_data A `data.frame` whose row names include all
#'   `colnames(g1bins_matrix)`; it is reordered to the matrix columns.
#' @param param_grid A `data.frame` with one row per combination; each column
#'   must be an argument of [infer_single_cell_states()] (e.g.
#'   `init_ratio_lower` and `final_ratio_lower`).
#' @param target_range Numeric length-2 range of `rep_percent` counted as
#'   on target. We use `c(0.001, 0.2)`.
#' @param infer_args Named list of fixed arguments of
#'   [infer_single_cell_states()]; columns of `param_grid` take precedence.
#'   We use `list(hmm_aggr_width = 15L)`.
#' @param assign_args Named list of further arguments of
#'   [assign_replication_states()] (e.g. `gap_max_bp`, `chunk_size`).
#' @param qc_args Named list of arguments of [compute_rep_percent()]
#'   (`min_rep_percent`, `max_sigma`, `neighbor_window`,
#'   `outlier_threshold`). We use
#'   `list(min_rep_percent = 0.01, max_sigma = 0.3, neighbor_window = 5,
#'   outlier_threshold = 0.2)`.
#' @param n_cells `NULL` (all cells) or the number of pilot cells sampled
#'   without replacement (matrix order is kept), e.g. `500`.
#' @param seed `NULL` or an integer seed for the pilot sampling.
#' @param keep_results Logical; also return the [compute_rep_percent()]
#'   results (e.g. to plot `rep_percent` histograms).
#' @param BPPARAM A `BiocParallel::BiocParallelParam`, passed to both
#'   [assign_replication_states()] and [compute_rep_percent()].
#' @param verbose Logical; print one message per combination. Default `FALSE`.
#'
#' @return A `data.frame`: the columns of `param_grid`, `status` (`"ok"` or
#'   `"error"`), `error` and the summary columns above. With
#'   `keep_results = TRUE`, the list of [compute_rep_percent()] results (`NULL`
#'   for failed rows) is attached as attribute `"results"`. The pilot cells
#'   are attached as attribute `"cells"`.
#'
#' @seealso [rank_state_parameters()]
#' @export
evaluate_state_parameters <- function(g1bins_matrix,
                                      g1_granges,
                                      genome_granges,
                                      meta_data,
                                      param_grid,
                                      target_range = c(0.001, 0.2),
                                      infer_args = list(),
                                      assign_args = list(),
                                      qc_args = list(),
                                      n_cells = NULL,
                                      seed = NULL,
                                      keep_results = FALSE,
                                      BPPARAM = BiocParallel::SerialParam(),
                                      verbose = FALSE) {
  # ---------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------
  if (is.null(colnames(g1bins_matrix))) {
    stop("colnames(g1bins_matrix) must not be NULL.")
  }
  if (!is.data.frame(meta_data)) {
    stop("'meta_data' must be a data.frame.")
  }
  if (!all(colnames(g1bins_matrix) %in% rownames(meta_data))) {
    stop("Cells in g1bins_matrix are missing from meta_data.")
  }
  if (!is.data.frame(param_grid) || nrow(param_grid) == 0L || ncol(param_grid) == 0L) {
    stop("'param_grid' must be a non-empty data.frame.")
  }
  infer_formals <- setdiff(names(formals(infer_single_cell_states)),
                           c("counts", "breakpoints"))
  bad <- setdiff(names(param_grid), infer_formals)
  if (length(bad) > 0L) {
    stop("Columns of 'param_grid' are not arguments of infer_single_cell_states(): ",
         paste(bad, collapse = ", "))
  }
  if (!is.numeric(target_range) || length(target_range) != 2L || anyNA(target_range) ||
      target_range[1] > target_range[2]) {
    stop("'target_range' must be an increasing numeric vector of length 2.")
  }
  for (nm in c("infer_args", "assign_args", "qc_args")) {
    a <- get(nm)
    if (!is.list(a) || (length(a) > 0L && (is.null(names(a)) || any(!nzchar(names(a)))))) {
      stop("'", nm, "' must be a named list.")
    }
  }
  bad <- setdiff(names(infer_args), infer_formals)
  if (length(bad) > 0L) {
    stop("Unknown infer_args: ", paste(bad, collapse = ", "))
  }
  reserved <- c("g1bins_matrix", "g1_granges", "genome_granges", "infer_args",
                "BPPARAM", "verbose")
  bad <- intersect(names(assign_args), reserved)
  if (length(bad) > 0L) {
    stop("'assign_args' must not contain: ", paste(bad, collapse = ", "))
  }
  bad <- setdiff(names(qc_args), c("min_rep_percent", "max_sigma",
                                   "neighbor_window", "outlier_threshold"))
  if (length(bad) > 0L) {
    stop("Unknown qc_args: ", paste(bad, collapse = ", "))
  }
  if (!is.null(n_cells) &&
      (!is.numeric(n_cells) || length(n_cells) != 1L || is.na(n_cells) || n_cells < 1)) {
    stop("'n_cells' must be NULL or a positive integer.")
  }

  # ---------------------------------------------------------------
  # Cell selection (pilot subset keeps the matrix order)
  # ---------------------------------------------------------------
  meta_data <- meta_data[colnames(g1bins_matrix), , drop = FALSE]
  cells <- colnames(g1bins_matrix)
  if (!is.null(n_cells)) {
    sample_cells <- function() {
      sampled <- sample(cells, size = min(n_cells, length(cells)), replace = FALSE)
      cells[cells %in% sampled]
    }
    cells <- if (is.null(seed)) sample_cells() else withr::with_seed(seed, sample_cells())
    g1bins_matrix <- g1bins_matrix[, cells, drop = FALSE]
    meta_data <- meta_data[cells, , drop = FALSE]
  }

  # ---------------------------------------------------------------
  # Evaluate each combination
  # ---------------------------------------------------------------
  results <- vector("list", nrow(param_grid))
  rows <- vector("list", nrow(param_grid))

  for (i in seq_len(nrow(param_grid))) {
    params <- as.list(param_grid[i, , drop = FALSE])
    .vmsg(verbose, "Evaluating ",
          paste(names(params), unlist(params), sep = " = ", collapse = ", "))
    ia <- utils::modifyList(infer_args, params)
    start_time <- Sys.time()

    rows[[i]] <- tryCatch({
      assigned <- do.call(assign_replication_states, c(
        list(g1bins_matrix = g1bins_matrix, g1_granges = g1_granges,
             genome_granges = genome_granges, infer_args = ia,
             BPPARAM = BPPARAM, verbose = FALSE),
        assign_args
      ))
      rep_result <- do.call(compute_rep_percent, c(
        list(genobins_rep_states = assigned$genobins_rep_states,
             genome_granges = genome_granges, meta_data = meta_data,
             BPPARAM = BPPARAM, verbose = FALSE),
        qc_args
      ))
      runtime_min <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
      if (isTRUE(keep_results)) results[[i]] <- rep_result
      .summarize_state_result(rep_result$meta_data, target_range, runtime_min)
    }, error = function(e) {
      .summarize_state_error(conditionMessage(e), ncol(g1bins_matrix),
                             as.numeric(difftime(Sys.time(), start_time, units = "mins")))
    })
  }

  out <- cbind(param_grid, do.call(rbind, rows))
  rownames(out) <- NULL
  if (isTRUE(keep_results)) attr(out, "results") <- results
  attr(out, "cells") <- cells
  out
}

#' Summary row of one successful evaluation
#' @keywords internal
#' @noRd
.summarize_state_result <- function(result_meta, target_range, runtime_min) {
  rep_percent <- result_meta$rep_percent
  pass_qc <- !is.na(result_meta$pass_qc) & result_meta$pass_qc

  in_target <-
    !is.na(rep_percent) &
    rep_percent >= target_range[1] &
    rep_percent <= target_range[2]

  n_total <- length(rep_percent)
  n_non_missing <- sum(!is.na(rep_percent))
  n_target <- sum(in_target)
  n_pass_qc <- sum(pass_qc)

  data.frame(
    status = "ok",
    error = NA_character_,
    n_total = n_total,
    n_non_missing = n_non_missing,
    n_target = n_target,
    n_pass_qc = n_pass_qc,
    target_prop_all_cells = n_target / n_total,
    target_prop_non_missing = ifelse(n_non_missing > 0, n_target / n_non_missing,
                                     NA_real_),
    target_prop_pass_qc = ifelse(n_pass_qc > 0, sum(in_target & pass_qc) / n_pass_qc,
                                 NA_real_),
    non_missing_prop = n_non_missing / n_total,
    pass_qc_prop = n_pass_qc / n_total,
    runtime_min = runtime_min,
    stringsAsFactors = FALSE
  )
}

#' Summary row of one failed evaluation
#' @keywords internal
#' @noRd
.summarize_state_error <- function(message, n_total, runtime_min) {
  data.frame(
    status = "error",
    error = message,
    n_total = n_total,
    n_non_missing = NA_integer_,
    n_target = NA_integer_,
    n_pass_qc = NA_integer_,
    target_prop_all_cells = NA_real_,
    target_prop_non_missing = NA_real_,
    target_prop_pass_qc = NA_real_,
    non_missing_prop = NA_real_,
    pass_qc_prop = NA_real_,
    runtime_min = runtime_min,
    stringsAsFactors = FALSE
  )
}

#' Rank evaluated replication-state parameters
#'
#' Keep the successful rows of [evaluate_state_parameters()] and order them by
#' decreasing `target_prop_all_cells`, `target_prop_non_missing` and
#' `pass_qc_prop`, then by increasing `distance_to_current` (sum of absolute
#' differences to `reference`) and finally by the parameter columns named in
#' `reference` (increasing). `NA` values sort last.
#'
#' @param results Output of [evaluate_state_parameters()].
#' @param reference Named numeric vector (or list) of the current parameter
#'   values. We use
#'   `c(init_ratio_lower = 1.2, final_ratio_lower = 1.4)`.
#'
#' @return `results` restricted to `status == "ok"`, with a
#'   `distance_to_current` column, in ranked order. An error is raised if no
#'   row succeeded.
#'
#' @export
rank_state_parameters <- function(results, reference) {
  needed <- c("status", "target_prop_all_cells", "target_prop_non_missing",
              "pass_qc_prop")
  if (!is.data.frame(results) || !all(needed %in% names(results))) {
    stop("'results' must be the output of evaluate_state_parameters().")
  }
  reference <- unlist(reference)
  if (!is.numeric(reference) || length(reference) == 0L || is.null(names(reference)) ||
      !all(names(reference) %in% names(results))) {
    stop("'reference' must be a named numeric vector of parameter columns of 'results'.")
  }

  ok <- results[!is.na(results$status) & results$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0L) {
    stop("All parameter combinations failed.")
  }

  distance <- 0
  for (nm in names(reference)) {
    distance <- distance + abs(ok[[nm]] - reference[[nm]])
  }
  ok$distance_to_current <- distance

  ord <- do.call(order, c(
    list(-ok$target_prop_all_cells, -ok$target_prop_non_missing, -ok$pass_qc_prop,
         ok$distance_to_current),
    unname(as.list(ok[names(reference)]))
  ))
  out <- ok[ord, , drop = FALSE]
  rownames(out) <- NULL
  out
}
