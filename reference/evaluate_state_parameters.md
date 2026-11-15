# Evaluate replication-state parameters on a grid

For every row of `param_grid`, run
[`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md)
followed by
[`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md)
and summarise how many cells end up with a replicated fraction inside
`target_range`. Optionally the evaluation is run on a random pilot
subset of `n_cells` cells. A failing combination is recorded with
`status = "error"` instead of stopping the whole grid. Use
[`rank_state_parameters()`](https://midimagi.github.io/scRTtools/reference/rank_state_parameters.md)
to order the result.

## Usage

``` r
evaluate_state_parameters(
  g1bins_matrix,
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
  verbose = FALSE
)
```

## Arguments

- g1bins_matrix, g1_granges, genome_granges:

  As in
  [`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md).

- meta_data:

  A `data.frame` whose row names include all `colnames(g1bins_matrix)`;
  it is reordered to the matrix columns.

- param_grid:

  A `data.frame` with one row per combination; each column must be an
  argument of
  [`infer_single_cell_states()`](https://midimagi.github.io/scRTtools/reference/infer_single_cell_states.md)
  (e.g. `init_ratio_lower` and `final_ratio_lower`).

- target_range:

  Numeric length-2 range of `rep_percent` counted as on target. We use
  `c(0.001, 0.2)`.

- infer_args:

  Named list of fixed arguments of
  [`infer_single_cell_states()`](https://midimagi.github.io/scRTtools/reference/infer_single_cell_states.md);
  columns of `param_grid` take precedence. We use
  `list(hmm_aggr_width = 15L)`.

- assign_args:

  Named list of further arguments of
  [`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md)
  (e.g. `gap_max_bp`, `chunk_size`).

- qc_args:

  Named list of arguments of
  [`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md)
  (`min_rep_percent`, `max_sigma`, `neighbor_window`,
  `outlier_threshold`). We use
  `list(min_rep_percent = 0.01, max_sigma = 0.3, neighbor_window = 5, outlier_threshold = 0.2)`.

- n_cells:

  `NULL` (all cells) or the number of pilot cells sampled without
  replacement (matrix order is kept), e.g. `500`.

- seed:

  `NULL` or an integer seed for the pilot sampling.

- keep_results:

  Logical; also return the
  [`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md)
  results (e.g. to plot `rep_percent` histograms).

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html),
  passed to both
  [`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md)
  and
  [`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md).

- verbose:

  Logical; print one message per combination. Default `FALSE`.

## Value

A `data.frame`: the columns of `param_grid`, `status` (`"ok"` or
`"error"`), `error` and the summary columns above. With
`keep_results = TRUE`, the list of
[`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md)
results (`NULL` for failed rows) is attached as attribute `"results"`.
The pilot cells are attached as attribute `"cells"`.

## Details

Summary columns, with `in_target` = `rep_percent` within `target_range`
(inclusive) and `pass` = `pass_qc` (`NA` counts as `FALSE`): `n_total`
(cells), `n_non_missing` (non-`NA` `rep_percent`), `n_target`,
`n_pass_qc`, `target_prop_all_cells` = `n_target / n_total`,
`target_prop_non_missing` = `n_target / n_non_missing`,
`target_prop_pass_qc` = `sum(in_target & pass) / n_pass_qc`,
`non_missing_prop`, `pass_qc_prop` and `runtime_min`.

## See also

[`rank_state_parameters()`](https://midimagi.github.io/scRTtools/reference/rank_state_parameters.md)
