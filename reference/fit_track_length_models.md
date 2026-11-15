# Fit track-length velocity models per group

Run
[`fit_track_velocity_linear()`](https://midimagi.github.io/scRTtools/reference/fit_track_velocity_linear.md)
separately for each level of `group_col` and stack the result tables.
Group `i` uses hold-out seed `seed + i`.

## Usage

``` r
fit_track_length_models(
  data,
  group_col = "state",
  levels = NULL,
  seed = 2026L,
  n_grid = 200L,
  ...
)
```

## Arguments

- data:

  Output of
  [`summarize_track_lengths()`](https://midimagi.github.io/scRTtools/reference/summarize_track_lengths.md)
  with a grouping column.

- group_col:

  Grouping column.

- levels:

  Group levels to fit (in this order); default: sorted unique values.

- seed:

  Base seed of the hold-out splits.

- n_grid:

  Number of prediction points per group (e.g. `200`).

- ...:

  Further arguments of
  [`fit_track_velocity_linear()`](https://midimagi.github.io/scRTtools/reference/fit_track_velocity_linear.md)
  (`level`, `quantile`, `holdout_prop`, `fit_seed`,
  `positive_velocity_only`, `min_cells`, `min_distinct`, `sort_by_rep`).

## Value

A list with `state_models` (per-group results), `has_quantreg`, and the
stacked tables `model_data`, `model_prediction`, `linear_residuals`,
`robust_residuals`, `validation`, `validation_metrics` and
`model_comparison`, each with `group_col` as a factor.
