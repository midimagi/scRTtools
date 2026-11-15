# Robust quadratic model of track length versus replication percentage

Fit `log_length ~ rep_percent + I(rep_percent^2)` with
[`robustbase::lmrob()`](https://rdrr.io/pkg/robustbase/man/lmrob.html)
(plus [`lm()`](https://rdrr.io/r/stats/lm.html) with a linear term and,
optionally, a median regression with
[`quantreg::rq()`](https://rdrr.io/pkg/quantreg/man/rq.html)), predict
on an even grid with confidence and prediction intervals, derive the
length velocity `exp(fit) * (b1 + 2 b2 x)` and validate on a random
hold-out split.

## Usage

``` r
fit_track_velocity_linear(
  data,
  n_grid = 100L,
  level = 0.95,
  quantile = TRUE,
  holdout_prop = 0.5,
  seed = 2026L,
  fit_seed = NULL,
  positive_velocity_only = TRUE,
  min_cells = 20L,
  min_distinct = 8L,
  sort_by_rep = TRUE
)
```

## Arguments

- data:

  Output of
  [`summarize_track_lengths()`](https://midimagi.github.io/scRTtools/reference/summarize_track_lengths.md)
  (`rep_percent`, `median_length`, `log_length`).

- n_grid:

  Number of prediction points (e.g. `100` or `200`).

- level:

  Confidence level of the intervals.

- quantile:

  Logical; also fit the median regression (needs quantreg).

- holdout_prop:

  Fraction of cells used for training in the hold-out validation
  (`floor(holdout_prop * n)`).

- seed:

  Seed of the hold-out split.

- fit_seed:

  `NULL` or a seed for the full-data fits (see Details).

- positive_velocity_only:

  Logical; drop prediction rows with `length_velocity <= 0`.

- min_cells, min_distinct:

  Minimum number of cells and of distinct `rep_percent` values (e.g.
  `20` and `8`).

- sort_by_rep:

  Logical; sort the cells by `rep_percent` before fitting (default);
  `FALSE` keeps the input order. The order changes the hold-out split
  and the row order of `data`, and moves the `lmrob()` fits within its
  convergence tolerance.

## Value

An object of class `scrt_track_velocity` (a list) with
`model = "linear"`, the fits `fit_linear`, `fit_lmrob`, `fit_rq`, `data`
(with `is_outlier` = zero robustness weight), `prediction`,
`linear_residuals`, `robust_residuals`, `validation`,
`validation_metrics` (`rmse`, `pseudo_r2`, `correlation`),
`model_comparison` (long table of the lmrob / rq predicted lengths),
`minimal_predicted_length` (`exp(b0 - b1^2 / (4 b2))`) and `rep_range`.

## Details

Cells with non-finite `rep_percent` or `median_length`, or
`median_length <= 0`, are dropped and (with `sort_by_rep = TRUE`) the
data are sorted by `rep_percent`. Intervals use
`qt((1 + level) / 2, df.residual)` with the lmrob standard errors and
residual scale.

`lmrob()` draws random subsamples for its initial estimate; with
`fit_seed = NULL` the current random-number state is used, otherwise the
fits run under `withr::with_seed(fit_seed)`.

## See also

[`fit_track_length_models()`](https://midimagi.github.io/scRTtools/reference/fit_track_length_models.md),
[`fit_track_velocity_ode()`](https://midimagi.github.io/scRTtools/reference/fit_track_velocity_ode.md),
[`predict_rep_percent()`](https://midimagi.github.io/scRTtools/reference/predict_rep_percent.md)
