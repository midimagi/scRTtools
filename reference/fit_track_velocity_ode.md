# Robust ODE model of track length versus replication percentage

Fit the closed-form solution of the ODE
`d log L / dx = alpha + beta * exp(gamma * (x - rep_center))`,
`log L(x) = log_L_ref + alpha * d + beta / gamma * expm1(gamma * d)`
with `d = x - rep_center`, by robust nonlinear regression
(`robustbase::nlrob(method = "M")`, `beta` and `gamma` on the log scale)
and optionally by median regression
([`quantreg::nlrq()`](https://rdrr.io/pkg/quantreg/man/nlrq.html),
started at the robust estimates). Starting values come from a robust
quadratic fit (`gamma = 0.03`, `beta = max(2 * b2 / gamma, 1e-4)`,
`alpha = b1 - beta`).

## Usage

``` r
fit_track_velocity_ode(
  data,
  rep_center = NULL,
  start = NULL,
  n_grid = 100L,
  level = 0.95,
  quantile = TRUE,
  holdout_prop = 0.5,
  seed = 2026L,
  fit_seed = NULL,
  maxiter = 500L
)
```

## Arguments

- data:

  Output of
  [`summarize_track_lengths()`](https://midimagi.github.io/scRTtools/reference/summarize_track_lengths.md)
  (`rep_percent`, `median_length`, `log_length`).

- rep_center:

  Centre of the parameterisation; default: median of `data$rep_percent`.

- start:

  `NULL` (computed) or a named list of starting values (`log_L_ref`,
  `alpha`, `log_beta`, `log_gamma`).

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

- maxiter:

  `maxiter` of
  [`stats::nls.control()`](https://rdrr.io/r/stats/nls.control.html) for
  `nlrob()`.

## Value

An object of class `scrt_track_velocity` (a list) with `model = "ode"`,
`rep_center`, `start`, `fit_ode`, `fit_nlrq`, `parameters`
(`rep_center`, `log_L_ref`, `alpha`, `beta`, `gamma`), `data` (with
`is_outlier`), `residuals`, `prediction` (all grid rows), `velocity`
(rows with positive velocity), `validation`, `validation_metrics`,
`rse`, `pseudo_r2`, `turning_rep_percent`, `minimal_predicted_length`
and `rep_range`.

## Details

Predictions on an even grid get delta-method confidence intervals and
prediction intervals with `qnorm((1 + level) / 2)`; the length velocity
is `L(x) * (alpha + beta * exp(gamma * d))`. A random hold-out split
refits the model (with `rep_center` of the full data) and reports RMSE,
pseudo-R2 and correlation on the test cells.
