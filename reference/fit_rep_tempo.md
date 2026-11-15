# Fit a replication tempo curve to every region

Fits one curve per row of a
[`rep_tempo()`](https://midimagi.github.io/scRTtools/reference/rep_tempo.md)
table (regional replication percentage against genome-wide replication
progress) and extracts timing features from it.

## Usage

``` r
fit_rep_tempo(
  df_rep_tempo,
  rep_percent = NULL,
  model_type = c("sslogis", "gam_low", "gam_fast"),
  n_boot = 0L,
  conf_level = 0.95,
  seed = 1L,
  gam_k = 10L,
  gam_method = c("GCV.Cp", "REML"),
  gam_fixed_df = FALSE,
  anchors = TRUE,
  grid_points = 1001L,
  curve_ids = NULL,
  verbose = FALSE,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- df_rep_tempo:

  Output of
  [`rep_tempo()`](https://midimagi.github.io/scRTtools/reference/rep_tempo.md):
  numeric regions by progress groups, values in 0-100, with unique row
  names.

- rep_percent:

  Optional sorted progress values in 0-100, one per column, for a plain
  matrix without the `rep_bins` attribute.

- model_type:

  `"sslogis"` (default), `"gam_low"` or `"gam_fast"`.

- n_boot:

  Number of bootstrap replicates; 0 (default) disables the bootstrap.
  200 replicates take hours genome-wide.

- conf_level:

  Confidence level of the percentile intervals.

- seed:

  Non-negative integer seed for the bootstrap.

- gam_k:

  Maximum basis dimension of the GAM smooth (reduced automatically for
  regions with few points).

- gam_method:

  Smoothing-parameter method, `"GCV.Cp"` or `"REML"`.

- gam_fixed_df:

  If `TRUE`, fit an unpenalised regression spline with fixed degrees of
  freedom (`fx = TRUE`).

- anchors:

  Add the two fixed end-point observations to the SSlogis fit.

- grid_points:

  Number of grid points on 0-100 used for curves, AUC and grid-based
  thresholds.

- curve_ids:

  Region IDs whose fitted curves are returned in the `curves` attribute;
  `NULL` means all regions.

- verbose:

  Print progress messages. Default `FALSE`.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html);
  regions are fitted independently. Worker processes are limited to one
  BLAS/OpenMP thread.

## Value

A data.frame with one row per region: `region_id`, `model_type`,
`fit_status`, `t10`/`t25`/`t50`/`t75`/`t90` with their status columns,
`window_10_90`, `max_slope`, `auc`, `auc_normalized`, fit diagnostics
(`rmse`, `r_squared`, `deviance`, residual summaries), SSlogis
parameters, QC flags and bootstrap intervals. Attributes: `fitted` and
`residuals` (percent, at the input progress values), `curves`,
`rep_bins`, `range_info`, `settings` and `units`.

## Details

`"sslogis"` is an unweighted three-parameter logistic NLS
(`Asym * plogis((x - xmid) / scal)`) with two optional fixed end-point
observations (0 at 0 percent, 1 at 100 percent) and a `plinear`
fallback. `"gam_low"` and `"gam_fast"` fit the same Gaussian `mgcv` GAM
`y ~ s(x)` and differ only in how thresholds are extracted: `"gam_low"`
refines them with `uniroot`, `"gam_fast"` interpolates on a dense
prediction grid. GAM fits are neither bounded nor monotonic; QC flags
record this.

Thresholds (`t10` ... `t90`) are real increasing crossings inside each
region's observed progress range; thresholds that are not crossed, or
are crossed more than once, are `NA`. They are in genome-wide percentage
points, as is `window_10_90 = t90 - t10`.

The bootstrap resamples the observed (progress, percentage) points with
replacement and refits the model. It is a group-level interval, not a
cell-level one. A CI is reported only when at least
`max(20, ceiling(0.7 * n_boot))` replicates give a valid estimate. Row
`i` uses seed `seed + i - 1`, so serial and parallel runs agree. The
caller's RNG state is restored.

## See also

[`rep_tempo()`](https://midimagi.github.io/scRTtools/reference/rep_tempo.md),
[`subset_rep_tempo()`](https://midimagi.github.io/scRTtools/reference/subset_rep_tempo.md)
