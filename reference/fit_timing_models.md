# Logistic timing model of each initiation zone

For each IZ, fit `usage ~ progress` by logistic regression
([`stats::glm.fit`](https://rdrr.io/r/stats/glm.html)) on cells with
finite values and report the intercept, `S50 = -intercept / slope`
(progress at 50% usage), slope,
`timing_heterogeneity = 1 / (1 + |slope|)` and the deviance-explained
quality. Flags: `always_0` / `always_1`, `no_observations` (no finite
usage), `fit_failed`, `non_converged`, `flat_slope`, `decreasing`,
`S50_out_of_range`, `separation` (\|coefficient\| \> 30) and
`non_monotonic` (a drop larger than `nonmonotonic_tolerance` in
`usage_probability`).

## Usage

``` r
fit_timing_models(
  usage,
  progress,
  usage_probability = NULL,
  nonmonotonic_tolerance = 0.1,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- usage:

  Cells x IZs usage matrix.

- progress:

  Replication progress of the cells.

- usage_probability:

  `NULL` or bins x IZs probabilities
  (`compute_usage_heterogeneity()$probability`).

- nonmonotonic_tolerance:

  Tolerated drop between consecutive bins.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  used over IZ chunks.

## Value

A `data.frame`, one row per IZ.
