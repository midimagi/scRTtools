# Co-usage of nearby initiation zones

For every pair of IZs on the same chromosome with centres at most
`max_distance` apart, correlate their usage across cells after removing
the replication-progress trend
([`correct_signal_matrix()`](https://midimagi.github.io/scRTtools/reference/correct_signal_matrix.md));
pairs with fewer than `min_complete` complete cells or constant
residuals get `NA`. P-values (t-test of the correlation) are
BH-adjusted; significant pairs are `"co_usage"` (positive) or
`"competition"` (negative), others `"neutral"`.

## Usage

``` r
compute_nearby_iz_associations(
  usage,
  progress,
  iz_metadata,
  max_distance = 1e+06,
  correction = c("bin", "spline"),
  n_progress_bins = 10L,
  spline_df = 4L,
  min_complete = 30L,
  fdr_threshold = 0.05,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- usage:

  Cells x IZs usage matrix with IZ ids as column names.

- progress:

  Replication progress of the cells.

- iz_metadata:

  Table with `IZ_id`, `chr` and `center`.

- max_distance:

  Maximum centre distance (bp).

- correction, n_progress_bins, spline_df:

  Trend-correction settings.

- min_complete:

  Minimum number of complete cells.

- fdr_threshold:

  FDR threshold.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  used over chromosomes.

## Value

A `data.frame` of pairs: `chr`, `IZ1`, `IZ2`, `distance`, `association`,
`p_value`, `n_cells`, `FDR`, `association_type` and `distance_bin`.
