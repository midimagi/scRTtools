# Initiation-position and fork-asymmetry heterogeneity per IZ

Per IZ, remove the replication-progress trend from the relative
initiation position and the fork asymmetry of its tracks
([`correct_signal_matrix()`](https://midimagi.github.io/scRTtools/reference/correct_signal_matrix.md);
plain centring with fewer than 4 tracks) and summarise means, medians
and the spread (SD, MAD, IQR, variance) of the residuals.

## Usage

``` r
summarize_track_features(
  track_intervals,
  method = c("spline", "bin"),
  n_progress_bins = 10L,
  spline_df = 4L,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- track_intervals:

  `prepare_iz_data()$track_intervals`.

- method, n_progress_bins, spline_df:

  Correction settings.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  used over IZs.

## Value

A list with `summary` (one row per IZ), `track_intervals` (with
`position_residual` and `fork_asymmetry_residual`) and
`correction_method`.
