# Track-shape heterogeneity per IZ

For IZs with at least `min_active_cells` tracks, the binary coverage
profile of each track around the IZ centre is compared with its
progress-expected profile
([`correct_signal_matrix()`](https://midimagi.github.io/scRTtools/reference/correct_signal_matrix.md)):
mean correlation (`mean_track_similarity`; `1 - RMSE` where undefined),
residual variance and the residual error after a rank-`pca_rank` SVD.

## Usage

``` r
summarize_track_shape(
  track_intervals,
  iz_metadata,
  method = c("spline", "bin"),
  n_progress_bins = 10L,
  spline_df = 4L,
  window_bp = 500000L,
  bin_size = 20000L,
  min_active_cells = 5L,
  pca_rank = 2L,
  compute_pca = TRUE,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- track_intervals, iz_metadata:

  From
  [`prepare_iz_data()`](https://midimagi.github.io/scRTtools/reference/prepare_iz_data.md).

- method, n_progress_bins, spline_df:

  Correction settings.

- window_bp, bin_size:

  Profile window half-width and step (bp).

- min_active_cells:

  Minimum number of tracks.

- pca_rank, compute_pca:

  SVD settings.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  used over IZs.

## Value

A `data.frame`, one row per IZ.
