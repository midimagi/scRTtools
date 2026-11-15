# Find local maxima in a binned genomic signal

Within each chromosome, bins are split into blocks at gaps between
consecutive starts larger than `max_gap` (default: 1.5 times the median
bin step). Runs of equal scores (within `plateau_tol`) form plateaus. A
plateau of at most `max_plateau_bins` bins is a peak if it has at least
`max(1, ceiling(flank_bins / 2))` non-`NA` bins among the `flank_bins`
bins on each side, is higher than all of them by more than
`plateau_tol`, is at least `min_score` (unless `NULL`) and exceeds the
median of the flanking bins by at least `min_delta`. The middle bin of
the plateau is reported.

## Usage

``` r
find_rt_peaks(
  granges,
  score_col = "score",
  flank_bins = 5,
  min_score = 0.2,
  min_delta = 0.02,
  max_plateau_bins = 5,
  max_gap = NULL,
  plateau_tol = 1e-08
)
```

## Arguments

- granges:

  `GRanges` of bins with a numeric score column (e.g. from
  [`calculate_pseudobulk_rt()`](https://midimagi.github.io/scRTtools/reference/calculate_pseudobulk_rt.md)).

- score_col:

  Name of the score metadata column.

- flank_bins:

  Positive integer; flank size on each side.

- min_score:

  `NULL` or the minimum peak score. We use `0.2` on RT scores and `10`
  on track-centre counts.

- min_delta:

  Minimum difference between the peak and the flank median.

- max_plateau_bins:

  Maximum plateau size.

- max_gap:

  `NULL` or the maximum start-to-start distance within a block.

- plateau_tol:

  Tolerance for equal scores.

## Value

A sorted `GRanges` of peak bins with additional metadata columns
`RT_value` (the score), `local_median`, `local_delta` and
`plateau_size`.
