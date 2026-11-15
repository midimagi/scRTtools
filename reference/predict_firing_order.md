# Firing-order consistency of initiation zones from track centres

Like
[`calculate_firing_order()`](https://midimagi.github.io/scRTtools/reference/calculate_firing_order.md),
but an IZ is fired in a cell when at least one of the cell's track
centres (`floor((start + end) / 2)`) falls inside it. Cells are selected
by `rep_range[1] <= rep_percent < rep_range[2]` and tracks optionally
filtered with `track_filter`.

## Usage

``` r
predict_firing_order(
  gr_izones,
  rep_tracks,
  rep_range = c(0, 1),
  track_filter = NULL
)
```

## Arguments

- gr_izones:

  `GRanges` of IZs with a numeric `rt` column (no `NA`).

- rep_tracks:

  Named `GRangesList` of tracks with metadata columns `cell_id` (equal
  to the element name), `rep_percent` and `num_bins`.

- rep_range:

  `rep_percent` range (left-closed, right-open).

- track_filter:

  Optional unquoted expression evaluated against each cell's track
  metadata (non-standard evaluation), e.g. `num_bins <= 10`.

## Value

A numeric matrix as in
[`calculate_firing_order()`](https://midimagi.github.io/scRTtools/reference/calculate_firing_order.md).
