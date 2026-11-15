# Prepare track, cell and IZ tables for heterogeneity analysis

Collect the IZ-assigned tracks (one row per cell and IZ; duplicated
cell-IZ assignments are merged: outer start/end, width-weighted centre,
maximum flag), the IZ metadata and the cell replication progress, and
derive per-track features relative to the IZ centre (`q50`): relative
initiation position, left/right/total track length and fork asymmetry.

## Usage

``` r
prepare_iz_data(
  tracks,
  gr_izones,
  cell_metadata = NULL,
  progress_col = "rep_percent"
)
```

## Arguments

- tracks:

  Named `GRangesList` of tracks (every cell with at least one track)
  with columns `cell_id`, `rep_percent`, `iz_id` and `flag` (e.g.
  `define_initiation_zones()$tracks`).

- gr_izones:

  `GRanges` of IZs with `iz_id`, `q5`-`q95`, `n_tracks` and optionally
  the firing range (e.g.
  [`calculate_firing_range()`](https://midimagi.github.io/scRTtools/reference/calculate_firing_range.md)).

- cell_metadata:

  `NULL` or a `data.frame` with `cell_id` and `progress_col`; if it has
  `NA`s, the track `rep_percent` is used.

- progress_col:

  Progress column of `cell_metadata`.

## Value

A list with `track_intervals`, `cell_metadata` (with
`replication_progress`), `iz_metadata` and `diagnostics`.
