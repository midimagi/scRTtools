# Per-cell median track length

Summarise a track table per cell (and optional extra grouping columns):
`median_length = median(width_bp, na.rm = TRUE)`,
`log_length = log(median_length)` and `rep_percent` = first
`rep_percent` of the cell times 100.

## Usage

``` r
summarize_track_lengths(track_lengths, by = character())
```

## Arguments

- track_lengths:

  A `data.frame` with columns `cell_id`, `width_bp` and `rep_percent`
  (fraction), e.g. from
  [`calc_rep_track_lengths()`](https://midimagi.github.io/scRTtools/reference/calc_rep_track_lengths.md).

- by:

  Extra grouping columns (e.g. `"state"`).

## Value

A `data.frame` with the grouping columns, `median_length`, `log_length`
and `rep_percent` (in percent), ordered by the grouping columns (C
locale).
