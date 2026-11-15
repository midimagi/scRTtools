# Adjacent IZ pairs with firing-range features

Order IZs by chromosome (seqlevel order), start and end, attach the
firing range of each IZ (matched by `iz_id`) and build all adjacent
pairs on the same chromosome with distances and firing-range overlaps.

## Usage

``` r
make_adjacent_pairs(izones, fire_ranges)
```

## Arguments

- izones:

  `GRanges` of IZs with an `iz_id` column.

- fire_ranges:

  Output of
  [`calculate_firing_range()`](https://midimagi.github.io/scRTtools/reference/calculate_firing_range.md)
  (`iz_id`, `earliest_rep_percent`, `latest_unrep_percent`, optionally
  `corroborated`), with the same coordinates for matching ids.

## Value

A list with `izones` (ordered), `features` (one row per IZ) and `pairs`
(one row per adjacent pair: indices into `features`, ids, coordinates,
`center_distance_bp`, `gap_bp`, `overlap_bp`, left/right firing-range
columns, `firing_range_overlap`, `earliest_rep_difference`).

## Details

A firing range is valid when both bounds are finite, `earliest >= 0`,
`latest <= 1` and `earliest <= latest`; `firing_range_overlap` is
computed for pairs with two valid ranges.
