# Firing-order consistency of initiation zones from replication states

IZs are ordered by decreasing `rt` (earliest first). For each cell with
`rep_percent` in `rep_range` (inclusive), the state of the bin
containing each IZ centre is read (cells with less than 90% non-`NA` IZ
states are set to `NA`); if `K` IZs are fired (state `4`), the first `K`
IZs in timing order are expected to be fired. Observed (`4` fired, `2`
unfired) and expected states are cross-tabulated per IZ.

## Usage

``` r
calculate_firing_order(
  izs,
  genobins_rep_states,
  meta_data,
  rep_range = c(0, 1)
)
```

## Arguments

- izs:

  `GRanges` of IZs with a numeric `rt` column (no `NA`).

- genobins_rep_states:

  Matrix of states, bins x cells, with row names `"chr-start-end"`;
  every IZ centre must lie in exactly one bin.

- meta_data:

  A `data.frame` with columns `cell_id` (covering all matrix columns)
  and `rep_percent` in `[0, 1]`.

- rep_range:

  Inclusive `rep_percent` range of the cells used.

## Value

A numeric matrix, one row per IZ in timing order (row names
`"chr:start-end"`), with the percentages
`observed_fired_expected_fired`, `observed_unfired_expected_unfired`,
`observed_fired_expected_unfired` and `observed_unfired_expected_fired`
(`NA` if an IZ has no observation).
