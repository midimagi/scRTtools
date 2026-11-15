# Co-firing of IZ pairs in replication-progress bins

An IZ is fired in a cell when one of the cell's tracks is assigned to it
(`iz_id`). Cells are grouped into `bin_num` equal bins of `rep_percent`
(cells without `rep_percent` form one `NA` bin); per pair and bin the
numbers of cells with the left / right IZ fired are summed, and per pair
the bins are classified by which members fired at least once.

## Usage

``` r
summarise_pair_cofiring_bins(pairs, tracks, meta_data, bin_num = 100L)
```

## Arguments

- pairs:

  Output of
  [`make_pair_init()`](https://midimagi.github.io/scRTtools/reference/make_pair_init.md).

- tracks:

  Named `GRangesList` of tracks with an `iz_id` column (e.g.
  `define_initiation_zones()$tracks`).

- meta_data:

  A `data.frame` with columns `cell_id` and `rep_percent`.

- bin_num:

  Number of `rep_percent` bins.

## Value

A list with `pair_bins` (per pair and bin: `bin_mid`, median
`rep_percent`, summed `left_iz_state` / `right_iz_state`) and
`pair_fired_info` (per pair: numbers of bins with neither / only left /
only right / both fired, `n_all` and the corresponding rates).
