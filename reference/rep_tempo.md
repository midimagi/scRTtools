# Regional replication percentage along genome-wide replication progress

For every target region, averages the replication state of the
overlapping bins within each cell and then averages cells within groups
ordered by genome-wide replication progress. Bins contribute equally
whenever they overlap a target (strand is ignored). Missing bins are
excluded within each cell; cells without an observed bin are excluded
from their group's mean.

## Usage

``` r
rep_tempo(
  rep_matrix,
  granges,
  ncells = NULL,
  bin_granges = NULL,
  cell_rep_percent = NULL,
  meta_data = NULL,
  var_use = "rep_percent",
  nbin_rep = NULL,
  unrep_state = 2,
  rep_state = 4
)
```

## Arguments

- rep_matrix:

  Numeric matrix or `Matrix`, bins by cells, holding `unrep_state`,
  `rep_state` and `NA` only. Row and column names must be unique bin and
  cell IDs.

- granges:

  `GRanges` of target regions. IDs are taken from the `range_id` column,
  `names(granges)`, the `iz_id` column, or the coordinates, in that
  order.

- ncells:

  Optional positive integer: consecutive cells per non-overlapping group
  (cells must then be sorted by progress). Mutually exclusive with
  `nbin_rep`. If neither is given each cell is its own group.

- bin_granges:

  Optional `GRanges` with a `bin_id` column or names matching the matrix
  row names. Otherwise the row names are parsed as `chr-start-end` or
  `chr:start-end`.

- cell_rep_percent:

  Optional numeric genome-wide progress in 0-100, named by cell ID or in
  matrix column order. Defaults to 100 times the column means of the
  recoded matrix. Do not combine with `meta_data`.

- meta_data:

  Optional data.frame whose row names are cell IDs; its `var_use` column
  (finite, in 0-1) defines genome-wide progress.

- var_use:

  Column of `meta_data` holding progress.

- nbin_rep:

  Optional positive integer: number of equal-width progress intervals
  spanning 0-1. Cells in an interval are merged with equal weight; empty
  intervals stay as `NA` columns. Requires `meta_data`. Genome-wide
  analyses typically use `nbin_rep = 100`.

- unrep_state, rep_state:

  Codes of the unreplicated and replicated states in `rep_matrix`. The
  defaults (2 and 4) are the codes written by
  [`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md).
  Use 0 and 1 for a matrix that is already binary.

## Value

A data.frame (regions by progress groups) of percentages in 0-100.
Attributes: `rep_bins`, `range_info`, `cell_info`, `n_valid_cells` and
`units`. Save it with [`saveRDS()`](https://rdrr.io/r/base/readRDS.html)
to keep the attributes, and subset it with
[`subset_rep_tempo()`](https://midimagi.github.io/scRTtools/reference/subset_rep_tempo.md).

## See also

[`fit_rep_tempo()`](https://midimagi.github.io/scRTtools/reference/fit_rep_tempo.md)
