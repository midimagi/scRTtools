# Usage heterogeneity of initiation zones

Per progress bin, the usage probability `p` of each IZ is the mean usage
of the cells in the bin (ignoring `NA`). The heterogeneity is the
cell-count-weighted mean of `4 p (1 - p)` over the bins with a finite
`p` (weights renormalised over those bins); without `NA` all bins are
used.

## Usage

``` r
compute_usage_heterogeneity(usage, progress, n_progress_bins = 10L)
```

## Arguments

- usage:

  Cells x IZs usage matrix (0/1/`NA`; base or sparse), e.g.
  [`call_iz_usage()`](https://midimagi.github.io/scRTtools/reference/call_iz_usage.md)
  or
  [`call_iz_usage_from_states()`](https://midimagi.github.io/scRTtools/reference/call_iz_usage_from_states.md).

- progress:

  Replication progress of the cells.

- n_progress_bins:

  Number of progress bins.

## Value

A list with `summary` (`IZ_id`, `mean_usage`, `usage_by_progress`,
`usage_heterogeneity`, `n_usage_observed`), `by_progress`, `probability`
(bins x IZs) and `progress_bin`.
