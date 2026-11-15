# Subset a `rep_tempo` table and keep its attributes

Row-subsetting a data.frame drops the custom attributes that
[`fit_rep_tempo()`](https://midimagi.github.io/scRTtools/reference/fit_rep_tempo.md)
needs; this helper carries them along.

## Usage

``` r
subset_rep_tempo(x, idx)
```

## Arguments

- x:

  Output of
  [`rep_tempo()`](https://midimagi.github.io/scRTtools/reference/rep_tempo.md).

- idx:

  Row index (integer, logical or region IDs).

## Value

The subset of `x` with `rep_bins`, `n_valid_cells` and `range_info`
subset to match.
