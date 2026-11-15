# NA-tolerant summaries for heatmap tables

Helpers used to build the iZone velocity and gene tempo heatmaps. All
three drop non-finite values before computing.

## Usage

``` r
cor_safe(x, y)

zscore_safe(x)

peak_bin(position, value)
```

## Arguments

- x, y:

  Numeric vectors of the same length.

- position:

  Numeric vector of positions (e.g. replication progress bins).

- value:

  Numeric vector of values at `position`.

## Value

`cor_safe()` and `peak_bin()`: a single number. `zscore_safe()`: a
numeric vector as long as `x`.

## Details

- `cor_safe()` returns the Pearson correlation of the finite pairs of
  `x` and `y`, or `NA_real_` when fewer than 2 pairs are left.

- `zscore_safe()` centres and scales the finite values of `x`.
  Non-finite entries stay `NA`. The whole result is `NA` when fewer than
  2 values are finite or their standard deviation is 0 or not finite.

- `peak_bin()` returns the `position` at which `value` is largest, using
  only pairs where both are finite (first one on ties), or `NA_real_`
  when there are none.

## Examples

``` r
cor_safe(c(1, 2, NA, 4), c(2, 4, 6, Inf))
#> [1] 1
zscore_safe(c(1, 2, 3, NA))
#> [1] -1  0  1 NA
peak_bin(c(0, 25, 50, 75), c(0.1, 0.8, NaN, 0.3))
#> [1] 25
```
