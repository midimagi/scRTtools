# Label cells or regions as early / middle / late by replication progress

Label cells or regions as early / middle / late by replication progress

## Usage

``` r
add_rep_cluster(df, breaks)
```

## Arguments

- df:

  data.frame with a `rep_percent` column.

- breaks:

  Four increasing boundaries in the unit of `rep_percent`, e.g.
  `c(0, 33, 66, 100)` or `c(0, 0.33, 0.66, 1)`. Intervals are closed; a
  value on a boundary goes to the earlier class.

## Value

`df` with a factor column `rep_cluster` (levels `early`, `middle`,
`late`; `NA` outside the breaks).

## Examples

``` r
add_rep_cluster(data.frame(rep_percent = c(5, 33, 50, 99, 120)),
                breaks = c(0, 33, 66, 100))
#>   rep_percent rep_cluster
#> 1           5       early
#> 2          33       early
#> 3          50      middle
#> 4          99        late
#> 5         120        <NA>
```
