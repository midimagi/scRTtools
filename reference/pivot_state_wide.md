# Long table to one column per state

Spreads `feature` of a long table (one row per `id_col` and `state`)
into one column per state.

## Usage

``` r
pivot_state_wide(
  df,
  id_col,
  feature,
  na = c("drop", "zero", "keep"),
  states = c("ctr", "ama")
)
```

## Arguments

- df:

  data.frame with the `id_col` column(s), a `state` column and
  `feature`.

- id_col:

  Name(s) of the identifier column(s).

- feature:

  Name of the value column.

- na:

  `"drop"`: remove rows with any `NA`; `"zero"`: remove rows where both
  states are `NA` and set the remaining `NA` to 0; `"keep"`: leave as
  is.

- states:

  The two state names, `c(reference, treatment)`.

## Value

A data.frame with `id_col` and one column per state found in `df$state`
(in order of first appearance), one row per identifier.

## Examples

``` r
long <- data.frame(id = rep(c("a", "b", "c"), 2),
                   state = rep(c("ctr", "ama"), each = 3),
                   usage = c(1, 2, NA, 2, 2, 5))
pivot_state_wide(long, "id", "usage")
#>   id ctr ama
#> 1  a   1   2
#> 2  b   2   2
pivot_state_wide(long, "id", "usage", na = "zero")
#>   id ctr ama
#> 1  a   1   2
#> 2  b   2   2
#> 3  c   0   5
```
