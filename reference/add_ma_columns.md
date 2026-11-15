# Add A / M columns comparing two states

Add A / M columns comparing two states

## Usage

``` r
add_ma_columns(
  df,
  type = c("log2", "relative", "sum"),
  states = c("ctr", "ama")
)
```

## Arguments

- df:

  Wide data.frame with one column per state, e.g. from
  [`pivot_state_wide()`](https://midimagi.github.io/scRTtools/reference/pivot_state_wide.md).

- type:

  `"log2"`: `A` = mean of the log2 values, `M` = log2 fold change
  (treatment / reference); `"relative"`: `A` = reference,
  `M = (treatment - reference) / reference`; `"sum"`: `A` = sum,
  `M = treatment - reference`.

- states:

  The two state names, `c(reference, treatment)`.

## Value

`df` with columns `A` and `M` added.

## Examples

``` r
add_ma_columns(data.frame(id = 1:3, ctr = c(1, 2, 4), ama = c(2, 2, 1)))
#>   id ctr ama   A  M
#> 1  1   1   2 0.5  1
#> 2  2   2   2 1.0  0
#> 3  3   4   1 1.0 -2
```
