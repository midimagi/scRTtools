# Crossing point of the velocity curves of two states

Interpolates both curves on a common grid of 1000 points over the shared
`rep_percent` range and returns the first sign change of their
difference.

## Usage

``` r
find_intersection(data, state1, state2)
```

## Arguments

- data:

  data.frame with columns `state`, `rep_percent` and `length_velocity`
  (bp), e.g. the prediction table of
  [`fit_track_length_models()`](https://midimagi.github.io/scRTtools/reference/fit_track_length_models.md).

- state1, state2:

  The two states to compare.

## Value

A one-row data.frame with `x` (`rep_percent` at the crossing) and `y`
(velocity there, in kb), or `NULL` when the curves do not cross.

## Examples

``` r
d <- data.frame(state = rep(c("ctr", "ama"), each = 11),
                rep_percent = rep(seq(0, 100, 10), 2),
                length_velocity = c(seq(1e4, 3e4, length.out = 11),
                                    seq(2.5e4, 1.5e4, length.out = 11)))
find_intersection(d, "ctr", "ama")
#>    x  y
#> 1 50 20
```
