# Find the coverage shrink factor that matches two conditions

Search for the factor `s` in `interval` that minimises the two-sample
Kolmogorov-Smirnov statistic between the coverage of the `target` cells
and `s` times the coverage of the `reference` cells, using
[`stats::optimize()`](https://rdrr.io/r/stats/optimize.html). The
returned `shrink_factor` is `round(optimum, digits) + offset` and is
meant to be used as the keep probability in
[`downsample_fragments()`](https://midimagi.github.io/scRTtools/reference/downsample_fragments.md)
for the reference samples.

## Usage

``` r
optimize_shrink_factor(
  meta_data,
  group_col = "state",
  target = "ama",
  reference = "ctr",
  coverage_col = "coverage",
  interval = c(0.1, 1),
  digits = 2L,
  offset = 0,
  tol = .Machine$double.eps^0.25
)
```

## Arguments

- meta_data:

  A `data.frame` of cell metadata.

- group_col:

  Column holding the condition labels.

- target:

  Label of the condition whose coverage is kept (e.g. `"ama"`).

- reference:

  Label of the condition whose coverage is scaled (e.g. `"ctr"`).

- coverage_col:

  Numeric coverage column (e.g. from
  [`calculate_mapd()`](https://midimagi.github.io/scRTtools/reference/calculate_mapd.md)).

- interval:

  Search interval for the factor, passed to
  [`stats::optimize()`](https://rdrr.io/r/stats/optimize.html).

- digits:

  Rounding digits applied to the optimum.

- offset:

  Number added after rounding; the default `0` returns the rounded
  optimum.

- tol:

  Tolerance passed to
  [`stats::optimize()`](https://rdrr.io/r/stats/optimize.html) (its
  default).

## Value

A list with `optimum` (unrounded minimiser), `objective` (KS statistic
at the optimum) and `shrink_factor` (`round(optimum, digits) + offset`).

## Examples

``` r
set.seed(1)
md <- data.frame(state = rep(c("ama", "ctr"), each = 300),
                 coverage = c(rlnorm(300, log(80), 0.3), rlnorm(300, log(240), 0.3)))
optimize_shrink_factor(md)
#> $optimum
#> [1] 0.3376999
#> 
#> $objective
#> [1] 0.04333333
#> 
#> $shrink_factor
#> [1] 0.34
#> 
```
