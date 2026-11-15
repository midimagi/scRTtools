# Combine the estimates of two samples

Fixed-effect (inverse-variance) combination of per-sample estimates, one
pair per combination of the `by` columns, with a normal p-value,
Cochran's Q heterogeneity test (1 df) and `I2`. `replication_p` is the
larger of the two per-sample p-values when both estimates have the same
non-zero sign and 1 otherwise: evidence that both samples support the
effect.

## Usage

``` r
combine_two_samples(tab, by, sample_col = "sample", failed = c("keep", "drop"))
```

## Arguments

- tab:

  data.frame of per-sample fits with columns `beta`, `se`, `p`, the `by`
  columns, the sample column and optionally `status`, e.g. the stacked
  output of
  [`fit_ols_hc3()`](https://midimagi.github.io/scRTtools/reference/fit_ols_hc3.md).
  Exactly two samples, each once per combination of `by`.

- by:

  Columns identifying one test.

- sample_col:

  Column holding the sample name.

- failed:

  `"keep"`: failed pairs stay in the result with `beta`/`se` `NA`,
  `p = 1`, `replication_p = 1` and `status = "fit_failed"`, so they
  count in later multiple-testing corrections; `"drop"` removes them.

## Value

A data.frame with one row per test (in order of first appearance): the
`by` columns, `beta`, `se`, `p`, `replication_p`, `same_direction`, `Q`,
`heterogeneity_p`, `I2`, `status`, and `beta_<sample>`, `se_<sample>`,
`p_<sample>` for both samples (sorted by name).

## Details

A pair is failed when a standard error is missing or not positive or, if
`tab` has a `status` column, when a status is not `"ok"`. Combining two
samples does not add biological replicates.

## Examples

``` r
fits <- data.frame(gene = rep(c("a", "b"), each = 2),
                   sample = rep(c("s1", "s2"), 2),
                   beta = c(0.5, 0.4, 0.3, -0.2), se = c(0.1, 0.12, 0.1, 0.1),
                   p = c(1e-6, 1e-3, 3e-3, 0.05))
combine_two_samples(fits, by = "gene")
#>   gene      beta         se            p replication_p same_direction
#> 1    a 0.4590164 0.07682213 2.300134e-09         0.001           TRUE
#> 2    b 0.0500000 0.07071068 4.795001e-01         1.000          FALSE
#>            Q heterogeneity_p   I2 status beta_s1 se_s1  p_s1 beta_s2 se_s2
#> 1  0.4098361     0.522052724 0.00     ok     0.5   0.1 1e-06     0.4  0.12
#> 2 12.5000000     0.000406952 0.92     ok     0.3   0.1 3e-03    -0.2  0.10
#>    p_s2
#> 1 0.001
#> 2 0.050
```
