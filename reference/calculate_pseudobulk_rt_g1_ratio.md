# Pseudobulk replication timing as a z-scored S/G1 ratio

On the bins that are mappable and have a finite, non-negative S signal
and a finite, positive G1 signal, computes
`(s / sum(s)) / (g1 / sum(g1))` and z-scores it. Unlike
[`calculate_pseudobulk_rt()`](https://midimagi.github.io/scRTtools/reference/calculate_pseudobulk_rt.md)
there is no GC / mappability model, no removal of extreme bins and no
smoothing: the G1 pseudobulk is the only baseline.

## Usage

``` r
calculate_pseudobulk_rt_g1_ratio(counts, genome_granges)
```

## Arguments

- counts:

  A `data.frame` with columns `bin_id`, `s_counts` and `g1_counts`, as
  returned by
  [`aggregate_pseudobulk_s_g1_counts()`](https://midimagi.github.io/scRTtools/reference/aggregate_pseudobulk_s_g1_counts.md).
  `bin_id` must be identical to `names(genome_granges)`.

- genome_granges:

  Named `GRanges` with the logical metadata column `is_mappable`.

## Value

A sorted `GRanges` of the bins used, with metadata columns `score`
(z-score of the ratio), `s_g1_ratio`, `s_counts` and `g1_counts`.
