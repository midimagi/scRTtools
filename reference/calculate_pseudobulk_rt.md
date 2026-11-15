# Calculate a pseudobulk replication-timing signal

Aggregate QC-passed S-phase cells (`pass_qc & !is_G1`) and G1 cells
(`is_G1`), correct the S-phase counts for GC content and mappability
using the G1 pseudobulk, normalise by the genome-wide mean, z-score, set
bins with `|z| > extreme_cutoff` to `NA`, smooth with
[`stats::smooth.spline()`](https://rdrr.io/r/stats/smooth.spline.html)
within mappable contigs, then z-score and filter the smoothed signal
again.

## Usage

``` r
calculate_pseudobulk_rt(
  genobins_matrix,
  genome_granges,
  meta_data,
  extreme_cutoff = 5,
  min.gapwidth = 40001L,
  min_bins_per_contig = 5,
  min_non_na_points = 2,
  smooth_spar = NULL
)
```

## Arguments

- genobins_matrix:

  `dgCMatrix` of raw counts, genome bins x cells, with row names equal
  to `names(genome_granges)` and column names equal to
  `rownames(meta_data)`.

- genome_granges:

  Named `GRanges` with metadata columns `is_mappable`, `gc_frequency`
  and `mappability` (e.g. from
  [`compute_bin_mappability()`](https://midimagi.github.io/scRTtools/reference/compute_bin_mappability.md)).

- meta_data:

  A `data.frame` with logical columns `is_G1` and `pass_qc`.

- extreme_cutoff:

  Absolute z-score above which bins are set to `NA`.

- min.gapwidth:

  Passed to
  [`GenomicRanges::reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html)
  to build contigs.

- min_bins_per_contig:

  Minimum number of bins of a smoothed contig.

- min_non_na_points:

  Minimum number of non-`NA` points for the spline.

- smooth_spar:

  `NULL` (automatic) or `spar` of
  [`stats::smooth.spline()`](https://rdrr.io/r/stats/smooth.spline.html).

## Value

A sorted `GRanges` of the bins with a non-`NA` smoothed signal, with the
single metadata column `score`.

## Details

GC/mappability correction (on mappable bins): with
`g1m = g1_counts / mappability`, `mGC` = median of `g1m` within bins of
the same `round(gc_frequency, 2)` and `m` = median of all `g1m`, the
corrected S signal is `(m / mGC) * s_counts / mappability` (`NA` where
`mGC == 0`). Contigs are mappable bins merged with
`GenomicRanges::reduce(min.gapwidth = min.gapwidth)`; contigs with fewer
than `min_bins_per_contig` bins or fewer than `min_non_na_points`
distinct positions are not smoothed.
