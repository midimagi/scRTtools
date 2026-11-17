# S/G1 ratio pseudobulk for a group of samples, optionally written to bigWig

Runs
[`aggregate_pseudobulk_s_g1_counts()`](https://midimagi.github.io/scRTtools/reference/aggregate_pseudobulk_s_g1_counts.md)
and
[`calculate_pseudobulk_rt_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/calculate_pseudobulk_rt_g1_ratio.md)
for one group of samples (for example the replicates of one condition).
The result is returned, and a bigWig is written only when `file` is
given.

## Usage

``` r
export_pseudobulk_g1_ratio(
  genobins,
  meta_data,
  genome_granges,
  file = NULL,
  ...
)
```

## Arguments

- genobins:

  A sparse bins x cells count matrix, or a list of them (one per sample;
  the names, if any, are the sample ids).

- meta_data:

  A `data.frame` with cell ids as row names and logical columns `is_G1`
  and `pass_qc`.

- genome_granges:

  Named `GRanges` with the metadata column `is_mappable`. Its names are
  the bin ids. Writing a bigWig needs chromosome lengths in its
  `seqinfo`.

- file:

  `NULL` (default, nothing is written) or the path of the bigWig that
  receives the `score` column (needs the rtracklayer package).

- ...:

  Passed to
  [`aggregate_pseudobulk_s_g1_counts()`](https://midimagi.github.io/scRTtools/reference/aggregate_pseudobulk_s_g1_counts.md),
  e.g. `verbose`.

## Value

The `GRanges` of
[`calculate_pseudobulk_rt_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/calculate_pseudobulk_rt_g1_ratio.md),
with the per-sample cell numbers in
`S4Vectors::metadata(x)$selected_cells`; returned invisibly when `file`
is given.
