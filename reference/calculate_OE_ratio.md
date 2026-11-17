# GC-correct read-count windows for one sample

Estimate per-GC-bin bias factors from the reads that fall inside the
retained windows `wins_filtered`, then compute expected counts and the
observed/expected ratio for every retained window.

## Usage

``` r
calculate_OE_ratio(
  wins_filtered,
  wins_gc,
  reads_list,
  coords_list,
  gc_bins_use = 81:321,
  verbose = FALSE
)
```

## Arguments

- wins_filtered:

  Named `GRangesList` of read-count windows with numeric metadata column
  `raw_counts` (e.g. output of
  [`filter_windows()`](https://midimagi.github.io/scRTtools/reference/filter_windows.md)).
  **Ranges must be named**: windows are matched to `wins_gc` by
  [`names()`](https://rdrr.io/r/base/names.html).

- wins_gc:

  Named `GRangesList` (same chromosome names) with named ranges and
  metadata columns `GC_bin_000` ... `GC_bin_400` (bases per GC bin) and
  `row_sum`.

- reads_list:

  Named list of numeric read-position vectors, one per chromosome (e.g.
  `count_reads()$reads_list`).

- coords_list:

  Genome positions per GC bin, either a named list (one element per
  chromosome, each a list of 401 numeric vectors) or a function taking a
  chromosome name and returning that list. The function form keeps only
  one chromosome in memory, e.g.
  `function(chr) readRDS(coords_files[[chr]])`.

- gc_bins_use:

  Integer indices (1-401) of GC bins used for correction.

- verbose:

  Logical; print progress messages. Default `FALSE`.

## Value

A list with

- `wins_ratio`:

  `wins_filtered` with added numeric metadata columns `expected_counts`
  and `oe_ratio`. `oe_ratio` is `NA` where `row_sum != win_size` or the
  expected count is 0.

- `gc_factors_per_bin`:

  named numeric vector (length 401) of GC bias factors, `NA` for GC bins
  without genome coverage.

## Details

Reads are assigned greedily to the 401 GC bins (GC 0.000-0.400 step
0.001 bins as produced by the tigeR genome processing), visiting the GC
bins with the most retained positions first. The bias factor of a GC bin
is `(reads fraction) / (genome fraction)`. The expected count of a
window is
`sum_{gc in gc_bins_use} bases_in_gc_bin * factor * median(raw_counts) / win_size / 2`,
where `win_size` is the median of `row_sum` in `wins_gc`.
