# Define G1-derived windows from a fragment file

Build variable-width genomic windows that each hold `target_reads` G1
fragments, and count fragments of all cells in those windows.

## Usage

``` r
define_G1_windows(
  genome_granges,
  frag_path,
  meta_data,
  g1_col = "is_G1",
  cov_col = "coverage",
  target_reads = 200,
  target_window_length = 20000,
  min.gapwidth = 100000L,
  verbose = FALSE
)
```

## Arguments

- genome_granges:

  A `GRanges` of retained genome bins (already filtered for gaps / CNV /
  mappability, e.g. the output of
  [`filter_windows()`](https://midimagi.github.io/scRTtools/reference/filter_windows.md),
  unlisted, or of
  [`synchronize_bins()`](https://midimagi.github.io/scRTtools/reference/synchronize_bins.md)).

- frag_path:

  Path to a bgzipped, tabix-indexed 10x fragment file (`<frag_path>.tbi`
  must exist). Columns: chr, start (0-based), end, barcode, ...

- meta_data:

  A `data.frame` whose unique row names are cell barcodes.

- g1_col:

  Column of `meta_data` flagging G1 cells (coerced to logical).

- cov_col:

  Numeric column of `meta_data` with per-cell coverage.

- target_reads:

  Positive integer; G1 fragments per window. We use `200`.

- target_window_length:

  Positive integer; target median window length (bp) and minimum contig
  width, e.g. `20000`.

- min.gapwidth:

  Positive integer passed to
  [`GenomicRanges::reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html).
  Default `100000L`.

- verbose:

  Logical; report progress. Default `FALSE`.

## Value

A list with

- `g1bins_counts`:

  `dgCMatrix` windows x cells (column names = `rownames(meta_data)`, row
  names = window IDs `"<chr>-<start>-<end>"`).

- `g1bins`:

  named `GRanges` of windows with metadata columns `midpoint`,
  `g1_read_count`, `is_masked_window` (G1 count differs from
  `target_reads` by more than 1) and `contig_id`. Its
  [`S4Vectors::metadata()`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html)
  holds `target_reads`, `target_window_length`, `pilot_chr`,
  `selected_g1_n`, `selected_g1_barcodes`, `selected_cov_cutoff`,
  `pilot_median_length`, `pilot_n_windows`.

## Details

Steps:

1.  `genome_granges` is merged with
    `GenomicRanges::reduce(min.gapwidth)` into contigs; contigs shorter
    than `target_window_length` and contigs on chromosomes absent from
    the tabix index are dropped.

2.  G1 cells (`meta_data[[g1_col]]` is `TRUE`, finite coverage) are
    ranked by decreasing coverage.

3.  On a pilot chromosome (largest total contig width), a binary search
    chooses the number `n` of top-ranked G1 cells whose windows have a
    median length closest to `target_window_length`.

4.  For every chromosome, fragment midpoints
    (`((start0 + 1) + end) %/% 2`) of the `n` selected G1 cells are
    sorted per contig and cut into consecutive blocks of `target_reads`;
    each window runs from the first fragment of its block to just before
    the first fragment of the next block (the last block ends at its
    last fragment; the incomplete remainder is dropped). If a chromosome
    has `>= 100` windows, windows outside the 1st-99th length
    percentiles are removed.

5.  Fragments of all cells in `rownames(meta_data)` are counted per
    window (by midpoint).
