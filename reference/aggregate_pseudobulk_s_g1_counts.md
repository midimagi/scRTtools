# Aggregate per-cell CP10K counts of S and G1 cells over several samples

Reads one raw genome-bin count matrix per sample, normalises every
selected cell to counts per 10,000 (the denominator is the cell's total
over all raw bins, before any mappability filter) and sums the
normalised counts per bin over the QC-passed S cells (`is_G1` is `FALSE`
and `pass_qc` is `TRUE`) and over the G1 cells (`is_G1` is `TRUE`). The
matrices are read one at a time.

## Usage

``` r
aggregate_pseudobulk_s_g1_counts(
  genobins_files,
  meta_data,
  bin_ids,
  cache_file = NULL,
  verbose = FALSE
)
```

## Arguments

- genobins_files:

  Named character vector of RDS files, one sparse bins x cells count
  matrix per sample; the names are the sample ids.

- meta_data:

  A `data.frame` with cell ids as row names and logical columns `is_G1`
  and `pass_qc`.

- bin_ids:

  Character vector of bin ids. Each matrix must have exactly these row
  names (rows are reordered if needed).

- cache_file:

  `NULL` or the path of an RDS cache. The cache is reused only when the
  files (path, size, modification time), the selected cells and
  `bin_ids` are unchanged; otherwise it is overwritten.

- verbose:

  Print one message per sample. Default `FALSE`.

## Value

A `data.frame` with columns `bin_id`, `s_counts` and `g1_counts` (summed
CP10K, not raw reads), and attributes `normalization`
(`"cp10k_per_cell_v1"`) and `selected_cells` (cells used per sample).

## Details

Cell names of a matrix that do not start with `"<sample>-"` get that
prefix before they are matched to `rownames(meta_data)`; every cell of
every matrix must be present in `meta_data`. Each sample must contribute
at least one S cell and one G1 cell.

## See also

[`calculate_pseudobulk_rt_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/calculate_pseudobulk_rt_g1_ratio.md),
[`export_pseudobulk_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/export_pseudobulk_g1_ratio.md)
