# Aggregate per-cell CP10K counts of S and G1 cells over several samples

Takes one raw genome-bin count matrix per sample, normalises every
selected cell to counts per 10,000 (the denominator is the cell's total
over all raw bins, before any mappability filter) and sums the
normalised counts per bin over the QC-passed S cells (`is_G1` is `FALSE`
and `pass_qc` is `TRUE`) and over the G1 cells (`is_G1` is `TRUE`).

## Usage

``` r
aggregate_pseudobulk_s_g1_counts(genobins, meta_data, bin_ids, verbose = FALSE)
```

## Arguments

- genobins:

  A sparse bins x cells count matrix, or a list of them (one per sample;
  the names, if any, are the sample ids).

- meta_data:

  A `data.frame` with cell ids as row names and logical columns `is_G1`
  and `pass_qc`.

- bin_ids:

  Character vector of bin ids. Each matrix must have exactly these row
  names (rows are reordered if needed).

- verbose:

  Print one message per sample. Default `FALSE`.

## Value

A `data.frame` with columns `bin_id`, `s_counts` and `g1_counts` (summed
CP10K, not raw reads), and attributes `normalization`
(`"cp10k_per_cell_v1"`) and `selected_cells` (cells used per sample).

## Details

When `genobins` is a named list, cell names of a matrix that do not
start with `"<sample>-"` get that prefix before they are matched to
`rownames(meta_data)`; a single matrix (or an unnamed list) is matched
as is. Every cell of every matrix must be present in `meta_data`, and
each matrix must contribute at least one S cell and one G1 cell.

To start from the per-sample RDS files, read them first, e.g.
`lapply(files, readRDS)` with a named `files` vector.

## See also

[`calculate_pseudobulk_rt_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/calculate_pseudobulk_rt_g1_ratio.md),
[`export_pseudobulk_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/export_pseudobulk_g1_ratio.md)
