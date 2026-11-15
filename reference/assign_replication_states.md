# Assign replication states to G1 and genome windows

Run
[`infer_single_cell_states()`](https://midimagi.github.io/scRTtools/reference/infer_single_cell_states.md)
on every cell of a G1-window count matrix (in chunks of `chunk_size`
cells, parallelised with `BPPARAM`), using the chromosome starts of
`g1_granges` as HMM breakpoints, then project the states onto genome
windows.

## Usage

``` r
assign_replication_states(
  g1bins_matrix,
  g1_granges,
  genome_granges,
  gap_max_bp = 50000L,
  BPPARAM = BiocParallel::SerialParam(),
  chunk_size = 100L,
  infer_args = list(hmm_aggr_width = 15L, init_ratio_lower = 1.2, final_ratio_lower =
    1.4),
  verbose = FALSE
)
```

## Arguments

- g1bins_matrix:

  A `dgCMatrix` of counts, G1 windows x cells, with row names identical
  to `names(g1_granges)` and non-`NULL` column names (e.g.
  `define_G1_windows()$g1bins_counts`).

- g1_granges:

  `GRanges` of the G1 windows, grouped by chromosome, with a logical
  `is_masked_window` metadata column.

- genome_granges:

  Named `GRanges` of genome windows, grouped by chromosome, with a
  logical `is_mappable` metadata column (e.g. from
  [`compute_bin_mappability()`](https://midimagi.github.io/scRTtools/reference/compute_bin_mappability.md)).

- gap_max_bp:

  Maximum total width (bp) of an unmappable gap to fill.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  used over cell chunks.

- chunk_size:

  Number of cells per chunk.

- infer_args:

  Named list of extra arguments for
  [`infer_single_cell_states()`](https://midimagi.github.io/scRTtools/reference/infer_single_cell_states.md).
  A list given here replaces the default as a whole; tune the values
  with
  [`evaluate_state_parameters()`](https://midimagi.github.io/scRTtools/reference/evaluate_state_parameters.md).

- verbose:

  Logical; print progress messages. Default `FALSE`.

## Value

A list with

- `genobins_rep_states`:

  Integer matrix, genome windows x cells, values `2`, `4` or `NA`; input
  for
  [`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md).

- `g1bins_rep_states`:

  Integer matrix, G1 windows x cells.

## Details

States of masked G1 windows (`mcols(g1_granges)$is_masked_window`, `NA`
counts as masked) are set to `NA`. Each genome window takes the state of
the unmasked G1 window on the same chromosome with the largest overlap;
unmappable genome windows (`mcols(genome_granges)$is_mappable` `FALSE`
or `NA`) are set to `NA`, and runs of unmappable windows of at most
`gap_max_bp` bp that are not at a chromosome end are filled when both
flanking windows have the same state.
