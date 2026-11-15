# Identify replication tracks from bin-level replication states

Identify contiguous segments with replication state `cn` (typically `4`)
in each cell of a bins x cells state matrix. R translation of the MATLAB
functions `define_replication_tracks()` and `list_cn_segments()`.

## Usage

``` r
define_replication_tracks(
  genobins_matrix,
  genome_granges,
  rep_percent,
  cn = 4,
  max_na_run = 5L,
  max_na_fraction = 0.5,
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
```

## Arguments

- genobins_matrix:

  Numeric matrix of replication states (`2`, `4` or `NA`), bins x cells,
  with row names equal to `names(genome_granges)` and unique column
  names (e.g. `compute_rep_percent()$genobins_rep_states_filtered`).

- genome_granges:

  Named `GRanges`, one range per row, with an `is_mappable` metadata
  column (`NA` counts as not mappable).

- rep_percent:

  Named numeric vector of replicated fractions; names must include all
  column names of `genobins_matrix` (no `NA`).

- cn:

  State to extract.

- max_na_run:

  Non-negative integer; a new segment starts after more than this many
  consecutive `NA` bins.

- max_na_fraction:

  Number in `[0, 1]`; maximum `NA` fraction of a kept segment.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html);
  chromosomes are processed in parallel.

- verbose:

  Logical; print progress messages. Default `FALSE`.

## Value

A `GRangesList` named by cell; each element holds the cell's segments
with metadata columns `cell_id`, `rep_percent`, `cn`, `num_bins`, `flag`
and `na_fraction`.

## Details

Within each chromosome and cell, bins that are `NA` in all cells or not
mappable are treated as `NA`. A new segment starts at the first bin,
when the state changes, at the first non-`NA` bin after an `NA` bin with
no previous state in the segment, and after more than `max_na_run`
consecutive `NA` bins. Trailing `NA` bins belong to the preceding
segment. Segments whose `NA` fraction exceeds `max_na_fraction` are
dropped, and only segments with state `cn` are returned.

The `flag` column (MATLAB fourth column) is `1` if a segment neighbours
a segment without a state (all `NA`) or is adjacent to another segment
with state `cn`.
