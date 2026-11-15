# Calculate per-cell coverage and MAPD from a tile count matrix

Compute per-cell coverage, MAPD (median absolute pairwise difference)
and scaled MAPD from a tile-by-cell count matrix.

## Usage

``` r
calculate_mapd(
  tile_matrix,
  genome_granges,
  meta_data,
  aggr_bins = 50,
  verbose = FALSE,
  BPPARAM = BiocParallel::SerialParam()
)
```

## Arguments

- tile_matrix:

  A tile-by-cell count matrix (base or Matrix sparse matrix) with row
  names equal to `names(genome_granges)` and column names equal to
  `rownames(meta_data)`, in the same order.

- genome_granges:

  A named `GRanges` with one range per row of `tile_matrix` and metadata
  columns `gc_frequency`, `mappability` and `is_mappable` (see
  [`compute_bin_mappability()`](https://midimagi.github.io/scRTtools/reference/compute_bin_mappability.md)).

- meta_data:

  A `data.frame` of cell metadata, row names = cell IDs.

- aggr_bins:

  Positive integer; number of consecutive tiles aggregated before
  computing MAPD. We use `50`.

- verbose:

  Logical; print progress messages. Default `FALSE`.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html);
  cells are processed in chunks in parallel. Results do not depend on
  the number of workers.

## Value

`meta_data` with added or updated columns `coverage`, `mapd`,
`scaled_mapd` (numeric) and `mapd_error` (character, `NA` when the cell
succeeded).

## Details

Steps:

1.  A per-tile GC correction value is computed from the pseudo-bulk (row
    sums): tiles are grouped by `round(gc_frequency, 1)`, and the value
    is
    `median(all reads/mappability) / median(group reads/mappability)`.
    Only mappable tiles with finite GC, positive mappability and finite
    correction are used further.

2.  For each cell, corrected reads are `reads / mappability * gc_corr`;
    `coverage = 1e6 * sum(corrected) / effective genome size` (sum of
    the widths of the used tiles).

3.  Corrected reads are summed in groups of `aggr_bins` consecutive used
    tiles per chromosome; MAPD is the median absolute deviation (around
    the median, unscaled) of the differences between adjacent groups on
    the same chromosome, and
    `scaled_mapd = round(mapd / sqrt(coverage), 3)`.

## Examples

``` r
set.seed(1)
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, by = 1000, length.out = 200),
                                                      width = 1000))
names(gr) <- paste0("bin", seq_along(gr))
gr$gc_frequency <- runif(200, 0.35, 0.5)
gr$mappability <- 1
gr$is_mappable <- TRUE
m <- matrix(rpois(200 * 5, 5), 200, 5, dimnames = list(names(gr), paste0("c", 1:5)))
m <- Matrix::Matrix(m, sparse = TRUE)
md <- data.frame(row.names = colnames(m), cell_id = colnames(m))
calculate_mapd(m, gr, md, aggr_bins = 10, verbose = FALSE)
#>    cell_id coverage      mapd scaled_mapd mapd_error
#> c1      c1 4967.292  6.250000       0.089       <NA>
#> c2      c2 5207.500  6.500000       0.090       <NA>
#> c3      c3 5125.000  3.708333       0.052       <NA>
#> c4      c4 5006.042 10.083333       0.143       <NA>
#> c5      c5 4924.167  5.500000       0.078       <NA>
```
