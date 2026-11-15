# Compute percent replication and QC cells from bin-level replication states

Run two QC passes on a bins x cells replication-state matrix (states `2`
= unreplicated, `4` = replicated, `NA` = unknown), compute the per-cell
replicated fraction and return the filtered matrix with updated
metadata.

## Usage

``` r
compute_rep_percent(
  genobins_rep_states,
  genome_granges,
  meta_data,
  min_rep_percent = 0.001,
  max_sigma = 0.4,
  neighbor_window = 5,
  outlier_threshold = 0.2,
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
```

## Arguments

- genobins_rep_states:

  Matrix-like bins x cells replication states; row names =
  `names(genome_granges)`, column names = `rownames(meta_data)` (same
  order).

- genome_granges:

  Named `GRanges` with one range per row.

- meta_data:

  A `data.frame` whose unique row names are cell IDs.

- min_rep_percent:

  Minimum replicated fraction (first pass). We use `0.01`.

- max_sigma:

  Maximum `sigma_chr` (first pass). We use `0.3`.

- neighbor_window:

  Neighbourhood half-width for the second pass. We use `5`.

- outlier_threshold:

  Maximum fraction of local correlation outliers (second pass). We use
  `0.2`.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html);
  chromosomes are processed in parallel in the second pass.

- verbose:

  Logical; print progress messages. Default `FALSE`.

## Value

A list with

- `genobins_rep_states_filtered`:

  numeric matrix, all bins x cells passing both QC passes (columns
  ordered by increasing `rep_percent`), with second-pass masked
  chromosome blocks set to `NA`.

- `meta_data`:

  input `meta_data` with added columns `rep_percent`, `sigma_chr`
  (numeric, all cells) and `pass_qc` (logical).

## Details

- First pass (cell level):

  `rep_percent` = fraction of non-`NA` bins in state 4; `sigma_chr` =
  population SD of the per-chromosome mean states. A cell fails if
  `rep_percent < min_rep_percent` or is `NA`, `sigma_chr > max_sigma`,
  or at least one chromosome is entirely 2 while another is entirely 4
  (`has_cnv`). Passing cells are ordered by increasing `rep_percent`.

- Second pass (per chromosome):

  Pairwise Pearson correlations of the cells' states on the chromosome
  are computed (values equal to 1 are set to `NA`); cells with a
  negative mean correlation are masked. For each cell with both states
  present on the chromosome, among its `neighbor_window` neighbours on
  each side in `rep_percent` order, the cell is masked if no neighbour
  correlates better with any column than it does, or if the fraction of
  MAD outliers (\|r - median\| \> 3 MAD, IQR rule if MAD = 0) in its row
  exceeds `outlier_threshold`. Masked cells get `NA` on that chromosome;
  a cell masked on at least `n_chromosomes - 1` chromosomes fails the
  second pass.

## Examples

``` r
set.seed(1)
gr <- GenomicRanges::GRanges(rep(c("chr1", "chr2", "chr3"), each = 100),
                             IRanges::IRanges(rep(seq(1, by = 1000, length.out = 100), 3),
                                              width = 1000))
names(gr) <- paste0("bin", seq_along(gr))
rp <- seq(0.1, 0.9, length.out = 20)
st <- sapply(rp, function(p) ifelse(rep(seq(0, 1, length.out = 100), 3) < p, 4, 2))
dimnames(st) <- list(names(gr), paste0("c", 1:20))
md <- data.frame(row.names = colnames(st), cell_id = colnames(st))
res <- compute_rep_percent(st, gr, md, min_rep_percent = 0.01, max_sigma = 0.3,
                           verbose = FALSE)
head(res$meta_data)
#>    cell_id rep_percent sigma_chr pass_qc
#> c1      c1        0.10         0    TRUE
#> c2      c2        0.15         0    TRUE
#> c3      c3        0.19         0    TRUE
#> c4      c4        0.23         0    TRUE
#> c5      c5        0.27         0    TRUE
#> c6      c6        0.31         0    TRUE
```
