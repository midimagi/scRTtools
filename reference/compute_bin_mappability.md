# Compute GC fraction and mappability for fixed genome windows

Compute the per-window GC fraction from a reference FASTA and the
per-window mappable fraction from a mappability track (e.g. a GenMap
bedGraph). Windows whose coordinates fall outside the reference sequence
are retained and get `NA` GC fraction.

## Usage

``` r
compute_bin_mappability(
  genome_bins,
  genome_mappability,
  th = 0.8,
  window_prop = 0.7,
  ref_path,
  BPPARAM = BiocParallel::SerialParam(),
  verbose = FALSE
)
```

## Arguments

- genome_bins:

  A `GRanges` of fixed windows, e.g. from
  [`make_fixed_bins()`](https://midimagi.github.io/scRTtools/reference/make_fixed_bins.md).

- genome_mappability:

  A bedGraph file path, or a `GRanges` with a numeric `score` metadata
  column. Importing a file requires rtracklayer.

- th:

  Numeric scalar in `(0, 1]`. Mappability score threshold.

- window_prop:

  Numeric scalar in `[0, 1]`. A window is mappable when its mappable
  fraction is strictly greater than this value.

- ref_path:

  Character scalar. Path to the reference genome FASTA (requires
  Biostrings). Sequence names may be with or without the `"chr"` prefix.

- BPPARAM:

  A
  [`BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html);
  GC content is computed per chromosome in parallel. We used
  `BiocParallel::MulticoreParam(workers = 8L)`.

- verbose:

  Logical; print progress messages. Default `FALSE`.

## Value

`genome_bins` with added metadata columns:

- `bin_id`:

  character, `"<chr>-<start>-<end>"`.

- `gc_frequency`:

  numeric GC fraction in `[0, 1]` (`NA` for windows outside the
  reference).

- `mappability`:

  numeric mappable fraction in `[0, 1]`.

- `is_mappable`:

  logical, `mappability > window_prop`.

If no mappability interval survives chromosome filtering or the `th`
threshold, the column `frac_mappable = 0` is added instead of
`mappability` and all windows are `is_mappable = FALSE`.

## Details

A base is counted as mappable when it lies in an interval with
`score >= th`; the mappable fraction of a window is the number of such
bases divided by the window width, and `is_mappable` is
`mappability > window_prop` (strictly greater).

## Examples

``` r
if (requireNamespace("Biostrings", quietly = TRUE)) {
  fa <- tempfile(fileext = ".fa")
  seqs <- Biostrings::DNAStringSet(c(
    chr1 = paste(rep("ACGTGGCCAT", 100), collapse = "")
  ))
  Biostrings::writeXStringSet(seqs, fa)
  si <- GenomeInfoDb::Seqinfo("chr1", 1000L)
  bins <- make_fixed_bins(si, 200L, use_chr = 1)
  map <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 501), c(500, 1000)),
                                score = c(1, 0.5))
  GenomeInfoDb::seqinfo(map) <- si
  compute_bin_mappability(bins, map, ref_path = fa)
}
#> GRanges object with 5 ranges and 4 metadata columns:
#>       seqnames    ranges strand |        bin_id gc_frequency mappability
#>          <Rle> <IRanges>  <Rle> |   <character>    <numeric>   <numeric>
#>   [1]     chr1     1-200      * |    chr1-1-200          0.6         1.0
#>   [2]     chr1   201-400      * |  chr1-201-400          0.6         1.0
#>   [3]     chr1   401-600      * |  chr1-401-600          0.6         0.5
#>   [4]     chr1   601-800      * |  chr1-601-800          0.6         0.0
#>   [5]     chr1  801-1000      * | chr1-801-1000          0.6         0.0
#>       is_mappable
#>         <logical>
#>   [1]        TRUE
#>   [2]        TRUE
#>   [3]       FALSE
#>   [4]       FALSE
#>   [5]       FALSE
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome
```
