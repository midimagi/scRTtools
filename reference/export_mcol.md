# Export a GRangesList metadata column to BigWig

Unlist a `GRangesList`, use one metadata column as the BigWig score,
drop `NA` values, attach the supplied `Seqinfo`, sort, and write a
BigWig file with rtracklayer.

## Usage

``` r
export_mcol(grl, mcol_name, file_path, seqinf)
```

## Arguments

- grl:

  A `GRangesList`.

- mcol_name:

  Character scalar. Metadata column to export (coerced with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html)).

- file_path:

  Character scalar. Output BigWig path; its directory must exist.

- seqinf:

  A `Seqinfo` with lengths for every chromosome present in `grl`.

## Value

Invisibly, `file_path`.

## Examples

``` r
if (requireNamespace("rtracklayer", quietly = TRUE)) {
  si <- GenomeInfoDb::Seqinfo("chr1", 10000L)
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 1001), width = 1000),
                               val = c(1.5, NA))
  f <- tempfile(fileext = ".bw")
  export_mcol(GenomicRanges::GRangesList(chr1 = gr), "val", f, si)
}
#> Warning: replacing previous import ‘S4Arrays::makeNindexFromArrayViewport’ by ‘DelayedArray::makeNindexFromArrayViewport’ when loading ‘SummarizedExperiment’
```
