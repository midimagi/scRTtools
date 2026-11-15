# Keep autosomal ranges only

Keep autosomal ranges only

## Usage

``` r
keep_autosomes(gr, chromosomes = paste0("chr", 1:19))
```

## Arguments

- gr:

  A `GRanges`.

- chromosomes:

  Chromosomes to keep; the mouse autosomes by default. The seqlevels of
  the result are exactly these.

## Value

The filtered `GRanges`.

## Examples

``` r
gr <- GenomicRanges::GRanges(c("chr1", "chrX", "chr19"),
                             IRanges::IRanges(1, 10))
keep_autosomes(gr)
#> GRanges object with 2 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      1-10      *
#>   [2]    chr19      1-10      *
#>   -------
#>   seqinfo: 19 sequences from an unspecified genome; no seqlengths
```
