# Define fixed-width genomic windows

Create fixed-width genomic windows for selected chromosomes using
chromosome information from a `Seqinfo` object. Chromosome names are
normalised to the `"chr"`-prefixed style, both in `seqinfo` and in
`use_chr`.

## Usage

``` r
make_fixed_bins(seqinfo, bin_width, use_chr = c(1:19, "X"))
```

## Arguments

- seqinfo:

  A
  [`GenomeInfoDb::Seqinfo`](https://rdrr.io/pkg/Seqinfo/man/Seqinfo-class.html)
  object containing chromosome names and lengths.

- bin_width:

  A single positive integer specifying the window size (bp). We use
  `20000L`.

- use_chr:

  Chromosomes to use, with or without the `"chr"` prefix. Defaults to
  `c(1:19, "X")` (mm10 autosomes + X); use `1:19` for autosomes only.

## Value

A `GRanges` of contiguous, non-overlapping windows
`[1 + (k-1)*bin_width, min(k*bin_width, seqlength)]` per chromosome,
with no metadata columns and with `seqinfo` set to the (chr-normalised)
input `seqinfo`. The last window of each chromosome may be shorter.

## Examples

``` r
si <- GenomeInfoDb::Seqinfo(c("chr1", "chr2"), c(105000L, 60000L))
make_fixed_bins(si, bin_width = 20000L, use_chr = 1:2)
#> GRanges object with 9 ranges and 0 metadata columns:
#>       seqnames        ranges strand
#>          <Rle>     <IRanges>  <Rle>
#>   [1]     chr1       1-20000      *
#>   [2]     chr1   20001-40000      *
#>   [3]     chr1   40001-60000      *
#>   [4]     chr1   60001-80000      *
#>   [5]     chr1  80001-100000      *
#>   [6]     chr1 100001-105000      *
#>   [7]     chr2       1-20000      *
#>   [8]     chr2   20001-40000      *
#>   [9]     chr2   40001-60000      *
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome
```
