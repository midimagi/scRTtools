# Smooth window values and compute z-scores

Smooth one metadata column of a named `GRangesList` without crossing
genome gaps or large data gaps. Within each gap-free fragment holding at
least `min_wins` non-`NA` windows, values are first averaged with a
centred rolling mean of width 10
(`zoo::rollapply(..., partial = TRUE)`), then a smoothing spline is
fitted to the rolling means against window midpoints. Genome-wide
z-scores of both the rolling mean and the spline fit are added.

## Usage

``` r
smooth_windows(
  wins_grl,
  mcols_id,
  gaps,
  min_wins = 20L,
  gap_th = 100000L,
  smooth_arg = list(spar = 0.2)
)
```

## Arguments

- wins_grl:

  Named `GRangesList` of genomic windows (names must equal the seqnames,
  e.g. `"chr1"`).

- mcols_id:

  Character scalar. Numeric metadata column to smooth, e.g. `"raw_cn"`.

- gaps:

  `GRanges` of genome gap regions.

- min_wins:

  Minimum number of non-`NA` windows required to smooth a fragment
  (integer `>= 2`).

- gap_th:

  Minimum gap size (bp) used to split smoothing fragments.

- smooth_arg:

  List of extra arguments for
  [`stats::smooth.spline()`](https://rdrr.io/r/stats/smooth.spline.html),
  e.g. `list(spar = 0.6)`, `list(cv = TRUE)` or `list(df = 10)`. We use
  `list(spar = 0.1)`.

## Value

`wins_grl` with four added numeric metadata columns: `{mcols_id}_ave`
(rolling mean), `{mcols_id}_smoothed` (spline fit), `{mcols_id}_ave_z`
and `{mcols_id}_smoothed_z` (genome-wide z-scores; all `NA` if the SD is
`NA` or 0). Windows outside smoothed fragments are `NA`.

## Details

Only genome gaps of width `>= gap_th` are used as breaks; a data gap is
declared between consecutive non-`NA` windows whose midpoints are
`>= gap_th` apart.

## Examples

``` r
set.seed(1)
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, by = 1000, length.out = 50),
                                                      width = 1000))
gr$raw_cn <- rnorm(50, 2, 0.2)
gaps <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1))
res <- smooth_windows(GenomicRanges::GRangesList(chr1 = gr), "raw_cn", gaps)
head(res$chr1)
#> GRanges object with 6 ranges and 5 metadata columns:
#>       seqnames    ranges strand |    raw_cn raw_cn_ave raw_cn_smoothed
#>          <Rle> <IRanges>  <Rle> | <numeric>  <numeric>       <numeric>
#>   [1]     chr1    1-1000      * |   1.87471         NA              NA
#>   [2]     chr1 1001-2000      * |   2.03673    2.03133         2.03149
#>   [3]     chr1 2001-3000      * |   1.83287    2.04795         2.04824
#>   [4]     chr1 3001-4000      * |   2.31906    2.05635         2.05412
#>   [5]     chr1 4001-5000      * |   2.06590    2.04330         2.04682
#>   [6]     chr1 5001-6000      * |   1.83591    2.06921         2.06659
#>       raw_cn_ave_z raw_cn_smoothed_z
#>          <numeric>         <numeric>
#>   [1]           NA                NA
#>   [2]     0.294957          0.307066
#>   [3]     0.821072          0.849256
#>   [4]     1.087023          1.039422
#>   [5]     0.674006          0.803076
#>   [6]     1.494067          1.443127
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
