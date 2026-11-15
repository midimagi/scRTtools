# Synchronise retained bins across samples

Apply a shared bin filter to several samples: for each sample the
*removed* bins are the bins of `gr_cn_list[[s]]` with no exactly
matching bin (same seqnames, start, end, strand) in
`gr_filtered_list[[s]]`; the removed bins of the samples in
`union_samples` are pooled (unique coordinates), and every sample's
retained bins that exactly match a pooled removed bin are dropped. A bin
is therefore kept only if it survived filtering in all `union_samples`
(among the bins present in each sample's unfiltered set).

## Usage

``` r
synchronize_bins(
  gr_cn_list,
  gr_filtered_list,
  union_samples = names(gr_cn_list)
)
```

## Arguments

- gr_cn_list:

  Named list (or `GRangesList`) of the unfiltered per-sample bins.

- gr_filtered_list:

  Named list (or `GRangesList`) of the per-sample retained bins
  (`gr_filtered.rds`, e.g. unlisted
  [`filter_windows()`](https://midimagi.github.io/scRTtools/reference/filter_windows.md)
  output), with the same names as `gr_cn_list`.

- union_samples:

  Names of the samples whose removed bins are pooled. Default: all
  samples.

## Value

A list with

- `filtered`:

  named list of `GRanges`: `gr_filtered_list` with the pooled removed
  bins dropped; metadata columns are kept.

- `removed`:

  named list of `GRanges`: per-sample removed bins.

- `removed_union`:

  `GRanges` of pooled removed bins (no metadata columns).

## Examples

``` r
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, by = 10, length.out = 6),
                                                      width = 10))
res <- synchronize_bins(
  list(a = gr, b = gr),
  list(a = gr[-2], b = gr[-5])
)
res$filtered$a   # bins 2 and 5 removed in both samples
#> GRanges object with 4 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      1-10      *
#>   [2]     chr1     21-30      *
#>   [3]     chr1     31-40      *
#>   [4]     chr1     51-60      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
