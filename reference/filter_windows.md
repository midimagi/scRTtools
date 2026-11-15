# Segment and filter copy-number windows

For each chromosome, split the windows at genome gaps, segment the
values of `mcols_id` within each gap-free fragment with
[`changepoint::cpt.mean()`](https://rdrr.io/pkg/changepoint/man/cpt.mean.html),
and replace each window value by its segment median (`seg_score`).
Windows are then removed if any of the following holds:

- the window is more than twice as long as the shortest window of the
  chromosome (`end - start` used as length);

- `|seg_score - genome median| > std_th * genome SD` (median and SD of
  the raw `mcols_id` values over all chromosomes);

- `|seg_score - chromosome median of retained seg_score| > std_th * genome SD`;

- `seg_score == 0` or `seg_score` is `NA`.

Only windows with non-`NA` values `<= 20` are segmented, and only
fragments with more than 20 such windows; all other windows get `NA` and
are removed. If `cpt.mean()` fails on a fragment, all its windows get
the fragment median.

## Usage

``` r
filter_windows(
  wins_grl,
  gaps,
  mcols_id,
  std_th = 2,
  cpt_arg = list(method = "PELT", penalty = "MBIC", test.stat = "Normal", minseglen = 5,
    class = TRUE)
)
```

## Arguments

- wins_grl:

  A `GRangesList` (one element per chromosome) containing the numeric
  metadata column `mcols_id`.

- gaps:

  A `GRanges` of genome gap regions.

- mcols_id:

  Character scalar. Metadata column to segment, e.g. `"raw_cn"`.

- std_th:

  Positive numeric scalar. SD multiplier for filtering. We use `0.5`
  (sometimes `0.35`).

- cpt_arg:

  List of arguments passed to
  [`changepoint::cpt.mean()`](https://rdrr.io/pkg/changepoint/man/cpt.mean.html).

## Value

A `GRangesList` with the same names as `wins_grl` holding only the
retained windows; metadata columns are inherited, plus a numeric column
named `paste0(mcols_id, "_seg_score")`.

## Examples

``` r
set.seed(1)
gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, by = 1000, length.out = 60),
                                                      width = 1000))
gr$raw_cn <- c(rnorm(30, 2, 0.1), rnorm(30, 4, 0.1))
gaps <- GenomicRanges::GRanges("chr1", IRanges::IRanges(10001, 12000))
res <- filter_windows(GenomicRanges::GRangesList(chr1 = gr), gaps, "raw_cn",
                      std_th = 2)
length(res$chr1)
#> [1] 48
```
