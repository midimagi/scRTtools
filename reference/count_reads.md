# Count reads per genomic window

Count reads from per-chromosome read-position files and add raw read
counts and raw copy-number estimates to a `GRangesList`.

## Usage

``` r
count_reads(
  wins_grl,
  data_dir,
  bad_coords = NULL,
  raw_cn_th = 6,
  verbose = FALSE
)
```

## Arguments

- wins_grl:

  A named `GRangesList`; each element holds the windows of one
  chromosome, named like `"chr1"`.

- data_dir:

  Directory containing one whitespace-separated read-position file per
  chromosome, named `readsChr<N>.txt` where `<N>` is the element name
  without the `"chr"` prefix.

- bad_coords:

  Optional named list (names = chromosome names, as in `wins_grl`) of
  read positions to drop before counting, or `NULL`.

- raw_cn_th:

  Positive numeric threshold for `raw_cn`.

- verbose:

  Logical; print one progress message per chromosome. Default `FALSE`.

## Value

A list with

- `wins_count`:

  `wins_grl` with metadata columns `raw_counts` (integer) and `raw_cn`
  (numeric).

- `reads_list`:

  named list of numeric read positions per chromosome after removing
  `bad_coords`.

## Details

Reads are assigned to windows with
[`findInterval()`](https://rdrr.io/r/base/findInterval.html) on the cut
points `c(start(windows), end(last window))`, so windows are assumed to
be sorted, contiguous and non-overlapping within each chromosome. Reads
outside the covered span are dropped.
`raw_cn = 2 * raw_counts / mean(raw_counts)`, where the mean is taken
over all windows of all chromosomes; values above `raw_cn_th` or equal
to 0 are set to `NA`.

## Examples

``` r
d <- tempfile(); dir.create(d)
writeLines(as.character(c(5, 15, 15, 25, 35)), file.path(d, "readsChr1.txt"))
wins <- GenomicRanges::GRangesList(
  chr1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 11, 21, 31), width = 10))
)
res <- count_reads(wins, d, verbose = FALSE)
res$wins_count$chr1
#> GRanges object with 4 ranges and 2 metadata columns:
#>       seqnames    ranges strand | raw_counts    raw_cn
#>          <Rle> <IRanges>  <Rle> |  <integer> <numeric>
#>   [1]     chr1      1-10      * |          1       1.6
#>   [2]     chr1     11-20      * |          2       3.2
#>   [3]     chr1     21-30      * |          1       1.6
#>   [4]     chr1     31-40      * |          1       1.6
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
