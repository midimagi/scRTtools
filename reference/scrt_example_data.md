# Example datasets

Loads one of the two example datasets of the package.

## Usage

``` r
scrt_example_data(
  which = c("small", "vignette"),
  cache_dir = tools::R_user_dir("scRTtools", "cache"),
  timeout = 600
)
```

## Arguments

- which:

  `"small"` or `"vignette"`.

- cache_dir:

  Directory holding the downloaded vignette data; default
  `tools::R_user_dir("scRTtools", "cache")`.

- timeout:

  Download timeout in seconds (the `timeout` option is raised to at
  least this value during the download).

## Value

A list with

- `bins`:

  `GRanges` of the 123,148 20 kb bins of chr1-chr19 with `gc_frequency`,
  `mappability` and `is_mappable`.

- `bin_counts`:

  Integer matrix of raw read counts, bins by cells.

- `g1_windows`:

  `GRanges` of the 101,297 variable-width G1 windows of chr1-chr19
  (about 200 G1 reads each) as produced by
  [`define_G1_windows()`](https://midimagi.github.io/scRTtools/reference/define_G1_windows.md).

- `window_counts`:

  Integer matrix of raw read counts, G1 windows by cells.

- `rep_states`:

  Integer matrix of replication states from the original analysis (2 =
  unreplicated, 4 = replicated, `NA` = no call), bins by the 28 S-phase
  cells.

- `rna_counts`:

  Sparse matrix of RNA counts of the 24,860 chr1-chr19 genes, genes by
  cells. Use `meta$nCount_rna` as the library size: these columns hold
  only the chr1-chr19 genes (about 87% of the counts).

- `genes`:

  data.frame with `gene`, `chr`, `start`, `end`, `strand`.

- `meta`:

  data.frame of cell metadata (row names = cell IDs): `is_G1`,
  `rep_percent` (genome-wide fraction replicated), `coverage`, `mapd`,
  RNA depth and QC flags.

- `info`:

  Sample, chromosomes, genome, bin width and creation date.

## Details

`which = "small"` (default): chromosomes 1 to 19 of 40 cells (12 G1 and
28 S-phase) of one mouse ESC 10x multiome sample (`scRT1`, mm10, 20 kb
bins). It is meant for examples, the vignette and tests, not for
biological conclusions. The script that built it is
`data-raw/make_example_data.R` in the package sources.

`which = "vignette"`: the same sample and bins for 500 cells (100 called
G1 and 400 called non-G1 by the original analysis, picked evenly along
their coverage and not filtered on QC, all with RNA), used by the
`figures` vignette. It has the same elements except `rep_states`, plus
`bulk_rt`: a `GRanges` of the 20 kb bins whose `score` is a bulk
Repli-seq replication timing (Loess smoothed; higher = earlier). Built
by `data-raw/make_vignette_data.R`.

The vignette data (49 MB) are not installed with the package. The first
call downloads them from the GitHub release of the package (option
`scRTtools.vignette_data_url`, see below), checks their MD5 sum and
keeps them in `cache_dir`; later calls read the cached copy.

## Examples

``` r
ex <- scrt_example_data()
dim(ex$bin_counts)
#> [1] 123148     40
table(ex$meta$is_G1)
#> 
#> FALSE  TRUE 
#>    28    12 
```
