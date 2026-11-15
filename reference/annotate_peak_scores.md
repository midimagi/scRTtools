# Annotate regions with scores from feature tracks

For each feature file (imported with `rtracklayer::import()`, must have
a `score` column), aggregate the scores of the feature ranges
overlapping each region: `"weighted_mean"` (overlap-width weighted
mean), `"sum"` (sum of score x overlap width) or `"max"`. Regions
without overlap get `NA`.

## Usage

``` r
annotate_peak_scores(
  regions = NULL,
  features,
  window = NULL,
  method = c("weighted_mean", "sum", "max")
)
```

## Arguments

- regions:

  A `GRanges`, or the path of a BED/TSV file with columns `chr`, `start`
  (0-based) and `end` (read with
  [`data.table::fread()`](https://rdrr.io/pkg/data.table/man/fread.html)).

- features:

  Named character vector (or list) of feature file paths; names become
  column names (file basenames are used for unnamed entries, duplicated
  names are made unique).

- window:

  `NULL` or a width (bp) to resize each region to, around its centre
  (then trimmed).

- method:

  Aggregation method.

## Value

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) of the
regions with one column per feature.
