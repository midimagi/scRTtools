# Convert a table with a `region_id` column to GRanges

Convert a table with a `region_id` column to GRanges

## Usage

``` r
region_id_to_granges(df)
```

## Arguments

- df:

  data.frame with a `region_id` column formatted `chr:start-end`.

## Value

A `GRanges` with all columns of `df` (including `region_id`) as
metadata.

## Examples

``` r
region_id_to_granges(data.frame(region_id = c("chr1:1-100", "chr2:50-80"),
                                score = 1:2))
#> GRanges object with 2 ranges and 2 metadata columns:
#>       seqnames    ranges strand |   region_id     score
#>          <Rle> <IRanges>  <Rle> | <character> <integer>
#>   [1]     chr1     1-100      * |  chr1:1-100         1
#>   [2]     chr2     50-80      * |  chr2:50-80         2
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
