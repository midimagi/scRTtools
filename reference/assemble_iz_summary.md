# Join IZ metadata and heterogeneity summaries

Left-join the summaries to `iz_metadata` by `IZ_id` (clashing column
names get `.x` / `.y` suffixes, as
[`dplyr::left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html))
and order by `IZ_id`.

## Usage

``` r
assemble_iz_summary(
  iz_metadata,
  usage_summary,
  timing_summary,
  feature_summary,
  shape_summary
)
```

## Arguments

- iz_metadata:

  `prepare_iz_data()$iz_metadata`.

- usage_summary, timing_summary, feature_summary, shape_summary:

  Per-IZ tables with an `IZ_id` column.

## Value

A `data.frame`.
