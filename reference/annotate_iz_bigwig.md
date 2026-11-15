# Mean bigWig signal around IZ centres

Mean bigWig signal around IZ centres

## Usage

``` r
annotate_iz_bigwig(iz_metadata, features, window_bp = 100000L)
```

## Arguments

- iz_metadata:

  Table with `IZ_id`, `chr` and `center`.

- features:

  Named character vector of bigWig files.

- window_bp:

  Window width centred on the IZ centre.

## Value

A `data.frame` with `IZ_id` and one column per feature.
