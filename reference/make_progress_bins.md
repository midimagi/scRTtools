# Replication-progress bins

Clamp `progress` to `[0, 1]` and cut it into `n_bins` equal bins
labelled `P01`, `P02`, ...

## Usage

``` r
make_progress_bins(progress, n_bins = 10L)
```

## Arguments

- progress:

  Numeric replication progress (fraction).

- n_bins:

  Number of bins.

## Value

A factor.
