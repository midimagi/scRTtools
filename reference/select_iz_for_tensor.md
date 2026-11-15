# Select informative initiation zones for the tensor

Score IZs by their mean absolute scaled metric plus relative
`log1p(n_active_cells)`, keep those with at least 10 active cells, and
take the top IZs per k-means cluster, then the top `n_iz` overall.

## Usage

``` r
select_iz_for_tensor(iz_summary, assignments, metrics, n_iz = 80L)
```

## Arguments

- iz_summary:

  Per-IZ table with the metrics and `n_active_cells`.

- assignments:

  `cluster_iz_heterogeneity()$assignments`.

- metrics:

  Metric column names.

- n_iz:

  Number of IZs.

## Value

Integer vector of `IZ_id`.
