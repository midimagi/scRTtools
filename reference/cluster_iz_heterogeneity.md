# Cluster initiation zones by heterogeneity metrics

Metrics are winsorised (1%-99%), median-imputed and scaled, then
clustered with k-means for every `k` in `k_range`; `k` defaults to the
value with the highest mean silhouette width (computed on at most
`silhouette_n` IZs). Ward hierarchical clustering cut at `k` and a PCA
are returned as well.

## Usage

``` r
cluster_iz_heterogeneity(
  iz_summary,
  metrics,
  k_range = 2:8,
  k = NULL,
  nstart = 50L,
  silhouette_n = 2000L,
  seed = 1L
)
```

## Arguments

- iz_summary:

  Per-IZ table with `IZ_id` and the metric columns (e.g.
  [`assemble_iz_summary()`](https://midimagi.github.io/scRTtools/reference/assemble_iz_summary.md)).

- metrics:

  Metric column names.

- k_range:

  Candidate numbers of clusters.

- k:

  `NULL` or the number of clusters to use (must be in `k_range`).

- nstart:

  `nstart` of [`stats::kmeans()`](https://rdrr.io/r/stats/kmeans.html).

- silhouette_n:

  Maximum number of IZs used for silhouettes.

- seed:

  Seed for the subsampling and k-means.

## Value

A list with `assignments`, `profiles` (metric means per cluster),
`diagnostics`, `z`, `centers`, `hclust`, `pca`, `k` and `metrics`.
