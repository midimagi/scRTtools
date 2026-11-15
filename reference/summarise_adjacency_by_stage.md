# Adjacent-IZ co-firing by replication stage

Run
[`summarise_pair_adjacency()`](https://midimagi.github.io/scRTtools/reference/summarise_pair_adjacency.md)
for all cells and for each stage, and summarise cells and pairs per
stage (means of finite values; median pair phi).

## Usage

``` r
summarise_adjacency_by_stage(
  fired,
  pairs,
  cell_stats,
  rep_stage,
  stages = c("Q1", "Q2", "Q3", "Q4")
)
```

## Arguments

- fired:

  Output of
  [`build_firing_matrix()`](https://midimagi.github.io/scRTtools/reference/build_firing_matrix.md)`$fired`.

- pairs:

  `make_adjacent_pairs()$pairs`.

- cell_stats:

  Per-cell table with `n_fired_izones`, `p11`, `phi`,
  `p_expected_chr_shuffle` (e.g.
  `cbind(cell_info, summarise_cell_adjacency(...))`), rows matching
  `fired`.

- rep_stage:

  Stage of each cell (e.g.
  [`assign_rep_stage()`](https://midimagi.github.io/scRTtools/reference/assign_rep_stage.md)).

- stages:

  Stages to summarise besides `"all"`.

## Value

A list with `pair_stats` (stacked per-pair statistics, `rep_stage` =
`"all"` or a stage) and `stage_summary` (one row per group).
