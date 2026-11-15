# Infer replication states of one cell

Infer per-window replication states (2 = unreplicated, 4 = replicated)
of a single cell from read counts in G1-derived windows, following the
steps of the original MATLAB-derived pipeline:

## Usage

``` r
infer_single_cell_states(
  counts,
  breakpoints,
  autocorr_aggr_width = 15L,
  autocorr_threshold = 0.2,
  cnv_aggr_width = 15L,
  cnv_fc_cutoff = 1.5,
  cnv_fdr_cutoff = 0.01,
  min_segment_length = 2L,
  hmm_aggr_width = 10L,
  init_ratio_lower = 1.5,
  init_ratio_upper = 2.5,
  min_switch_prob = 1e-04,
  tolerance = 1e-06,
  maxiter = 500L,
  final_ratio_lower = 1.5,
  final_ratio_upper = 2.2
)
```

## Arguments

- counts:

  Numeric vector of read counts per window (names are kept).

- breakpoints:

  Numeric vector of 1-based indices where a new chromosome starts.

- autocorr_aggr_width:

  Aggregation width (windows) for the autocorrelation filter.

- autocorr_threshold:

  Maximum lag-20 autocorrelation.

- cnv_aggr_width:

  Aggregation width for CNV filtering.

- cnv_fc_cutoff:

  Minimum fold change of a CNV segment (must be \> 1).

- cnv_fdr_cutoff:

  FDR cutoff for CNV segments, in (0, 1).

- min_segment_length:

  Minimum segment length (aggregated windows) for CNV segmentation (must
  be \>= 2).

- hmm_aggr_width:

  Aggregation width for HMM initialisation.

- init_ratio_lower, init_ratio_upper:

  Accepted range of the ratio of the two initial Poisson means.

- min_switch_prob:

  Minimum state-switch probability in `[0, 0.5)`.

- tolerance:

  Convergence tolerance of Baum-Welch.

- maxiter:

  Maximum number of Baum-Welch iterations.

- final_ratio_lower, final_ratio_upper:

  Accepted range of the final state-4 / state-2 read-depth ratio.

## Value

A list with `filtered_counts` (numeric, `NaN` at missing or CNV-masked
windows) and `states` (numeric, `2`, `4` or `NA`), both of length
`length(counts)`.

## Details

1.  **CNV filtering**: counts are aggregated by `cnv_aggr_width`
    windows; two-change-point segmentation (linear statistic, BIC
    penalty) plus Poisson tests (BH-adjusted) iteratively mask segments
    whose mean differs from the background by at least
    `cnv_fc_cutoff`-fold with FDR below `cnv_fdr_cutoff`. Masked windows
    become `NaN` in `filtered_counts`.

2.  **Autocorrelation filter**: if the lag-20 autocorrelation of counts
    aggregated by `autocorr_aggr_width` exceeds `autocorr_threshold`,
    all states are `NA`.

3.  **Initialisation**: a two-component Poisson mixture on counts
    aggregated by `hmm_aggr_width`; if the ratio of the two means is
    outside `[init_ratio_lower, init_ratio_upper]`, all states are set
    to `2`.

4.  **Baum-Welch training** (up to `maxiter` iterations, convergence
    tolerance `tolerance`, switch probabilities at least
    `min_switch_prob`); a cell that does not converge gets all states
    `NA`.

5.  **Viterbi decoding**, restarting at `breakpoints` (chromosome
    starts).

6.  **Post-filter**: if the mean-count ratio between state-4 and state-2
    segments is outside `[final_ratio_lower, final_ratio_upper]`, all
    states are `NA`.

## See also

[`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md)
to run this over all cells.
