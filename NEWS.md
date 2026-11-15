# scRTtools 0.0.0.9000

First development version: 111 exported functions for single-cell DNA
replication timing analysis of 10x multiome data. Plotting code is not
included.

## Modules

* Genome bins and counting: `make_fixed_bins()`, `compute_bin_mappability()`,
  `count_reads()`, `filter_windows()`, `calculate_OE_ratio()`,
  `smooth_windows()`, `export_mcol()`.
* Cell QC and G1 windows: `calculate_mapd()`, `find_g1_cells()`,
  `define_G1_windows()`, `synchronize_bins()`.
* Fragment downsampling: `optimize_shrink_factor()`, `downsample_fragments()`
  (wraps `inst/scripts/downsample_fragments_fast.sh`).
* Replication states: `infer_single_cell_states()` (HMM with three Rcpp
  kernels), `assign_replication_states()`, `compute_rep_percent()`,
  `evaluate_state_parameters()`, `rank_state_parameters()`.
* Tracks: `define_replication_tracks()`, `calc_rep_track_lengths()`,
  `calc_rep_track_gaps()`, `count_track_centers()`,
  `summarize_track_lengths()`, `fit_track_velocity_linear()`,
  `fit_track_length_models()`, `fit_track_velocity_ode()`,
  `predict_rep_percent()`.
* Pseudobulk replication timing: `calculate_pseudobulk_rt()`,
  `find_rt_peaks()`, `annotate_peak_scores()`; the S/G1 ratio variant on
  per-cell CP10K, `aggregate_pseudobulk_s_g1_counts()`,
  `calculate_pseudobulk_rt_g1_ratio()` and `export_pseudobulk_g1_ratio()`.
* Initiation zones: `define_initiation_zones()`, `calculate_firing_range()`,
  `calculate_firing_order()`, `predict_firing_order()`,
  `fit_all_izone_velocities()`; eight adjacent-zone functions; 24 zone
  heterogeneity functions (usage, replication tensor, clustering,
  non-negative CP decomposition, epigenomic models).
* Tempo: `rep_tempo()`, `fit_rep_tempo()`, `subset_rep_tempo()`.
* RNA: `SC2P_gene_unit()`, `gex_tempo()`, and 12 functions for S-phase
  progression markers.
* Statistics: 11 group comparison functions; HC3 regressions
  (`fit_ols_hc3()`, `fit_poisson_hc3()`, `combine_two_samples()`,
  `full_rank_model_matrix()`, `add_interaction_terms()`,
  `permute_within_bins()`, `make_bin()`).
* Example data: `scrt_example_data()` (sample scRT1, chr1-chr19, 40 cells).
* Vignette `scRTtools`: the workflow end to end on the example data.

## Design

* Parallel functions take `BPPARAM` (default `BiocParallel::SerialParam()`).
* All inputs are explicit arguments; no function reads a config file or a
  fixed project path.
* Functions return objects and write files only when that is their purpose
  (`export_pseudobulk_g1_ratio(file = )`, `downsample_fragments()`).
* Progress messages go through `verbose`, which defaults to `FALSE`.
* Functions that draw random numbers take a `seed` argument and leave the
  caller's random number state untouched.
* Results are plain `data.frame`s (base R and data.table; no tidyverse
  dependency).
* The HMM and change-point kernels are compiled C++ (Rcpp).

## Testing

* Every module is tested on synthetic data with known answers.
* The functions were validated against the in-house analysis code they were
  derived from, on synthetic data and on a 21-step chain over the example
  data (cell QC to HC3 regressions). These comparison tests need files that
  are not distributed and skip when they are missing.
