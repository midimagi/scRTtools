# Changelog

## scRTtools 0.0.0.9001

- New vignette `figures` and 500-cell vignette data downloaded from a
  GitHub release (see below).
- [`scrt_example_data()`](https://midimagi.github.io/scRTtools/reference/scrt_example_data.md)
  loads the GenomicRanges namespace, so the restored `GRanges` work
  (e.g. [`names()`](https://rdrr.io/r/base/names.html)) in a fresh
  session.

## scRTtools 0.0.0.9000

First development version: 111 exported functions for single-cell DNA
replication timing analysis of 10x multiome data. Plotting code is not
included.

### Modules

- Genome bins and counting:
  [`make_fixed_bins()`](https://midimagi.github.io/scRTtools/reference/make_fixed_bins.md),
  [`compute_bin_mappability()`](https://midimagi.github.io/scRTtools/reference/compute_bin_mappability.md),
  [`count_reads()`](https://midimagi.github.io/scRTtools/reference/count_reads.md),
  [`filter_windows()`](https://midimagi.github.io/scRTtools/reference/filter_windows.md),
  [`calculate_OE_ratio()`](https://midimagi.github.io/scRTtools/reference/calculate_OE_ratio.md),
  [`smooth_windows()`](https://midimagi.github.io/scRTtools/reference/smooth_windows.md),
  [`export_mcol()`](https://midimagi.github.io/scRTtools/reference/export_mcol.md).
- Cell QC and G1 windows:
  [`calculate_mapd()`](https://midimagi.github.io/scRTtools/reference/calculate_mapd.md),
  [`find_g1_cells()`](https://midimagi.github.io/scRTtools/reference/find_g1_cells.md),
  [`define_G1_windows()`](https://midimagi.github.io/scRTtools/reference/define_G1_windows.md),
  [`synchronize_bins()`](https://midimagi.github.io/scRTtools/reference/synchronize_bins.md).
- Fragment downsampling:
  [`optimize_shrink_factor()`](https://midimagi.github.io/scRTtools/reference/optimize_shrink_factor.md),
  [`downsample_fragments()`](https://midimagi.github.io/scRTtools/reference/downsample_fragments.md)
  (wraps `inst/scripts/downsample_fragments_fast.sh`).
- Replication states:
  [`infer_single_cell_states()`](https://midimagi.github.io/scRTtools/reference/infer_single_cell_states.md)
  (HMM with three Rcpp kernels),
  [`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md),
  [`compute_rep_percent()`](https://midimagi.github.io/scRTtools/reference/compute_rep_percent.md),
  [`evaluate_state_parameters()`](https://midimagi.github.io/scRTtools/reference/evaluate_state_parameters.md),
  [`rank_state_parameters()`](https://midimagi.github.io/scRTtools/reference/rank_state_parameters.md).
- Tracks:
  [`define_replication_tracks()`](https://midimagi.github.io/scRTtools/reference/define_replication_tracks.md),
  [`calc_rep_track_lengths()`](https://midimagi.github.io/scRTtools/reference/calc_rep_track_lengths.md),
  [`calc_rep_track_gaps()`](https://midimagi.github.io/scRTtools/reference/calc_rep_track_gaps.md),
  [`count_track_centers()`](https://midimagi.github.io/scRTtools/reference/count_track_centers.md),
  [`summarize_track_lengths()`](https://midimagi.github.io/scRTtools/reference/summarize_track_lengths.md),
  [`fit_track_velocity_linear()`](https://midimagi.github.io/scRTtools/reference/fit_track_velocity_linear.md),
  [`fit_track_length_models()`](https://midimagi.github.io/scRTtools/reference/fit_track_length_models.md),
  [`fit_track_velocity_ode()`](https://midimagi.github.io/scRTtools/reference/fit_track_velocity_ode.md),
  [`predict_rep_percent()`](https://midimagi.github.io/scRTtools/reference/predict_rep_percent.md).
- Pseudobulk replication timing:
  [`calculate_pseudobulk_rt()`](https://midimagi.github.io/scRTtools/reference/calculate_pseudobulk_rt.md),
  [`find_rt_peaks()`](https://midimagi.github.io/scRTtools/reference/find_rt_peaks.md),
  [`annotate_peak_scores()`](https://midimagi.github.io/scRTtools/reference/annotate_peak_scores.md);
  the S/G1 ratio variant on per-cell CP10K,
  [`aggregate_pseudobulk_s_g1_counts()`](https://midimagi.github.io/scRTtools/reference/aggregate_pseudobulk_s_g1_counts.md),
  [`calculate_pseudobulk_rt_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/calculate_pseudobulk_rt_g1_ratio.md)
  and
  [`export_pseudobulk_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/export_pseudobulk_g1_ratio.md).
- Initiation zones:
  [`define_initiation_zones()`](https://midimagi.github.io/scRTtools/reference/define_initiation_zones.md),
  [`calculate_firing_range()`](https://midimagi.github.io/scRTtools/reference/calculate_firing_range.md),
  [`calculate_firing_order()`](https://midimagi.github.io/scRTtools/reference/calculate_firing_order.md),
  [`predict_firing_order()`](https://midimagi.github.io/scRTtools/reference/predict_firing_order.md),
  [`fit_all_izone_velocities()`](https://midimagi.github.io/scRTtools/reference/fit_all_izone_velocities.md);
  eight adjacent-zone functions; 24 zone heterogeneity functions (usage,
  replication tensor, clustering, non-negative CP decomposition,
  epigenomic models).
- Tempo:
  [`rep_tempo()`](https://midimagi.github.io/scRTtools/reference/rep_tempo.md),
  [`fit_rep_tempo()`](https://midimagi.github.io/scRTtools/reference/fit_rep_tempo.md),
  [`subset_rep_tempo()`](https://midimagi.github.io/scRTtools/reference/subset_rep_tempo.md).
- RNA:
  [`SC2P_gene_unit()`](https://midimagi.github.io/scRTtools/reference/SC2P_gene_unit.md),
  [`gex_tempo()`](https://midimagi.github.io/scRTtools/reference/gex_tempo.md),
  and 12 functions for S-phase progression markers.
- Statistics: 11 group comparison functions; HC3 regressions
  ([`fit_ols_hc3()`](https://midimagi.github.io/scRTtools/reference/fit_ols_hc3.md),
  [`fit_poisson_hc3()`](https://midimagi.github.io/scRTtools/reference/fit_poisson_hc3.md),
  [`combine_two_samples()`](https://midimagi.github.io/scRTtools/reference/combine_two_samples.md),
  [`full_rank_model_matrix()`](https://midimagi.github.io/scRTtools/reference/full_rank_model_matrix.md),
  [`add_interaction_terms()`](https://midimagi.github.io/scRTtools/reference/add_interaction_terms.md),
  [`permute_within_bins()`](https://midimagi.github.io/scRTtools/reference/permute_within_bins.md),
  [`make_bin()`](https://midimagi.github.io/scRTtools/reference/make_bin.md)).
- Example data:
  [`scrt_example_data()`](https://midimagi.github.io/scRTtools/reference/scrt_example_data.md)
  (sample scRT1, chr1-chr19, 40 cells); `scrt_example_data("vignette")`
  loads a larger set of 500 cells (not QC filtered, with RNA and a bulk
  Repli-seq RT track), built by `data-raw/make_vignette_data.R`. These
  data (49 MB) are downloaded on first use from the `data-v1` GitHub
  release, checked by MD5 and cached in
  `tools::R_user_dir("scRTtools", "cache")`.
- Vignette `scRTtools`: the workflow end to end on the example data.
- Vignette `figures`: the main figures (MAPD, QC, rep_percent,
  pseudobulk vs bulk RT, single-cell states, UMAP, track-length
  velocity, initiation zones, gene replication vs expression) recomputed
  with ggplot2 on the 500-cell data. It takes about 20-35 minutes to
  build, most of it in
  [`assign_replication_states()`](https://midimagi.github.io/scRTtools/reference/assign_replication_states.md),
  so its code only runs when the environment variable
  `SCRTTOOLS_BUILD_FIGURES=true` is set (e.g. before a release).

### Design

- Parallel functions take `BPPARAM` (default
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)).
- All inputs are explicit arguments; no function reads a config file or
  a fixed project path.
- Functions return objects and write files only when that is their
  purpose (`export_pseudobulk_g1_ratio(file = )`,
  [`downsample_fragments()`](https://midimagi.github.io/scRTtools/reference/downsample_fragments.md)).
- Count data are passed in as objects:
  [`aggregate_pseudobulk_s_g1_counts()`](https://midimagi.github.io/scRTtools/reference/aggregate_pseudobulk_s_g1_counts.md)
  and
  [`export_pseudobulk_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/export_pseudobulk_g1_ratio.md)
  take a sparse matrix or a named list of per-sample matrices (no file
  reading, no `cache_file`);
  [`export_pseudobulk_g1_ratio()`](https://midimagi.github.io/scRTtools/reference/export_pseudobulk_g1_ratio.md)
  takes `genome_granges` only as a `GRanges`. Likewise
  [`load_and_merge_rna_counts()`](https://midimagi.github.io/scRTtools/reference/load_and_merge_rna_counts.md)
  takes a named list of Seurat objects or count matrices,
  [`prepare_metadata()`](https://midimagi.github.io/scRTtools/reference/prepare_metadata.md)
  takes `meta_data` only as a `data.frame`, and
  [`calculate_OE_ratio()`](https://midimagi.github.io/scRTtools/reference/calculate_OE_ratio.md)
  takes `coords_list` as a named list or a `function(chr)` loader. No
  function calls [`readRDS()`](https://rdrr.io/r/base/readRDS.html) on
  its inputs.
- Progress messages go through `verbose`, which defaults to `FALSE`.
- Functions that draw random numbers take a `seed` argument and leave
  the caller’s random number state untouched.
- Results are plain `data.frame`s (base R and data.table; no tidyverse
  dependency).
- The HMM and change-point kernels are compiled C++ (Rcpp).

### Testing

- Every module is tested on synthetic data with known answers.
- The functions were validated against the in-house analysis code they
  were derived from, on synthetic data and on a 21-step chain over the
  example data (cell QC to HC3 regressions). These comparison tests need
  files that are not distributed and skip when they are missing.
