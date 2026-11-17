# scRTtools

<!-- badges: start -->
[![R-CMD-check](https://github.com/midimagi/scRTtools/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/midimagi/scRTtools/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/midimagi/scRTtools/actions/workflows/pkgdown.yaml/badge.svg)](https://midimagi.github.io/scRTtools/)
<!-- badges: end -->

Documentation website: <https://midimagi.github.io/scRTtools/>

Single-cell DNA replication timing (scRT) analysis for 10x multiome data.

scRTtools covers the analysis from genome binning and cell QC to replication
states, replication tracks, pseudobulk replication timing, initiation zones,
replication tempo and replication-expression statistics. Plotting is not part
of the package: every function returns an object (`GRanges`, matrix,
`data.frame` or list) and writes no files unless that is its stated purpose.

The functions were validated against the in-house analysis code they were
derived from: on synthetic inputs for every module, and on the bundled
example data (chr1 to chr19, 40 cells) for a 21-step chain that runs from
cell QC to the HC3 regressions.

## Installation

The package needs R >= 4.1 and a C++ compiler; it is developed and tested
with R 4.3.3.

Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("midimagi/scRTtools")
```

or from a local clone of the repository:

```bash
R CMD INSTALL scRTtools
```

To keep a Conda environment untouched, install into a library of your own
(`R CMD INSTALL -l ~/R/scRTtools-lib scRTtools`) and load it with
`library(scRTtools, lib.loc = "~/R/scRTtools-lib")`. During development,
`devtools::load_all()` in the package directory loads the package without
installing it.

### Dependencies

- Imports: BiocGenerics, BiocParallel, changepoint, data.table, GenomeInfoDb,
  GenomicRanges, IRanges, Matrix, Rcpp, rlang, Rsamtools, S4Vectors, withr,
  zoo, plus the base packages methods, splines, stats, tools and utils.
- Suggests, each needed by a few functions only and checked when the
  function is called:

  | Package | Needed by |
  |---|---|
  | mgcv | `fit_rep_tempo(model_type = "gam_low" / "gam_fast")`, `fit_all_izone_velocities`, the S-phase progression marker functions |
  | robustbase, quantreg | track length and velocity models (`fit_track_velocity_linear`, `fit_track_length_models`, `fit_track_velocity_ode`) |
  | rtracklayer | `annotate_peak_scores`, `annotate_iz_bigwig`, bigWig input in general |
  | Biostrings | GC content in `make_fixed_bins` |
  | SeuratObject | `load_and_merge_rna_counts` |
  | cluster | silhouette widths in `cluster_iz_heterogeneity` |
  | DNAcopy | `annotate_track_segments` |
  | RhpcBLASctl | limits BLAS threads inside parallel `fit_rep_tempo` workers |

- `downsample_fragments()` runs the bundled shell script
  `inst/scripts/downsample_fragments_fast.sh` and needs `bgzip` and `tabix`
  on the `PATH`.

## Getting started

```r
library(scRTtools)

ex <- scrt_example_data()   # sample scRT1, chr1-chr19, 12 G1 + 28 S-phase cells
names(ex)

vignette("scRTtools")       # the whole workflow on the example data
vignette("figures")         # the main figures, recomputed on 500 cells
```

The `scRTtools` vignette runs the pipeline end to end on the example data and
shows the object each step returns. The `figures` vignette recomputes the main
figures (QC, pseudobulk vs bulk RT, single-cell states, UMAP, track-length
velocity, initiation zones, gene replication vs expression) on a larger set of
500 cells. Those data (49 MB) are not installed with the package:
`scrt_example_data("vignette")` downloads them once from the
[`data-v1` release](https://github.com/midimagi/scRTtools/releases/tag/data-v1)
and caches them. The vignette takes 20-35 minutes, so its code only runs when
`SCRTTOOLS_BUILD_FIGURES=true` is set.

## Overview of the workflow

```
genome            make_fixed_bins -> compute_bin_mappability
                         |
per sample        count_reads / filter_windows -> calculate_mapd -> find_g1_cells
                  -> define_G1_windows -> synchronize_bins (across samples)
                  (optional) optimize_shrink_factor -> downsample_fragments
                         |
states            assign_replication_states (HMM: infer_single_cell_states)
                  -> compute_rep_percent
                  (tuning) evaluate_state_parameters -> rank_state_parameters
                         |
tracks            define_replication_tracks -> calc_rep_track_lengths / _gaps
                  -> summarize_track_lengths -> fit_track_velocity_linear / _ode
                         |
pseudobulk RT     calculate_pseudobulk_rt -> find_rt_peaks -> annotate_peak_scores
                  (S/G1 ratio variant: aggregate_pseudobulk_s_g1_counts ->
                   calculate_pseudobulk_rt_g1_ratio; export_pseudobulk_g1_ratio does both)
                  count_track_centers
                         |
initiation zones  define_initiation_zones -> calculate_firing_range
                  -> calculate_firing_order / predict_firing_order
                  -> fit_all_izone_velocities, adjacent-zone statistics,
                     zone heterogeneity (usage, tensor, clustering, CP, epigenome)
                         |
tempo             rep_tempo -> fit_rep_tempo (sslogis / gam_low / gam_fast)
                         |
RNA               SC2P_gene_unit, gex_tempo, S-phase progression markers
                         |
statistics        group comparison (pivot_state_wide ... call_diff_groups),
                  HC3 regressions (fit_ols_hc3, fit_poisson_hc3,
                  combine_two_samples, permute_within_bins)
```

Conventions used throughout:

- Bins and windows are `GRanges`; per-cell tracks are a `GRangesList` or a
  list of `GRanges`.
- State matrices are bins x cells with copy-number states 2 (unreplicated)
  and 4 (replicated), `NA` where no call was made.
- Functions that loop over cells take `BPPARAM` (default
  `BiocParallel::SerialParam()`).
- Functions that draw random numbers take `seed`.
- Progress messages are controlled by `verbose`, which defaults to `FALSE`
  everywhere; pass `verbose = TRUE` to see progress of the long-running steps
  (counting, G1 calling, state inference, tracks, initiation zones).

## Functions by analysis step

### 1. Bins, cells and replication states

| Step | Package functions |
|---|---|
| fixed-size genome bins, mappability | `make_fixed_bins`, `compute_bin_mappability` |
| per-cell bin counts, cell QC, G1 calling, G1 windows | `count_reads`, `filter_windows`, `calculate_OE_ratio`, `smooth_windows`, `export_mcol`, `calculate_mapd`, `find_g1_cells`, `define_G1_windows` |
| match coverage between conditions by downsampling fragments | `optimize_shrink_factor`, `downsample_fragments` |
| common bins across samples | `synchronize_bins` |
| HMM replication states, per-cell replicated fraction | `assign_replication_states`, `infer_single_cell_states`, `compute_rep_percent` |
| parameter scan for state inference | `evaluate_state_parameters`, `rank_state_parameters` |
| merging several samples | `find_sample_file`, `add_sample_prefix` |

### 2. Tracks

| Step | Package functions |
|---|---|
| replication tracks per cell | `define_replication_tracks`, `calc_rep_track_lengths`, `calc_rep_track_gaps` |
| track length statistics | `summarize_track_lengths` |
| linear length model, inverse prediction | `fit_track_velocity_linear`, `predict_rep_percent` |
| ODE length model | `fit_track_velocity_ode` |
| robust length models by group | `fit_track_length_models` |

### 3. Pseudobulk replication timing and initiation zones

| Step | Package functions |
|---|---|
| pseudobulk RT, RT peaks | `calculate_pseudobulk_rt`, `find_rt_peaks`, `annotate_peak_scores` |
| pseudobulk RT as z-scored S/G1 ratio of per-cell CP10K, over several samples | `aggregate_pseudobulk_s_g1_counts`, `calculate_pseudobulk_rt_g1_ratio`, `export_pseudobulk_g1_ratio` |
| track centre density | `count_track_centers` |
| initiation zones, tracks assigned to zones | `define_initiation_zones` |

### 4. Initiation zone analyses and tempo

| Step | Package functions |
|---|---|
| signal around zones | `annotate_peak_scores` |
| S-phase window in which a zone fires | `calculate_firing_range` |
| firing order | `calculate_firing_order`, `predict_firing_order` |
| co-firing of neighbouring zones | `make_adjacent_pairs`, `make_pair_init`, `build_firing_matrix`, `assign_rep_stage`, `summarise_pair_cofiring_bins`, `summarise_cell_adjacency`, `summarise_pair_adjacency`, `summarise_adjacency_by_stage` |
| fork velocity per zone | `fit_all_izone_velocities`; heatmap helpers `cor_safe`, `zscore_safe`, `peak_bin` |
| zone usage, replication tensor, heterogeneity | `make_progress_bins`, `prepare_iz_data`, `call_iz_usage`, `call_iz_usage_from_states`, `build_replication_tensor`, `correct_signal_matrix`, `correct_replication_tensor`, `compute_usage_heterogeneity`, `fit_timing_models`, `summarize_track_features`, `summarize_track_shape`, `assemble_iz_summary` |
| clustering, neighbour associations, CP decomposition, epigenomic models | `cluster_iz_heterogeneity`, `compute_nearby_iz_associations`, `sample_cells_by_progress`, `select_iz_for_tensor`, `split_signed_residual`, `fit_nonnegative_cp`, `select_cp_rank`, `cp_factor_tables`, `annotate_iz_bigwig`, `correlate_epigenomic_features`, `fit_epigenomic_models`, `fit_epigenomic_classifier` |
| replication tempo curves of genes or zones and their features | `rep_tempo`, `fit_rep_tempo`, `subset_rep_tempo` |

### 5. RNA and replication-expression statistics

| Step | Package functions |
|---|---|
| genes that change along S phase (GAM) | `load_and_merge_rna_counts`, `prepare_metadata`, `normalize_sparse_counts`, `prepare_group_inputs`, `fit_gene_gams`, `build_gene_models`, `summarize_gene_models`, `predict_marginal_curve`, `classify_progression_markers`, `compute_gene_curves`, `scale_rows`, `analyze_one_group` |
| SC2P expression phases | `SC2P_gene_unit` |
| expression along replication progress | `gex_tempo` |
| paired DNA-RNA regressions with HC3 standard errors | `make_bin`, `full_rank_model_matrix`, `add_interaction_terms`, `fit_ols_hc3`, `fit_poisson_hc3`, `combine_two_samples`, `permute_within_bins` |
| comparison of two conditions | `pivot_state_wide`, `add_ma_columns`, `add_rep_cluster`, `call_diff_groups`, `summarise_by_diff_group`, `combine_rep_tempo_feature`, `region_id_to_granges`, `keep_autosomes`, `annotate_track_segments`, `find_intersection`, `format_bp` |

Not included: plotting, the Seurat WNN/UMAP workflow (standard Seurat
steps; the vignette shows how Seurat objects connect to the package) and
transcriptional bursting models.

## Defaults worth knowing

- `fit_rep_tempo`: `n_boot` defaults to 0; bootstrap confidence intervals
  (e.g. `n_boot = 200`) take hours genome-wide.
- `assign_replication_states` uses `infer_args = list(hmm_aggr_width = 15,
  init_ratio_lower = 1.2, final_ratio_lower = 1.4)`; tune these on new data
  with `evaluate_state_parameters`.
- `define_G1_windows(min.gapwidth = 100000)` and
  `define_initiation_zones(min_usage = 3)` are tuned for 20 kb bins and a few
  hundred cells per sample; lower `min_usage` for small datasets.
- `fit_ols_hc3`: a fit whose largest leverage reaches `max_leverage = 0.99`
  returns `NA` with `status = "high_leverage"`; pass `Inf` to never reject a
  fit.
- `combine_two_samples`: with `failed = "keep"` failed fits stay in the result
  (`p = 1`, `status = "fit_failed"`) and count in later multiple-testing
  corrections; `"drop"` removes them.
- Group comparison functions take `states = c("ctr", "ama")`; set them to
  your own condition names.

`NEWS.md` lists the changes between versions.

## Development

```r
devtools::document()
devtools::test()
devtools::check()               # R CMD check --as-cran, builds the vignette
```

A few tests compare the package with the original in-house code and a
real-data regression reference. These files are not part of the repository;
the tests find them through environment variables (see
`tests/testthat/helper-local-paths.R`) and skip when the variables are unset,
as on GitHub Actions. `data-raw/` holds the scripts that build the example
data and the regression reference.

## License

MIT © 2026 Hongjie Ji. See `LICENSE.md`.
