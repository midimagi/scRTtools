# Match cell metadata to RNA cells and flag eligible S-phase cells

Aligns the merged cell metadata with the RNA cell IDs (prefixing
metadata row names with `"<sample>-"` when that matches more cells),
optionally derives a `state` column from sample names, and adds the
flags `rna_pass_qc`, `dna_pass_qc` and `s_phase_pass`. A cell passes
when it passed DNA QC, is not G1, has a finite `rep_percent` within
`[min_rep_percent, max_rep_percent]` and a non-missing group.

## Usage

``` r
prepare_metadata(
  meta_data,
  rna_cell_ids,
  group_by_col = "state",
  min_rep_percent = 0.01,
  max_rep_percent = 0.99,
  state_patterns = c(ama = "-ama-", ctr = "mESC-rep")
)
```

## Arguments

- meta_data:

  A data.frame of cell metadata with cell IDs as row names and columns
  `sample`, `is_G1`, `rep_percent`, `pass_qc`; or the path of an `.rds`
  file holding one.

- rna_cell_ids:

  Cell IDs of the RNA count matrix.

- group_by_col:

  Metadata column defining the analysis groups.

- min_rep_percent, max_rep_percent:

  Range of replication progress kept.

- state_patterns:

  Named character vector of regular expressions used only when
  `group_by_col` is absent: a `state` column is created whose value is
  the name of the first pattern matching `sample` (`NA` when none
  matches). `NULL` disables this.

## Value

`meta_data` restricted to (and ordered like) the cells shared with
`rna_cell_ids`, with the three flag columns added.
