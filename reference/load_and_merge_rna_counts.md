# Load and merge RNA count matrices of several samples

Reads one Seurat object per sample from
`file.path(data_dir, sample, file_name)`, extracts the raw counts of
`assay`, prefixes cell IDs with `"<sample>-"` (unless already prefixed),
keeps the genes common to all samples and column-binds the matrices.

## Usage

``` r
load_and_merge_rna_counts(
  sample_names,
  data_dir,
  file_name = "wnn_rna_filtered.rds",
  assay = "rna"
)
```

## Arguments

- sample_names:

  Character vector of sample names (sub-directories of `data_dir`).

- data_dir:

  Directory holding one sub-directory per sample.

- file_name:

  Name of the Seurat `.rds` file inside each sample directory.

- assay:

  Name of the RNA assay.

## Value

A gene-by-cell count matrix (sparse, as stored in the objects).
