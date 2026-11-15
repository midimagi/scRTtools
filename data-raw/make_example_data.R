# Build the small example dataset shipped in inst/extdata.
#
# Source: sample scRT1 of the scRT results directory (mm10, 20 kb bins),
# genome-wide (all autosomal bins of the 20 kb bin set, chr1-chr19) for 40 cells
# (12 G1 + 28 S-phase). No random numbers are used, so re-running this script on the same inputs gives the same file.
#
# Run from the package root:
#   conda run -n new_sc_envs Rscript data-raw/make_example_data.R
#
# Environment variables (required):
#   SCRTTOOLS_RESULTS  the scRT results directory
#   SCRTTOOLS_GTF      genes.gtf of refdata-cellranger-arc-mm10-2020-A-2.0.0

t_start <- proc.time()[["elapsed"]]
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Matrix)
})

results_dir <- Sys.getenv("SCRTTOOLS_RESULTS")
gtf_file <- Sys.getenv("SCRTTOOLS_GTF")
if (!nzchar(results_dir) || !nzchar(gtf_file)) {
  stop("set SCRTTOOLS_RESULTS and SCRTTOOLS_GTF")
}
sample_name <- "scRT1"
n_g1 <- 12L
n_s <- 28L
out_file <- file.path("inst", "extdata", "scRT1_example.rds")

sample_dir <- file.path(results_dir, "1.data", sample_name)
merge_dir <- file.path(results_dir, "1.data", "merge")

# ---- cells -----------------------------------------------------------------
meta <- as.data.frame(readRDS(file.path(merge_dir, "meta_data.rds")))
rownames(meta) <- meta$cell_id
usable <- meta$sample == sample_name & meta$pass_qc %in% TRUE & !(meta$has_cnv %in% TRUE)

# evenly spaced picks along an ordering, without random numbers
pick_evenly <- function(ids, n) ids[unique(round(seq(1, length(ids), length.out = n)))]

g1_pool <- meta[usable & meta$is_G1 %in% TRUE, ]
g1_cells <- pick_evenly(g1_pool$cell_id[order(g1_pool$coverage, g1_pool$cell_id)], n_g1)
s_pool <- meta[usable & meta$is_G1 %in% FALSE & is.finite(meta$rep_percent), ]
s_cells <- pick_evenly(s_pool$cell_id[order(s_pool$rep_percent, s_pool$cell_id)], n_s)
cells <- c(g1_cells, s_cells)
stopifnot(length(cells) == n_g1 + n_s, !anyDuplicated(cells))

# ---- genome bins and raw counts ---------------------------------------------
bins <- readRDS(file.path(results_dir, "0.mm10", "genobins_20kb_granges.rds"))
chromosomes <- seqlevels(bins)
message("bins: ", length(bins), " on ", paste(chromosomes, collapse = ","))

subset_matrix <- function(file, row_ids, col_ids) {
  x <- readRDS(file)
  message(basename(file), ": ", class(x)[1], " ", paste(dim(x), collapse = " x "))
  # per-sample files name cells by bare barcode, merged files by "<sample>-<barcode>"
  lookup <- col_ids
  if (!all(lookup %in% colnames(x))) lookup <- sub(paste0("^", sample_name, "-"), "", col_ids)
  if (!all(row_ids %in% rownames(x)) || !all(lookup %in% colnames(x))) {
    stop(basename(file), ": requested rows or cells are missing; e.g. rows ",
         paste(utils::head(rownames(x), 2), collapse = ", "), "; cells ",
         paste(utils::head(colnames(x), 2), collapse = ", "))
  }
  x <- x[row_ids, lookup, drop = FALSE]
  colnames(x) <- col_ids
  if (inherits(x, "Matrix") || inherits(x, "IterableMatrix")) x <- as.matrix(x)
  x
}

bin_counts <- subset_matrix(file.path(sample_dir, "genobins_20kb_raw_matrix.rds"),
                            names(bins), cells)
storage.mode(bin_counts) <- "integer"
gc()

g1_windows <- readRDS(file.path(sample_dir, "g1bins_granges.rds"))
message("g1bins_granges: ", class(g1_windows)[1], " n=", length(g1_windows),
        " mcols: ", paste(names(mcols(g1_windows)), collapse = ","))
stopifnot(all(seqlevels(g1_windows) %in% chromosomes))
window_counts <- subset_matrix(file.path(sample_dir, "g1bins_raw_matrix.rds"),
                               names(g1_windows), cells)
storage.mode(window_counts) <- "integer"
gc()

# ---- replication states from the original pipeline (reference) --------------
rep_states <- subset_matrix(file.path(merge_dir, "genobins_rep_states.rds"),
                            names(bins), s_cells)
storage.mode(rep_states) <- "integer"

# ---- RNA counts of the genes on the bin chromosomes -------------------------
rna <- readRDS(file.path(merge_dir, "merge_rna_counts_filtered.csv"))  # an RDS despite its name
message("rna: ", class(rna)[1], " ", paste(dim(rna), collapse = " x "))
gtf <- data.table::fread(
  cmd = paste("awk -F '\t'", shQuote('$3 == "gene"'),
              shQuote(gtf_file)),
  header = FALSE, sep = "\t", quote = "", showProgress = FALSE
)
genes <- data.frame(
  gene = sub('.*gene_name "([^"]+)".*', "\\1", gtf$V9),
  chr = gtf$V1, start = as.integer(gtf$V4), end = as.integer(gtf$V5), strand = gtf$V7,
  stringsAsFactors = FALSE
)
genes <- genes[genes$chr %in% chromosomes, ]
genes <- genes[!duplicated(genes$gene) & genes$gene %in% rownames(rna), ]
genes <- genes[order(match(genes$chr, chromosomes), genes$start, genes$gene), ]
rownames(genes) <- NULL
rna_cells <- intersect(cells, colnames(rna))
# the merged RNA table is a data.frame of counts (genes x cells)
rna_counts <- Matrix::Matrix(as.matrix(rna[genes$gene, rna_cells, drop = FALSE]), sparse = TRUE)
rna_counts <- methods::as(rna_counts, "generalMatrix")

# ---- metadata ---------------------------------------------------------------
keep_cols <- c("cell_id", "sample", "nCount_rna", "nFeature_rna", "coverage", "mapd",
               "scaled_mapd", "is_G1", "rep_percent", "sigma_chr", "has_cnv", "pass_qc")
cell_meta <- meta[cells, keep_cols]
cell_meta$has_rna <- cells %in% rna_cells

example <- list(
  bins = bins,
  bin_counts = bin_counts,
  g1_windows = g1_windows,
  window_counts = window_counts,
  rep_states = rep_states,
  rna_counts = rna_counts,
  genes = genes,
  meta = cell_meta,
  info = list(sample = sample_name, chromosomes = chromosomes, genome = "mm10",
              bin_width = 20000L, created = as.character(Sys.Date()))
)
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(example, out_file, compress = "xz")

message("\nWrote ", out_file, " (", round(file.size(out_file) / 1024), " kB)")
for (nm in names(example)) {
  x <- example[[nm]]
  message(sprintf("  %-14s %-12s %s  %s kB in memory", nm, class(x)[1],
                  paste(if (is.null(dim(x))) length(x) else dim(x), collapse = " x "),
                  round(as.numeric(object.size(x)) / 1024)))
}
print(table(G1 = cell_meta$is_G1, has_rna = cell_meta$has_rna))
print(summary(cell_meta$rep_percent[!cell_meta$is_G1]))
print(summary(colSums(bin_counts)))
print(table(rep_states, useNA = "a"))
print(g1_windows[1:3])
print(summary(Matrix::colSums(rna_counts)))
print(table(genes$chr)[chromosomes])
message(sprintf("elapsed: %.1f s", proc.time()[["elapsed"]] - t_start))
