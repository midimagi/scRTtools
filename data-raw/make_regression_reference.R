# Reference results for tests/testthat/test-regression-real.R.
#
# Runs the ORIGINAL scripts (rcode/processing, rcode/analysis; the marker
# script and retreat_figure6_helpers.R of the scRT_learning tree; the
# functions of the 6.6 / 6.7 analysis scripts) on the example dataset and
# stores what they return. The file is kept outside the package, in
# $SCRTTOOLS_REFDATA/regression_scRT1_genome.rds. The original scripts are
# found through the environment variables listed in
# tests/testthat/helper-local-paths.R.
#
# Run from the package root:
#   conda run -n new_sc_envs Rscript data-raw/make_regression_reference.R

if (!nzchar(Sys.getenv("SCRTTOOLS_REFDATA"))) stop("set SCRTTOOLS_REFDATA")
suppressMessages(devtools::load_all(".", quiet = TRUE))
for (f in list.files("tests/testthat", pattern = "^helper-.*\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

reference <- regression_pipeline(regression_impl_original())
attr(reference, "created") <- as.character(Sys.time())
attr(reference, "example_md5") <- unname(tools::md5sum(
  system.file("extdata", "scRT1_example.rds", package = "scRTtools")))

out_file <- regression_ref_file()
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
saveRDS(reference, out_file, compress = "xz")
message("Wrote ", out_file, " (", round(file.size(out_file) / 1024), " kB)")
for (nm in names(reference)) {
  x <- reference[[nm]]
  message(sprintf("  %-16s %s", nm, paste(class(x)[1], if (is.null(dim(x))) length(x) else
    paste(dim(x), collapse = " x "))))
}
