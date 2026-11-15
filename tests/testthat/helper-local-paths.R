# Directories with the original scripts / notebooks / reference data used by
# the comparison tests. They are only read from environment variables; when a
# variable is unset the path points to a non-existent directory and the tests
# that need it skip.
#
#   SCRTTOOLS_ORIG_PROCESSING      scRT_src rcode/processing (rcode = its parent)
#   SCRTTOOLS_ORIG_LEARNING_RCODE  scRT_learning rcode
#   SCRTTOOLS_NOTEBOOKS            scRT_src jupyter notebooks
#   SCRTTOOLS_NOTEBOOKS_DS         scRT_learning jupyter_downsample_frag_quantile
#   SCRTTOOLS_REFDATA              regression reference (regression_scRT1_genome.rds)
local_dir <- function(var) {
  d <- Sys.getenv(var)
  if (nzchar(d)) d else file.path(tempdir(), paste0("unset_", var))
}
