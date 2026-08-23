# ============================================================================
# common.R -- sourced by every entry script after PROJECT_ROOT is resolved.
# Do not source this file directly from other utils.
# ============================================================================

if (!exists("PROJECT_ROOT", inherits = FALSE)) {
  stop("PROJECT_ROOT must be defined before sourcing utils/common.R")
}

for (.module in c("cli.R", "logging.R", "io.R", "config.R", "validation.R", "gene_id.R")) {
  .p <- file.path(PROJECT_ROOT, "utils", .module)
  if (!file.exists(.p)) stop("utility module not found: ", .p)
  source(.p)
}

PDA_PIPELINE_VERSION <- "0.1.0"
