#!/usr/bin/env Rscript
# ============================================================================
# TCGA validate -- structural checks on a processed cohort directory
#
# Usage:
#   Rscript datasets/tcga/validate.R --config config/datasets/tcga.yaml --cohort SKCM
#   Rscript datasets/tcga/validate.R --processed-dir /path/to/processed/TCGA/SKCM
#
# Exit code 0 = all checks pass; 1 = at least one FAIL.
# ============================================================================

PROJECT_ROOT <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(fa)) stop("must be launched with Rscript")
  d <- dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/"))
  repeat {
    if (file.exists(file.path(d, "utils", "common.R"))) break
    p <- dirname(d); if (identical(p, d)) stop("project root not found")
    d <- p
  }
  d
})
source(file.path(PROJECT_ROOT, "utils", "common.R"))

main <- function() {
  spec <- list(
    config        = list(help = "dataset YAML (used with --cohort)", takes_value = TRUE,
                         required = FALSE, default = NULL),
    cohort        = list(help = "cohort code (used with --config)", takes_value = TRUE,
                         required = FALSE, default = NULL),
    `processed-dir` = list(help = "explicit path to processed/<SOURCE>/<cohort>",
                           takes_value = TRUE, required = FALSE, default = NULL),
    `max-na`      = list(help = "maximum allowed overall NA fraction",
                         takes_value = TRUE, required = FALSE, default = "0.5")
  )
  args <- cli_parse("Validate standardized TCGA processed outputs.", spec)

  dir_path <- if (!is.null(args[["processed-dir"]])) {
    resolve_path(args[["processed-dir"]])
  } else {
    if (is.null(args$config) || is.null(args$cohort)) {
      fail("provide either --processed-dir or both --config and --cohort")
    }
    cfg <- read_yaml_file(resolve_path(args$config), "dataset")
    resolve_path(file.path(cfg$paths$processed_root, args$cohort))
  }
  max_na <- as.numeric(args[["max-na"]])

  log_init(NULL)
  log_section("TCGA processed-data validation")
  log_inputs(list(processed_dir = dir_path))

  results <- list()
  check <- function(name, expr) {
    ok <- tryCatch({ expr; TRUE }, error = function(e) { message("     -> ", conditionMessage(e)); FALSE })
    results[[length(results) + 1L]] <<- c(name = name, pass = ok)
    cat(sprintf("[%s] %s\n", if (ok) "PASS" else "FAIL", name))
    invisible(ok)
  }

  req_files <- c(expression.rds = "expression.rds",
                 sample_metadata.rds = "sample_metadata.rds",
                 feature_metadata.rds = "feature_metadata.rds",
                 clinical.rds = "clinical.rds",
                 manifest.yaml = "manifest.yaml")

  check("required files exist and are non-empty", {
    for (f in req_files) require_file(file.path(dir_path, f), f)
  })

  expr <- NULL
  check("expression.rds loads as numeric matrix with unique dimnames", {
    expr <<- read_rds_input(file.path(dir_path, "expression.rds"), "expression")
    pda_assert(is.matrix(expr) && is.numeric(expr), "not a numeric matrix")
    pda_assert(nrow(expr) > 0 && ncol(expr) > 0, "empty matrix")
    pda_assert(!anyDuplicated(rownames(expr)), "duplicated feature labels")
    pda_assert(!anyDuplicated(colnames(expr)), "duplicated sample IDs")
  })

  check("overall NA fraction within limit", {
    pda_assert(sum(is.na(expr)) / length(expr) <= max_na,
               sprintf("NA fraction %.3f exceeds %.2f",
                       sum(is.na(expr)) / length(expr), max_na))
  })

  md <- cl <- fm <- NULL
  check("sample_metadata loads; sample_id matches expression columns exactly", {
    md <<- read_rds_input(file.path(dir_path, "sample_metadata.rds"), "sample metadata")
    validate_sample_metadata(md, c("sample_id"), context = "sample metadata")
    validate_sample_matching(md$sample_id, colnames(expr),
                             context = "metadata vs expression", strict = TRUE)
  })
  check("feature_metadata loads and covers all expression rows", {
    fm <<- read_rds_input(file.path(dir_path, "feature_metadata.rds"), "feature metadata")
    validate_sample_metadata(fm, c("gene_symbol"), context = "feature metadata")
    pda_assert(setequal(fm$gene_symbol, rownames(expr)),
               "feature_metadata gene_symbol set differs from expression rownames")
  })
  check("clinical.rds loads with sample_id/patient_id present", {
    cl <<- read_rds_input(file.path(dir_path, "clinical.rds"), "clinical")
    validate_sample_metadata(cl, c("sample_id", "patient_id"), context = "clinical")
    pda_assert(setequal(cl$sample_id, colnames(expr)),
               "clinical sample set differs from expression columns")
  })
  check("manifest.yaml parses and records fingerprint + pipeline_version", {
    man <- read_yaml_file(file.path(dir_path, "manifest.yaml"), "manifest")
    pda_assert(!is.null(man$fingerprint) && nzchar(man$fingerprint), "missing fingerprint")
    pda_assert(!is.null(man$pipeline_version), "missing pipeline_version")
  })

  n_fail <- sum(!vapply(results, function(r) isTRUE(r[["pass"]]), logical(1)))
  cat(sprintf("\nSummary: %d checks, %d failed\n", length(results), n_fail))
  log_finalize(if (n_fail == 0) "PASS" else "FAIL")
  if (n_fail > 0) .pda_exit(1)
}

run_main(main)
