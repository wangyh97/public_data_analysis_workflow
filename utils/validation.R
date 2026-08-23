# ============================================================================
# validation.R -- fail-loudly structural checks shared by all layers
#
# Errors stop the run. Data-quality issues that do not invalidate the analysis
# (e.g. a single missing gene) are warnings, per project convention.
# ============================================================================

pda_assert <- function(condition, msg) {
  if (!isTRUE(all(condition) && length(condition) > 0)) fail(msg)
  invisible(TRUE)
}

validate_expression <- function(m, context = "expression", max_na_overall = 0.5,
                                max_na_per_feature = 1.0, verbose = TRUE) {
  pda_assert(is.matrix(m), paste0(context, ": input is not a matrix"))
  pda_assert(is.numeric(m), paste0(context, ": matrix is not numeric (mode=", mode(m), ")"))
  pda_assert(nrow(m) > 0 && ncol(m) > 0,
             paste0(context, ": matrix is empty (dim ", nrow(m), "x", ncol(m), ")"))
  pda_assert(!is.null(colnames(m)) && !anyNA(colnames(m)) && !any(!nzchar(colnames(m))),
             paste0(context, ": sample (column) names are missing or contain NA/empty values"))
  pda_assert(!is.null(rownames(m)) && !anyNA(rownames(m)) && !any(!nzchar(rownames(m))),
             paste0(context, ": feature (row) names are missing or contain NA/empty values"))
  pda_assert(!any(duplicated(colnames(m))), paste0(context, ": duplicated sample IDs present"))
  pda_assert(!any(duplicated(rownames(m))), paste0(context, ": duplicated feature IDs present"))

  n_na <- sum(is.na(m))
  frac <- n_na / length(m)
  stats <- list(n_features = nrow(m), n_samples = ncol(m),
                n_na = n_na, na_fraction = round(frac, 4))
  if (verbose) {
    log_info(sprintf("%s: %d features x %d samples | overall NA fraction = %.4f",
                     context, nrow(m), ncol(m), frac))
  }
  if (!anyNA(m)) return(invisible(stats))

  if (frac > max_na_overall) {
    fail(sprintf("%s: overall NA fraction %.3f exceeds allowed maximum %.2f",
                 context, frac, max_na_overall))
  }
  feat_na_frac <- rowMeans(is.na(m))
  bad <- sum(feat_na_frac >= max_na_per_feature)
  if (bad > 0) {
    fail(sprintf("%s: %d features have NA fraction >= %.2f (e.g. %s)",
                 context, bad, max_na_per_feature,
                 paste(head(names(feat_na_frac)[feat_na_frac >= max_na_per_feature], 5),
                       collapse = ", ")))
  }
  if (frac > 0.05) log_warn(sprintf("%s: overall NA fraction is %.3f (> 0.05)", context, frac))
  invisible(stats)
}

validate_sample_metadata <- function(md, required = c("sample_id"), context = "sample metadata") {
  pda_assert(is.data.frame(md), paste0(context, ": not a data.frame"))
  missing_cols <- setdiff(required, colnames(md))
  if (length(missing_cols)) {
    fail(paste0(context, ": missing required column(s): ", paste(missing_cols, collapse = ", "),
                ". Available columns: ", paste(colnames(md), collapse = ", ")))
  }
  ids <- as.character(md[["sample_id"]])
  pda_assert(all(nzchar(ids)), paste0(context, ": empty sample_id values present"))
  dup <- unique(ids[duplicated(ids)])
  if (length(dup)) {
    fail(paste0(context, ": duplicated sample_id values: ", paste(head(dup, 10), collapse = ", ")))
  }
  invisible(TRUE)
}

validate_clinical <- function(clinical, standard_fields, context = "clinical") {
  validate_sample_metadata(clinical, c("sample_id", "patient_id"), context)
  absent <- setdiff(standard_fields, colnames(clinical))
  for (col in absent) clinical[[col]] <- NA
  available <- standard_fields[vapply(standard_fields, function(f)
    sum(!is.na(clinical[[f]])) > 0, logical(1))]
  missing_now <- setdiff(standard_fields, available)
  if (length(missing_now)) {
    log_warn(paste0(context, ": standard field(s) entirely unavailable in this dataset: ",
                    paste(missing_now, collapse = ", "), " (recorded as NA, never fabricated)"))
  }
  invisible(list(available = available, unavailable = missing_now))
}

# strict=TRUE : sets must be identical (prepared-data contract)
# strict=FALSE: intersection with loud accounting of dropped samples
validate_sample_matching <- function(meta_ids, expr_ids, context = "sample matching",
                                     strict = FALSE) {
  meta_ids <- as.character(meta_ids); expr_ids <- as.character(expr_ids)
  only_meta <- setdiff(meta_ids, expr_ids)
  only_expr <- setdiff(expr_ids, meta_ids)
  keep <- intersect(meta_ids, expr_ids)

  if (length(keep) == 0) {
    fail(paste0(context, ": zero samples match between metadata and expression matrix"))
  }
  if (strict && (length(only_meta) || length(only_expr))) {
    fail(paste0(context, ": strict mode requires identical sample sets; ",
                length(only_meta), " metadata-only, ", length(only_expr), " expression-only"))
  }
  if (length(only_expr)) {
    log_warn(sprintf("%s: %d expression samples lack metadata and were dropped (e.g. %s)",
                     context, length(only_expr), paste(utils::head(only_expr, 5), collapse = ", ")))
  }
  if (length(only_meta)) {
    log_warn(sprintf("%s: %d metadata samples lack expression and were dropped (e.g. %s)",
                     context, length(only_meta), paste(utils::head(only_meta, 5), collapse = ", ")))
  }
  invisible(keep)
}

validate_gene_symbols <- function(genes, context = "gene symbols") {
  suspicious <- genes[!grepl("^[A-Za-z][A-Za-z0-9._@:-]*$", genes)]
  if (length(suspicious)) {
    log_warn(paste0(context, ": entries with unusual characters may not match identifiers: ",
                    paste(head(suspicious, 10), collapse = ", ")))
  }
  invisible(TRUE)
}

# Returns genes actually present; warns about missing ones; errors only when none remain.
require_genes_present <- function(feature_names, requested, role = "requested") {
  requested <- unique(as.character(requested))
  found <- intersect(requested, feature_names)
  missing <- setdiff(requested, found)
  if (length(missing)) {
    log_warn(sprintf("%s gene(s) absent from data (%d of %d): %s",
                     role, length(missing), length(requested),
                     paste(utils::head(missing, 20), collapse = ", ")))
  }
  if (!length(found)) {
    fail(paste0("none of the ", role, " genes are present in the data; cannot continue"))
  }
  invisible(found)
}
