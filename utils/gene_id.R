# ============================================================================
# gene_id.R -- gene identifier normalization and duplicate-row resolution
#
# Policies (explicit, never hidden):
#   * Ensembl version stripping: remove trailing ".N" only
#   * unmapped Ensembl IDs fall back to the ID itself (flagged in metadata)
#   * multiple rows mapping to one symbol are collapsed by the configured policy
#     (default: max IQR, tie-break lexical identifier) with a full audit table
# ============================================================================

strip_ensembl_version <- function(ids) {
  sub("\\.[0-9]+$", "", as.character(ids))
}

detect_id_type <- function(ids) {
  ids <- as.character(head(ids, 1000))
  if (sum(grepl("^ENSG[0-9]{11}", ids)) > length(ids) * 0.5) "ensembl" else "symbol"
}

row_iqr <- function(m) apply(m, 1L, function(x) stats::IQR(x, na.rm = TRUE))
row_var <- function(m) apply(m, 1L, stats::var, na.rm = TRUE)

# Build an ensembl(stripped) -> symbol vector from a probemap data.frame.
build_symbol_map <- function(probemap, strip_versions = TRUE) {
  id_col  <- intersect(c("id", "probe", "gene_id"), colnames(probemap))[1]
  sym_col <- intersect(c("gene", "geneSymbol", "gene_symbol", "symbol"), colnames(probemap))[1]
  if (is.na(id_col) || is.na(sym_col)) {
    fail(paste0("probemap lacks recognizable id/gene columns; found: ",
                paste(colnames(probemap), collapse = ", ")))
  }
  ids <- as.character(probemap[[id_col]])
  syms <- trimws(as.character(probemap[[sym_col]]))
  ok <- nzchar(syms) & !is.na(syms)
  ids <- ids[ok]; syms <- syms[ok]
  if (strip_versions) ids <- strip_ensembl_version(ids)
  keep_first <- !duplicated(ids)   # deterministic: first occurrence wins
  m <- syms[keep_first]; names(m) <- ids[keep_first]
  m
}

# Collapse matrix rows that share the same feature label (e.g. gene symbol).
# method: max_iqr | max_variance | first_lexical
collapse_duplicate_features <- function(mat, labels, method = "max_iqr") {
  pda_assert(length(labels) == nrow(mat), "labels length must match matrix rows")
  labels <- as.character(labels)
  groups <- split(seq_len(nrow(mat)), factor(labels, levels = unique(labels)))
  sel_rows <- integer(0); audit <- list()
  score_fun <- switch(method,
                      max_variance = function(idx) row_var(mat[idx, , drop = FALSE]),
                      max_iqr      = function(idx) row_iqr(mat[idx, , drop = FALSE]),
                      first_lexical = NULL)
  for (sym in names(groups)) {
    idx <- groups[[sym]]
    if (!is.null(score_fun)) {
      sc <- score_fun(idx)
      sc[is.na(sc)] <- -Inf
      ties <- idx[sc == max(sc)]
      best <- ties[[which.min(rownames(mat)[ties])]]
    } else {
      best <- idx[[order(rownames(mat)[idx])[1]]]
    }
    sel_rows <- c(sel_rows, best)
    audit[[length(audit) + 1L]] <- list(
      feature_label = sym,
      selected_row_id = rownames(mat)[[best]],
      source_row_ids = paste(rownames(mat)[idx], collapse = ";"),
      source_row_count = length(idx),
      selection_method = method
    )
  }
  out <- mat[sel_rows, , drop = FALSE]
  new_order <- order(labels[sel_rows], method = "radix")
  out <- out[new_order, , drop = FALSE]
  rownames(out) <- labels[sel_rows][new_order]
  audit_df <- do.call(rbind, lapply(audit, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  list(matrix = out, audit = audit_df)
}

# Map ensembl rows to symbols; returns list(matrix=..., audit=data.frame).
map_and_collapse_to_symbols <- function(mat, symbol_map, strip_versions = TRUE,
                                        method = "max_iqr") {
  raw_ids <- rownames(mat)
  lookup_ids <- if (strip_versions) strip_ensembl_version(raw_ids) else raw_ids
  symbols <- unname(symbol_map[lookup_ids])
  fallback <- is.na(symbols) | !nzchar(symbols)
  if (any(fallback)) {
    log_warn(sprintf("%d of %d features have no symbol mapping; keeping stripped Ensembl ID as label",
                     sum(fallback), length(fallback)))
    symbols[fallback] <- if (strip_versions) lookup_ids[fallback] else raw_ids[fallback]
  }
  res <- collapse_duplicate_features(mat, symbols, method = method)
  audit <- res$audit
  sel_lookup <- if (strip_versions) strip_ensembl_version(audit$selected_row_id) else audit$selected_row_id
  audit$was_ensembl_fallback <- !(sel_lookup %in% names(symbol_map))
  res$audit <- audit
  res
}
