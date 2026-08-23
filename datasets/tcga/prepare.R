#!/usr/bin/env Rscript
# ============================================================================
# TCGA prepare -- raw files -> standardized processed-data contract
#
# Usage:
#   Rscript datasets/tcga/prepare.R --config config/datasets/tcga.yaml [--force]
#                                   [--cohorts SKCM,LUAD]
#
# Outputs (per cohort, under processed_root/<cohort>/):
#   expression.rds         numeric matrix, feature x sample, log2(TPM+0.001)
#   sample_metadata.rds/.tsv
#   feature_metadata.rds/.tsv
#   clinical.rds/.tsv      standardized fields; unavailable ones are NA
#   manifest.yaml          provenance + preprocessing-parameter fingerprint
#
# Reuse rule: if manifest fingerprint matches current parameters and inputs,
#             the cohort directory is reused unless --force.
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

# ---- helpers ----------------------------------------------------------------

.fingerprint_of <- function(lst) {
  txt <- paste(deparse(lst), collapse = "\n")
  tmp <- tempfile(fileext = ".txt")
  writeLines(txt, tmp)
  on.exit(unlink(tmp))
  file_sha256(tmp)$value
}

.gz_read_lines <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con)); readLines(con, warn = FALSE)
}

#' Read a small Xena-style table whose header line may start with '#'.
.gz_read_xena_table <- function(path, label = "table") {
  lines <- .gz_read_lines(path)
  hdr_idx <- which(startsWith(lines, "#"))
  if (length(hdr_idx)) h <- hdr_idx[[1]] else h <- 1L
  hdr <- strsplit(sub("^#\\s*", "", lines[[h]]), "\t", fixed = TRUE)[[1]]
  body <- lines[(h + 1L):length(lines)]
  body <- body[nzchar(trimws(body))]
  body <- body[!startsWith(body, "#")]
  if (!length(body)) fail(paste0(label, ": no data rows found in ", path))
  read.table(text = paste(body, collapse = "\n"), sep = "\t", header = FALSE,
             col.names = hdr, check.names = FALSE, quote = "",
             comment.char = "", stringsAsFactors = FALSE)
}

.cohort_code_from_study <- function(study) {
  m <- regmatches(study, regexpr("\\(([A-Za-z0-9]+)\\)\\s*$", study))
  out <- rep(NA_character_, length(study))
  hit <- lengths(m) > 0
  out[hit] <- gsub("[()]", "", m[hit])
  out
}

.barcode_patient <- function(sample_ids) {
  ids <- as.character(sample_ids)
  out <- ifelse(nchar(ids) >= 12, substr(ids, 1L, 12L), ids)
  out
}

main <- function() {
  spec <- list(
    config  = list(help = "path to dataset YAML", takes_value = TRUE,
                   required = TRUE, default = NULL),
    force   = list(help = "rebuild even when fingerprint matches",
                   takes_value = FALSE, required = FALSE, default = FALSE),
    cohorts = list(help = "comma-separated cohort override, e.g. SKCM,COAD",
                   takes_value = TRUE, required = FALSE, default = NULL)
  )
  args <- cli_parse("Prepare standardized TCGA processed data.", spec)

  cfg_path <- resolve_path(args$config)
  cfg <- read_yaml_file(cfg_path, "dataset")
  require_fields(cfg, c("paths.raw_root", "paths.processed_root", "download.files",
                        "prepare"), context = "tcga.yaml")

  raw_root <- resolve_path(cfg$paths$raw_root)
  proc_root <- resolve_path(cfg$paths$processed_root)

  prep_cfg <- cfg$prepare
  cohorts_req <- if (!is.null(args$cohorts)) {
    strsplit(args$cohorts, ",", fixed = TRUE)[[1]]
  } else as.character(get_field(prep_cfg, "cohorts", list()) %||% c())

  gene_identifier <- get_field(prep_cfg, "gene_identifier", "symbol")
  dup_policy <- get_field(prep_cfg, "duplicate_gene_policy", "max_iqr")
  strip_ver <- isTRUE(get_field(prep_cfg, "strip_ensembl_version", TRUE))
  one_per_patient <- isTRUE(get_field(prep_cfg, "one_sample_per_patient", TRUE))
  priority <- as.character(get_field(prep_cfg, "sample_type_priority", character(0)))
  sample_types_wanted <- as.character(get_field(prep_cfg, "sample_types", character(0)))
  na_drop_threshold <- as.numeric(get_field(prep_cfg, "drop_features_with_na_above", 1.0))

  ensure_dir(proc_root, "processed_root")
  log_init(file.path(proc_root, "logs", "prepare.log"))

  log_section("TCGA prepare")
  log_inputs(list(config = cfg_path,
                  raw_root = raw_root,
                  processed_root = proc_root))
  log_params(list(cohorts_requested = (if (length(cohorts_req)) cohorts_req else "<all>"),
                  sample_types = sample_types_wanted,
                  one_sample_per_patient = one_per_patient,
                  gene_identifier = gene_identifier,
                  duplicate_gene_policy = dup_policy,
                  strip_ensembl_version = strip_ver))

  # ---- 1. locate required raw inputs ---------------------------------------
  files_spec <- cfg$download$files
  role_of <- c(expression_tpm = "expression", phenotype = "phenotype",
               survival = "survival", probemap = "probemap")
  raw_paths <- list(); input_hashes <- list()
  for (key in names(role_of)) {
    if (is.null(files_spec[[key]])) fail(paste0("tcga.yaml download.files lacks entry '", key, "'"))
    pth <- file.path(raw_root, files_spec[[key]]$filename)
    require_file(pth, paste0("raw ", role_of[[key]]))
    raw_paths[[role_of[[key]]]] <- pth
    chk <- file_sha256(pth)
    input_hashes[[key]] <- list(filename = files_spec[[key]]$filename,
                                sha256 = tolower(chk$value),
                                algorithm = chk$algorithm)
  }

  fingerprint <- .fingerprint_of(list(
    pipeline_version = PDA_PIPELINE_VERSION,
    cohorts = sort(cohorts_req),
    sample_types = sample_types_wanted,
    one_sample_per_patient = one_per_patient,
    priority = priority,
    gene_identifier = gene_identifier,
    duplicate_gene_policy = dup_policy,
    strip_ensembl_version = strip_ver,
    drop_features_with_na_above = na_drop_threshold,
    inputs = lapply(input_hashes, `[[`, "sha256")
  ))
  log_info(sprintf("preprocessing parameter fingerprint: %s", fingerprint))

  # ---- 2. probemap -> symbol map --------------------------------------------
  pm <- read_tsv(raw_paths[["probemap"]], label = "probemap")
  symbol_map <- build_symbol_map(pm, strip_versions = strip_ver)
  log_info(sprintf("symbol map: %d Ensembl IDs mapped", length(symbol_map)))

  # ---- 3. phenotype table ----------------------------------------------------
  pheno <- .gz_read_xena_table(raw_paths[["phenotype"]], "phenotype")
  names(pheno)[[1]] <- "sample_id"
  pheno$sample_id <- as.character(pheno$sample_id)
  get_col <- function(df, candidates) {
    hit <- intersect(candidates, colnames(df))[1]
    if (is.na(hit)) return(rep(NA_character_, nrow(df)))
    as.character(df[[hit]])
  }
  pheno$study        <- get_col(pheno, c("_study", "study"))
  pheno$sample_type  <- get_col(pheno, c("_sample_type", "sample_type"))
  pheno$primary_site <- get_col(pheno, c("_primary_site", "primary_site"))
  pheno$sex          <- get_col(pheno, c("_gender", "gender"))
  pheno$disease      <- get_col(pheno, c("primary disease or tissue", "disease"))
  pheno$cohort       <- .cohort_code_from_study(pheno$study)

  is_tcga <- grepl("^TCGA\\b", pheno$study) & !is.na(pheno$study)
  if (!isTRUE(get_field(prep_cfg, "include_non_tcga", FALSE))) pheno <- pheno[is_tcga, , drop = FALSE]
  if (!nrow(pheno)) fail("no TCGA samples found in phenotype table")

  # ---- 4. cohort/sample filtering ---------------------------------------------
  avail <- sort(unique(na.omit(pheno$cohort)))
  wanted <- if (length(cohorts_req)) {
    unknown <- setdiff(cohorts_req, avail)
    if (length(unknown)) {
      fail(paste0("requested cohort(s) not present in data: ",
                  paste(unknown, collapse = ", "), ". Available: ",
                  paste(avail, collapse = ", ")))
    }
    cohorts_req
  } else avail

  sel <- pheno[pheno$cohort %in% wanted &
                 (pheno$sample_type %in% sample_types_wanted), , drop = FALSE]
  absent_types <- setdiff(sample_types_wanted, unique(sel$sample_type))
  if (length(absent_types)) {
    log_warn(paste0("requested sample type(s) with zero samples: ",
                    paste(absent_types, collapse = ", ")))
  }
  # keep only cohorts that survived type filtering
  wanted <- wanted[wanted %in% unique(sel$cohort)]
  if (!length(wanted)) fail("zero samples remain after cohort/sample-type filtering")

  sel$patient_id <- .barcode_patient(sel$sample_id)
  if (one_per_patient && length(priority)) {
    prio <- match(sel$sample_type, priority)
    prio[is.na(prio)] <- length(priority) + 1L
    ord <- order(sel$patient_id, prio, sel$sample_id, method = "radix")
    sel <- sel[ord, ]
    before <- nrow(sel)
    sel <- sel[!duplicated(sel$patient_id), , drop = FALSE]
    log_info(sprintf("one-sample-per-patient dedup: %d -> %d samples", before, nrow(sel)))
  }

  # ---- 5. read expression columns for kept samples ------------------------------
  kept_samples <- sel$sample_id
  expr_df <- read_gz_tsv_columns(raw_paths[["expression"]], keep_cols = kept_samples)
  first_col <- colnames(expr_df)[[1]]
  mat_raw <- as.matrix(expr_df[, -1L, drop = FALSE])
  storage.mode(mat_raw) <- "double"
  rownames(mat_raw) <- as.character(expr_df[[first_col]])
  rm(expr_df); invisible(gc(verbose = FALSE))

  miss_expr <- setdiff(kept_samples, colnames(mat_raw))
  if (length(miss_expr)) {
    log_warn(sprintf("%d filtered samples lack an expression column (dropped)", length(miss_expr)))
    keep_rows <- !(sel$sample_id %in% miss_expr)
    sel <- sel[keep_rows, , drop = FALSE]
  }
  if (!nrow(sel)) fail("zero samples with expression data remain")
  mat_raw <- mat_raw[, sel$sample_id, drop = FALSE]

  # ---- 6. map symbols & collapse duplicates -------------------------------------
  res <- map_and_collapse_to_symbols(mat_raw, symbol_map,
                                     strip_versions = strip_ver, method = dup_policy)
  expr <- res$matrix
  audit <- res$audit
  rm(mat_raw); invisible(gc(verbose = FALSE))

  # optional sparse-feature pruning
  if (na_drop_threshold < 1.0 && anyNA(expr)) {
    frac_na <- rowMeans(is.na(expr))
    dropped <- sum(frac_na > na_drop_threshold)
    if (dropped > 0) {
      log_warn(sprintf("dropping %d features with NA fraction > %.2f", dropped, na_drop_threshold))
      expr <- expr[frac_na <= na_drop_threshold, , drop = FALSE]
    }
  }
  expr_stats <- validate_expression(expr, context = "prepared expression")

  # feature metadata
  feature_meta <- data.frame(
    gene_symbol = rownames(expr),
    selected_row_id = audit$selected_row_id[match(rownames(expr), audit$feature_label)],
    source_row_count = audit$source_row_count[match(rownames(expr), audit$feature_label)],
    selection_method = audit$selection_method[match(rownames(expr), audit$feature_label)],
    stringsAsFactors = FALSE
  )

  # ---- 7. standardized clinical ---------------------------------------------------
  surv_lines <- .gz_read_lines(raw_paths[["survival"]])
  surv_hdr_idx <- which(startsWith(surv_lines, "#"))
  if (length(surv_hdr_idx)) {
    surv_body <- surv_lines[-seq_len(max(surv_hdr_idx))]
    surv_txt <- c(sub("^#\\s*", "", surv_lines[[max(surv_hdr_idx)]]), surv_body)
  } else surv_txt <- surv_lines
  surv <- read.table(text = paste(surv_txt[nzchar(trimws(surv_txt))], collapse = "\n"),
                     sep = "\t", header = FALSE, check.names = FALSE, quote = "",
                     comment.char = "", stringsAsFactors = FALSE,
                     fill = TRUE, na.strings = c("", "NA"))
  # locate key columns tolerantly (first body row was the header line -> drop it)
  hdr_tokens <- strsplit(surv_txt[[1]], "\t", fixed = TRUE)[[1]]
  colnames(surv) <- hdr_tokens[seq_along(hdr_tokens)]
  if (nrow(surv)) surv <- surv[-1L, , drop = FALSE]
  pick <- function(candidates) intersect(candidates, colnames(surv))[1]
  pid_col  <- pick(c("_PATIENT", "_patient_id", "patient"))
  os_t_col <- pick(c("OS.time", "OS_time", "os.time"))
  os_e_col <- pick(c("OS", "OS.event", "os"))
  if (is.na(pid_col)) fail("survival file lacks a patient identifier column")
  surv_key <- if (grepl("_PATIENT|patient", pid_col)) surv[[pid_col]] else .barcode_patient(surv[[pid_col]])

  std_fields <- as.character(get_field(cfg, "clinical.standard_fields",
                                       c("sample_id", "patient_id")))
  clin_base <- data.frame(
    sample_id     = sel$sample_id,
    patient_id    = sel$patient_id,
    cohort        = sel$cohort,
    cancer_type   = sel$study,
    sample_type   = sel$sample_type,
    stringsAsFactors = FALSE
  )
  if ("sex" %in% std_fields) clin_base$sex <- sel$sex
  if (!is.na(os_t_col)) {
    m <- match(clin_base$patient_id, surv_key)
    if ("OS_time" %in% std_fields) clin_base$OS_time <- suppressWarnings(as.numeric(surv[[os_t_col]][m]))
    if ("OS_event" %in% std_fields) clin_base$OS_event <- suppressWarnings(as.numeric(surv[[os_e_col]][m]))
  }
  clinical <- data.frame(row.names = seq_len(nrow(clin_base)), stringsAsFactors = FALSE)
  for (f in std_fields) clinical[[f]] <- if (f %in% colnames(clin_base)) clin_base[[f]] else NA
  clin_avail <- validate_clinical(clinical, std_fields)

  # ---- 8. per-cohort emission ------------------------------------------------------
  log_section("Per-cohort output")
  summary_rows <- list()
  for (coh in sort(wanted)) {
    idx <- which(sel$cohort == coh)
    stopifnot(length(idx) > 0)
    s_ids <- sel$sample_id[idx]
    m_coh <- expr[, s_ids, drop = FALSE]
    md_coh <- sel[idx, , drop = FALSE]
    cl_coh <- clinical[idx, , drop = FALSE]
    fm_coh <- feature_meta[rownames(m_coh), , drop = FALSE]

    out_dir <- file.path(proc_root, coh)
    manifest_path <- file.path(out_dir, "manifest.yaml")

    if (!isTRUE(args$force) && file.exists(manifest_path)) {
      old <- tryCatch(read_yaml_file(manifest_path, "old manifest"), error = function(e) NULL)
      outs_exist <- all(vapply(c("expression.rds", "clinical.rds", "manifest.yaml"),
                               function(f) file.exists(file.path(out_dir, f)) &&
                                 file_size_bytes(file.path(out_dir, f)) > 0, logical(1)))
      if (!is.null(old) && identical(old$fingerprint, fingerprint) && outs_exist) {
        log_info(sprintf("[reuse ] %-8s fingerprint unchanged and outputs verified", coh))
        summary_rows[[length(summary_rows) + 1L]] <-
          data.frame(cohort = coh, status = "reused", n_samples = old$counts$n_samples,
                     n_features = old$counts$n_features, stringsAsFactors = FALSE)
        next
      }
    }

    validate_expression(m_coh, context = sprintf("expression[%s]", coh))
    validate_sample_metadata(md_coh, "sample_id",
                             context = sprintf("metadata[%s]", coh))
    keep_ok <- validate_sample_matching(md_coh$sample_id, colnames(m_coh),
                                        context = sprintf("matching[%s]", coh), strict = TRUE)
    pda_assert(identical(sort(keep_ok), sort(colnames(m_coh))),
               sprintf("matching[%s]: internal inconsistency", coh))

    ensure_dir(out_dir)
    paths <- list(
      expression = file.path(out_dir, "expression.rds"),
      sample_metadata_rds = file.path(out_dir, "sample_metadata.rds"),
      sample_metadata_tsv = file.path(out_dir, "sample_metadata.tsv"),
      feature_metadata_rds = file.path(out_dir, "feature_metadata.rds"),
      feature_metadata_tsv = file.path(out_dir, "feature_metadata.tsv"),
      clinical_rds = file.path(out_dir, "clinical.rds"),
      clinical_tsv = file.path(out_dir, "clinical.tsv")
    )
    write_rds_atomic(m_coh, paths$expression)
    write_rds_atomic(md_coh, paths$sample_metadata_rds)
    write_tsv_atomic(md_coh, paths$sample_metadata_tsv)
    write_rds_atomic(fm_coh, paths$feature_metadata_rds)
    write_tsv_atomic(fm_coh, paths$feature_metadata_tsv)
    write_rds_atomic(cl_coh, paths$clinical_rds)
    write_tsv_atomic(cl_coh, paths$clinical_tsv)

    manifest <- list(
      schema_version = "1.0",
      dataset = get_field(cfg, "dataset.source", "TCGA"),
      cohort = coh,
      created_at_utc = pda_timestamp_utc(),
      pipeline_version = PDA_PIPELINE_VERSION,
      fingerprint = fingerprint,
      preprocessing = list(
        sample_types = sample_types_wanted,
        one_sample_per_patient = one_per_patient,
        sample_type_priority = priority,
        gene_identifier = gene_identifier,
        duplicate_gene_policy = dup_policy,
        strip_ensembl_version = strip_ver,
        expression_unit = get_field(cfg, "expression.unit", "log2(TPM+0.001)"),
        orientation = get_field(cfg, "expression.orientation_after_prepare",
                                "feature x sample")
      ),
      inputs = input_hashes,
      counts = list(n_samples = ncol(m_coh), n_features = nrow(m_coh)),
      clinical_field_availability = list(
        available = clin_avail$available,
        unavailable_as_na = clin_avail$unavailable
      ),
      outputs = lapply(paths, function(p) list(bytes = file_size_bytes(p),
                                               sha256 = file_sha256(p)$value))
    )
    write_yaml_atomic(manifest, manifest_path)

    log_info(sprintf("[write ] %-8s %5d samples x %5d features -> %s",
                     coh, ncol(m_coh), nrow(m_coh), out_dir))
    summary_rows[[length(summary_rows) + 1L]] <-
      data.frame(cohort = coh, status = "written", n_samples = ncol(m_coh),
                 n_features = nrow(m_coh), stringsAsFactors = FALSE)
  }

  log_outputs(list(processed_root = proc_root,
                   fingerprint = fingerprint))
  print(do.call(rbind, summary_rows))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
run_main(main)
