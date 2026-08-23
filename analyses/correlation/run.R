#!/usr/bin/env Rscript
# ============================================================================
# Correlation analysis -- standardized processed data -> long-format result table
#
# Usage:
#   Rscript analyses/correlation/run.R --config config/analyses/correlation_example.yaml
#                                      [--output-root /path] [--force]
#
# This module is dataset-agnostic by construction: it knows NOTHING about GDC,
# Xena, TCGAbiolinks, or any API. It only reads the standardized processed-data
# contract (expression.rds + sample_metadata.rds + manifest.yaml).
#
# INPUT : expression matrices via registry layout {processed_root}/{source}/{cohort}
# OUTPUT: results/<analysis.name>/tables/correlation_results.tsv   (long format)
#         results/<analysis.name>/metadata/{analysis_metadata.yaml, config_used.yaml}
#         results/<analysis.name>/logs/run.log, sessionInfo.txt
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

.git_commit_if_available <- function() {
  out <- tryCatch(
    system2("git", c("-C", shQuote(PROJECT_ROOT), "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NULL)
  if (length(out) == 1L && grepl("^[0-9a-f]{7,40}$", out[[1]])) out[[1]] else NULL
}

.load_processed_block <- function(source, cohort, registry, override_dir = NULL) {
  dir_path <- if (!is.null(override_dir)) resolve_path(override_dir)
              else processed_dir_for(source, cohort, registry)
  require_file(file.path(dir_path, "manifest.yaml"), paste0(source, "/", cohort, " manifest"))
  expr <- read_rds_input(file.path(dir_path, "expression.rds"), "expression")
  md   <- read_rds_input(file.path(dir_path, "sample_metadata.rds"), "sample metadata")
  man  <- read_yaml_file(file.path(dir_path, "manifest.yaml"), "processed manifest")
  list(dir = dir_path, expression = expr, metadata = md, manifest = man)
}

.pairwise_correlation_rows <- function(expr, targets, comparisons, method, min_n,
                                        source, cohort) {
  rows <- list()
  k <- 0L
  for (tg in targets) {
    x_full <- expr[tg, ]
    for (cg in comparisons) {
      y_full <- expr[cg, ]
      ok <- !is.na(x_full) & !is.na(y_full)
      n <- sum(ok)
      row <- list(dataset = source, cohort = cohort,
                  target_gene = tg, comparison_gene = cg,
                  n = n, correlation = NA_real_, p_value = NA_real_,
                  p_adjust = NA_real_, method = method,
                  status = NA_character_)
      if (n < min_n) {
        row$status <- "insufficient_n"
        log_warn(sprintf("%s/%s: %s vs %s n=%d < min_n=%d -> insufficient_n",
                         source, cohort, tg, cg, n, min_n))
      } else {
        ct <- suppressWarnings(stats::cor.test(x_full[ok], y_full[ok],
                                               method = method, exact = FALSE))
        row$correlation <- unname(ct$estimate)
        row$p_value <- ct$p.value
        row$status <- "ok"
      }
      k <- k + 1L
      rows[[k]] <- row
    }
  }
  do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
}

.main_impl <- function(args) {
  cfg_path <- resolve_path(args$config)
  cfg_raw <- read_yaml_file(cfg_path, "analysis")
  require_fields(cfg_raw, c("analysis.name", "analysis.type",
                            "datasets", "targets", "comparison"),
                 context = "analysis config")

  analysis_name <- as.character(get_field(cfg_raw, "analysis.name"))
  analysis_type <- tolower(as.character(get_field(cfg_raw, "analysis.type")))
  if (analysis_type != "correlation") {
    fail(paste0("this module handles type=correlation; got '", analysis_type, "'"))
  }

  # ---- output locations -------------------------------------------------------
  out_root <- if (!is.null(args[["output-root"]])) resolve_path(args[["output-root"]])
              else {
                cfg_out <- get_field(cfg_raw, "output.root", "${PDA_RESULTS_ROOT}")
                resolve_path(cfg_out)
              }
  out_dir <- file.path(out_root, analysis_name)
  tables_dir <- file.path(out_dir, "tables")
  meta_dir <- file.path(out_dir, "metadata")
  logs_dir <- file.path(out_dir, "logs")
  for (d in c(tables_dir, meta_dir, logs_dir)) ensure_dir(d)

  log_init(file.path(logs_dir, "run.log"))
  log_section(paste0("Correlation analysis: ", analysis_name))
  log_inputs(list(analysis_config = cfg_path))
  log_params(list(datasets = vapply(get_field(cfg_raw, "datasets"), function(d)
                     paste0(d$source, "/", d$cohort), character(1)),
                   method = get_field(cfg_raw, "correlation.method", "spearman"),
                   adjust_method = get_field(cfg_raw, "correlation.adjust_method", "BH"),
                   adjust_family = get_field(cfg_raw, "correlation.adjust_family", "target_gene"),
                   min_n = get_field(cfg_raw, "correlation.min_n", 10)))

  started_at <- pda_timestamp_utc()
  t_start <- Sys.time()

  # ---- gene specs ---------------------------------------------------------------
  targets_spec <- get_field(cfg_raw, "targets", list())
  comparison_spec <- get_field(cfg_raw, "comparison", list())
  target_genes <- resolve_gene_spec(targets_spec, role = "target genes")
  comparison_genes <- resolve_gene_spec(comparison_spec, role = "comparison genes")
  overlap <- intersect(target_genes, comparison_genes)
  if (length(overlap)) {
    log_warn(paste0("gene(s) appear in both target and comparison lists; ",
                    "self-pairs are skipped: ", paste(overlap, collapse = ", ")))
  }
  validate_gene_symbols(c(target_genes, comparison_genes))

  method <- tolower(as.character(get_field(cfg_raw, "correlation.method", "spearman")))
  if (!method %in% c("spearman", "pearson")) {
    fail(paste0("correlation.method must be spearman|pearson; got ", method))
  }
  adj_method <- toupper(as.character(get_field(cfg_raw, "correlation.adjust_method", "BH")))
  adj_method_flag <- switch(adj_method, BH = "BH", HOLM = "holm", NONE = "none",
                            fail(paste0("unsupported adjust_method: ", adj_method)))
  min_n <- as.integer(get_field(cfg_raw, "correlation.min_n", 10))

  # ---- per-dataset blocks ---------------------------------------------------------
  registry <- load_registry()
  datasets_cfg <- get_field(cfg_raw, "datasets")
  all_blocks <- list(); block_i <- 0L

  for (entry in datasets_cfg) {
    source <- as.character(entry$source); cohort <- as.character(entry$cohort)
    lookup_dataset(registry, source)   # loud failure if adapter not active
    blk <- .load_processed_block(source, cohort, registry,
                                 override_dir = entry[["processed_dir"]])
    md_ids <- as.character(blk$metadata$sample_id)
    keep <- validate_sample_matching(md_ids, colnames(blk$expression),
                                     context = paste0(source, "/", cohort))
    ord <- match(keep, md_ids)
    expr_sub <- blk$expression[, keep, drop = FALSE]
    log_info(sprintf("%s/%s: %d samples x %d features | unit='%s' | fingerprint=%s",
                     source, cohort, ncol(expr_sub), nrow(expr_sub),
                     get_field(blk$manifest, "preprocessing.expression_unit", "?"),
                     substr(get_field(blk$manifest, "fingerprint", "?"), 1, 12)))

    tg_ok <- require_genes_present(rownames(expr_sub), target_genes, role = "target")
    cg_ok <- require_genes_present(rownames(expr_sub), comparison_genes, role = "comparison")

    block <- .pairwise_correlation_rows(expr_sub, tg_ok, cg_ok, method, min_n,
                                        source, cohort)
    # drop self-pairs explicitly
    self_pair <- block$target_gene == block$comparison_gene
    if (any(self_pair)) block <- block[!self_pair, , drop = FALSE]
    block$processed_fingerprint <- get_field(blk$manifest, "fingerprint", NA_character_)
    block_i <- block_i + 1L
    all_blocks[[block_i]] <- block
  }
  if (!length(all_blocks)) fail("no dataset blocks produced rows")
  results <- do.call(rbind, all_blocks)

  # ---- multiple-testing correction --------------------------------------------------
  fam_col <- tolower(as.character(get_field(cfg_raw, "correlation.adjust_family", "target_gene")))
  family_key <- switch(fam_col,
                       target_gene = paste(results$dataset, results$cohort, results$target_gene, sep = "\r"),
                       global = rep("all", nrow(results)),
                       fail(paste0("adjust_family must be 'target_gene' or 'global'; got ", fam_col)))
  results$p_adjust <- NA_real_
  tested <- !is.na(results$p_value)
  for (fam in unique(family_key[tested])) {
    idx <- which(tested & family_key == fam)
    if (adj_method_flag == "none") {
      results$p_adjust[idx] <- results$p_value[idx]
    } else {
      results$p_adjust[idx] <- stats::p.adjust(results$p_value[idx], method = adj_method_flag)
    }
  }

  # ---- write outputs -------------------------------------------------------------------
  results <- results[, c("dataset", "cohort", "target_gene", "comparison_gene",
                         "n", "correlation", "p_value", "p_adjust", "method",
                         "status", "processed_fingerprint")]
  results_path <- file.path(tables_dir, "correlation_results.tsv")
  write_tsv_atomic(results, results_path)

  save_config_used(cfg_raw, file.path(meta_dir, "config_used.yaml"))

  t_end <- Sys.time()
  pkg_ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  metadata <- list(
    pipeline_version = PDA_PIPELINE_VERSION,
    git_commit = .git_commit_if_available(),
    analysis = list(name = analysis_name, type = analysis_type,
                    method = method, adjust_method = adj_method,
                    adjust_family = fam_col, min_n = min_n),
    datasets = lapply(datasets_cfg, function(d) d),
    input_files = list(analysis_config = cfg_path),
    output_files = list(correlation_results = results_path),
    counts = list(rows = nrow(results),
                  tested = sum(tested),
                  insufficient_n = sum(results$status == "insufficient_n")),
    software = list(R = R.version.string,
                    yaml = pkg_ver("yaml"),
                    data.table = pkg_ver("data.table"),
                    ggplot2 = pkg_ver("ggplot2")),
    timestamps = list(started_at_utc = started_at,
                      finished_at_utc = pda_timestamp_utc(),
                      duration_s = round(as.numeric(t_end - t_start), 2))
  )
  write_yaml_atomic(metadata, file.path(meta_dir, "analysis_metadata.yaml"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

  log_outputs(list(result_table = results_path,
                   analysis_metadata = file.path(meta_dir, "analysis_metadata.yaml"),
                   sessionInfo = file.path(out_dir, "sessionInfo.txt"),
                   figures_dir_note = "figures are produced separately by plotting/correlation/plot.R"))
  cat(sprintf("\n%d result rows (%d tested) -> %s\n",
              nrow(results), sum(tested), results_path))
}

main <- function() {
  spec <- list(
    config       = list(help = "path to analysis YAML", takes_value = TRUE,
                        required = TRUE, default = NULL),
    `output-root` = list(help = "override output root (else config/env)",
                         takes_value = TRUE, required = FALSE, default = NULL),
    force        = list(help = "accepted for CLI symmetry; results are recomputed on every run",
                        takes_value = FALSE, required = FALSE, default = FALSE)
  )
  args <- cli_parse("Run correlation analysis on standardized processed data.", spec)
  .main_impl(args)
}

run_main(main)
