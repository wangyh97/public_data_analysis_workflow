#!/usr/bin/env Rscript
# ============================================================================
# Survival analysis -- standardized processed data -> survival result tables
#
# Usage:
#   Rscript analyses/survival/run.R --config config/analyses/survival_example.yaml
#                                   [--output-root /path]
#
# Dataset-agnostic: reads ONLY the processed contract (expression.rds +
# clinical.rds + manifest.yaml). No knowledge of TCGA/Xena/GDC.
#
# INPUT : clinical.rds must contain <ENDPOINT>_time and <ENDPOINT>_event
# OUTPUT: results/<name>/tables/survival_results.tsv   (long-format statistics)
#         results/<name>/tables/km_curves.tsv          (sample-level + step data)
#         results/<name>/metadata/{config_used,analysis_metadata}.yaml
#         results/<name>/logs/run.log ; sessionInfo.txt
#
# Requires the `survival` package (ships with base R as a recommended package).
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

require_pkg("survival", "Kaplan-Meier, log-rank and Cox models")

# ---- helpers ----------------------------------------------------------------

.git_commit_if_available <- function() {
  out <- tryCatch(
    system2("git", c("-C", shQuote(PROJECT_ROOT), "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NULL)
  if (length(out) == 1L && grepl("^[0-9a-f]{7,40}$", out[[1]])) out[[1]] else NULL
}

#' Build the analysis frame for one block: expression of targets + time/event.
.build_block_frame <- function(expr, md, clinical, target_genes_present, endpoint,
                               source, cohort, min_total) {
  t_col <- paste0(endpoint, "_time"); e_col <- paste0(endpoint, "_event")
  missing_cols <- setdiff(c(t_col, e_col), colnames(clinical))
  if (length(missing_cols)) {
    fail(paste0(source, "/", cohort, ": clinical lacks ", paste(missing_cols, collapse = ", "),
                ". Available clinical columns: ", paste(colnames(clinical), collapse = ", "),
                ". This dataset cannot support endpoint '", endpoint, "'."))
  }

  df <- data.frame(
    sample_id  = as.character(clinical$sample_id),
    patient_id = as.character(clinical$patient_id),
    time  = suppressWarnings(as.numeric(clinical[[t_col]])),
    event = suppressWarnings(as.numeric(clinical[[e_col]])),
    stringsAsFactors = FALSE
  )
  for (g in target_genes_present) df[[g]] <- as.numeric(expr[g, match(df$sample_id, colnames(expr))])

  before <- nrow(df)
  df <- df[order(df$patient_id, df$sample_id, method = "radix"), ]
  df <- df[!duplicated(df$patient_id), , drop = FALSE]
  if (nrow(df) < before) {
    log_warn(sprintf("%s/%s: kept first sample for %d multi-sample patients",
                     source, cohort, before - nrow(df)))
  }

  bad_time <- sum(is.na(df$time)); bad_event <- sum(is.na(df$event))
  invalid <- sum(!is.na(df$time) & df$time <= 0, na.rm = TRUE)
  bad_code <- sum(!is.na(df$event) & !df$event %in% c(0, 1), na.rm = TRUE)
  df <- df[!is.na(df$time) & !is.na(df$event) &
             df$time > 0 & df$event %in% c(0, 1), , drop = FALSE]
  dropped_txt <- sprintf("na_time=%d na_event=%d nonpositive_time=%d bad_event_code=%d",
                         bad_time, bad_event, invalid, bad_code)
  if (nrow(df) < before) {
    log_warn(sprintf("%s/%s: excluded %d samples (%s)",
                     source, cohort, before - nrow(df), dropped_txt))
  }
  n_events <- sum(df$event == 1)
  log_info(sprintf("%s/%s [%s]: %d patients, %d events (excluded: %s)",
                   source, cohort, endpoint, nrow(df), n_events, dropped_txt))
  if (nrow(df) < min_total || n_events < 3) {
    fail(sprintf("%s/%s: after cleaning only %d patients / %d events remain (<%d patients or <3 events); refusing to fit",
                 source, cohort, nrow(df), n_events, min_total))
  }
  df
}

.split_groups <- function(df, gene, split_cfg, context) {
  method <- tolower(as.character(split_cfg$method))
  x <- df[[gene]]
  if (method == "fixed") {
    cp <- as.numeric(split_cfg$cutoff)
    if (is.na(cp)) fail(paste0(context, ": split_method=fixed requires survival.cutoff"))
    grp <- ifelse(x >= cp, "High", "Low")
  } else if (method == "median") {
    cp <- stats::median(x, na.rm = TRUE)
    grp <- ifelse(x > cp, "High", "Low")     # strictly above median => deterministic halves
  } else {
    fail(paste0("split_method must be 'median' or 'fixed'; got ", method))
  }
  grp[is.na(x)] <- NA
  df$.expr_raw <- x
  df$.z <- as.numeric(scale(x))
  df$.group <- factor(grp, levels = c("Low", "High"))
  attr(df, "cutpoint") <- cp
  df
}

.logrank_row <- function(df, source, cohort, gene, endpoint, min_per_group) {
  tab <- table(df$.group)
  res <- data.frame(dataset = source, cohort = cohort, target_gene = gene,
                    endpoint = endpoint, model = "logrank",
                    term = "High_vs_Low", estimate_type = "chi2",
                    n = nrow(df), n_events = sum(df$event == 1),
                    estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                    p_value = NA_real_, p_adjust = NA_real_,
                    status = "ok", stringsAsFactors = FALSE)
  if (any(tab < min_per_group)) {
    res$status <- "insufficient_n"
    log_warn(sprintf("%s/%s/%s: log-rank group sizes (%s) below min_per_group=%d",
                     source, cohort, gene, paste(tab, collapse = "/"), min_per_group))
    return(res)
  }
  sd_fit <- survival::survdiff(survival::Surv(time, event) ~ .group, data = df)
  chi2 <- sd_fit$chisq
  p <- stats::pchisq(chi2, df = length(sd_fit$n) - 1L, lower.tail = FALSE)
  res$estimate <- unname(chi2); res$p_value <- unname(p)
  res
}

.cox_rows <- function(df, source, cohort, gene, endpoint, formula_terms, model_name,
                      min_n, min_events) {
  vars_needed <- unique(unlist(formula_terms))
  cc <- stats::complete.cases(df[, c("time", "event", vars_needed), drop = FALSE])
  fit_df <- df[cc, , drop = FALSE]
  n_ev <- sum(fit_df$event == 1)
  base <- data.frame(dataset = source, cohort = cohort, target_gene = gene,
                     endpoint = endpoint, model = model_name,
                     term = NA_character_, estimate_type = "hazard_ratio",
                     n = nrow(fit_df), n_events = n_ev,
                     estimate = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
                     p_value = NA_real_, p_adjust = NA_real_,
                     status = "ok", stringsAsFactors = FALSE)
  out <- list()
  if (nrow(fit_df) < min_n || n_ev < min_events) {
    log_warn(sprintf("%s/%s/%s %s: n=%d events=%d below thresholds -> insufficient_n",
                     source, cohort, gene, model_name, nrow(fit_df), n_ev))
    base$status <- "insufficient_n"
    for (t in formula_terms) {
      r <- base; r$term <- t; out[[length(out) + 1L]] <- r
    }
    return(do.call(rbind, out))
  }
  fml <- stats::as.formula(paste("survival::Surv(time, event) ~",
                                 paste(formula_terms, collapse = " + ")))
  fit <- tryCatch(survival::coxph(fml, data = fit_df, ties = "efron"),
                  error = function(e) {
                    log_warn(sprintf("coxph failed for %s/%s/%s (%s): %s",
                                     source, cohort, gene, model_name, conditionMessage(e)))
                    NULL })
  if (is.null(fit)) {
    for (t in formula_terms) { r <- base; r$term <- t; r$status <- "model_failure"; out[[length(out)+1L]] <- r }
    return(do.call(rbind, out))
  }
  s <- summary(fit)$coefficients
  ci <- summary(fit)$conf.int
  for (t in rownames(s)) {
    r <- base; r$term <- t
    r$estimate <- unname(ci[t, "exp(coef)"])
    r$ci_low   <- unname(ci[t, "lower .95"])
    r$ci_high  <- unname(ci[t, "upper .95"])
    r$p_value  <- unname(s[t, grep("^Pr", colnames(s))[1]])
    out[[length(out) + 1L]] <- r
  }
  do.call(rbind, out)
}

.km_curve_rows <- function(df, source, cohort, gene, endpoint, split_method, cutpoint) {
  sf <- survival::survfit(survival::Surv(time, event) ~ .group, data = df)
  strata_levels <- sub("^[^=]*=", "", as.character(sf$strata))
  rows <- list(); k <- 0L
  for (i in seq_along(strata_levels)) {
    g <- strata_levels[[i]]
    idx <- which(df$.group == g)
    n_g <- length(idx); ev_g <- sum(df$event[idx] == 1)
    times <- c(0, sf$time[sf$strata == i]); surv <- c(1, sf$surv[sf$strata == i])
    nrisk <- c(n_g, sf$n.risk[sf$strata == i]); nevt <- c(0, sf$n.event[sf$strata == i])
    for (j in seq_along(times)) {
      k <- k + 1L
      rows[[k]] <- data.frame(dataset = source, cohort = cohort, target_gene = gene,
                              endpoint = endpoint, group = g,
                              time = times[j], surv = surv[j],
                              n_risk = nrisk[j], n_event = nevt[j],
                              n_group = n_g, events_group = ev_g,
                              split_method = split_method,
                              cutpoint = cutpoint,
                              stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

.main_impl <- function(args) {
  cfg_path <- resolve_path(args$config)
  cfg_raw <- read_yaml_file(cfg_path, "analysis")
  require_fields(cfg_raw, c("analysis.name", "analysis.type", "datasets",
                            "targets", "survival"), context = "analysis config")

  analysis_name <- as.character(get_field(cfg_raw, "analysis.name"))
  if (tolower(as.character(get_field(cfg_raw, "analysis.type"))) != "survival") {
    fail("this module handles type=survival")
  }

  sur_cfg <- get_field(cfg_raw, "survival")
  endpoint <- toupper(as.character(get_field(sur_cfg, "endpoint", "OS")))
  split_cfg <- list(method = get_field(sur_cfg, "split_method", "median"),
                    cutoff = get_field(sur_cfg, "cutoff", NULL))
  covariates <- as.character(get_field(sur_cfg, "covariates", character(0)))
  models_wanted <- toupper(as.character(get_field(sur_cfg, "models", c("logrank", "cox_univariate"))))
  min_per_group <- as.integer(get_field(sur_cfg, "min_samples_per_group", 5))
  min_n <- as.integer(get_field(sur_cfg, "min_n", 10))
  min_total <- as.integer(get_field(sur_cfg, "min_patients", 20))

  out_root <- if (!is.null(args[["output-root"]])) resolve_path(args[["output-root"]])
              else resolve_path(get_field(cfg_raw, "output.root", "${PDA_RESULTS_ROOT}"))
  out_dir <- file.path(out_root, analysis_name)
  tables_dir <- file.path(out_dir, "tables"); meta_dir <- file.path(out_dir, "metadata")
  logs_dir <- file.path(out_dir, "logs")
  for (d in c(tables_dir, meta_dir, logs_dir)) ensure_dir(d)

  log_init(file.path(logs_dir, "run.log"))
  log_section(paste0("Survival analysis: ", analysis_name))
  log_inputs(list(analysis_config = cfg_path))
  log_params(list(endpoint = endpoint, split = split_cfg$method,
                   covariates = (if (length(covariates)) covariates else "<none>"),
                   models = models_wanted,
                   min_per_group = min_per_group, min_n = min_n,
                   min_patients = min_total))

  started_at <- pda_timestamp_utc(); t_start <- Sys.time()

  target_genes <- resolve_gene_spec(get_field(cfg_raw, "targets", list()), role = "target genes")
  validate_gene_symbols(target_genes)

  registry <- load_registry()
  result_rows <- list(); curve_rows <- list(); ri <- 0L; ci_ <- 0L

  for (entry in get_field(cfg_raw, "datasets")) {
    source <- as.character(entry$source); cohort <- as.character(entry$cohort)
    lookup_dataset(registry, source)
    dir_path <- if (!is.null(entry[["processed_dir"]])) resolve_path(entry[["processed_dir"]])
                else processed_dir_for(source, cohort, registry)
    require_file(file.path(dir_path, "manifest.yaml"), paste0(source, "/", cohort, " manifest"))
    expr <- read_rds_input(file.path(dir_path, "expression.rds"), "expression")
    clinical <- read_rds_input(file.path(dir_path, "clinical.rds"), "clinical")
    man <- read_yaml_file(file.path(dir_path, "manifest.yaml"), "processed manifest")

    tg_ok <- require_genes_present(rownames(expr), target_genes, role = "target")
    df <- .build_block_frame(expr, entry, clinical, tg_ok, endpoint, source, cohort, min_total)

    for (gene in tg_ok) {
      dfg <- .split_groups(df, gene, split_cfg,
                           context = paste0(source, "/", cohort, "/", gene))
      n_na_grp <- sum(is.na(dfg$.group))
      if (n_na_grp) {
        log_warn(sprintf("%s/%s/%s: %d samples have NA expression; excluded", source, cohort, gene, n_na_grp))
        dfg <- dfg[!is.na(dfg$.group), , drop = FALSE]
      }
      if ("LOGRANK" %in% models_wanted) {
        ri <- ri + 1L; result_rows[[ri]] <- .logrank_row(dfg, source, cohort, gene, endpoint, min_per_group)
      }
      if ("COX_UNIVARIATE" %in% models_wanted) {
        ri <- ri + 1L
        result_rows[[ri]] <- .cox_rows(dfg, source, cohort, gene, endpoint,
                                       list(".z"), "cox_univariate", min_n, 3L)
      }
      if ("COX_ADJUSTED" %in% models_wanted) {
        terms <- c(".z"); dropped <- character(0)
        for (cv in covariates) {
          if (!cv %in% colnames(clinical)) {
            log_warn(sprintf("%s/%s: covariate '%s' absent from clinical; dropped", source, cohort, cv))
            next
          }
          vals <- clinical[[cv]][match(dfg$sample_id, clinical$sample_id)]
          if (is.numeric(vals) || all(vals %in% c("NA", NA))) {
            num <- suppressWarnings(as.numeric(vals))
            if (all(is.na(num))) { dropped <- c(dropped, cv); next }
            dfg[[paste0("cv_", cv)]] <- num; terms <- c(terms, paste0("cv_", cv))
          } else {
            fac <- factor(trimws(as.character(vals)))
            if (nlevels(fac) < 2 || all(is.na(vals))) { dropped <- c(dropped, cv); next }
            dfg[[paste0("cv_", cv)]] <- droplevels(fac); terms <- c(terms, paste0("cv_", cv))
          }
        }
        if (length(dropped)) log_warn(sprintf("covariates dropped (empty/unusable): %s",
                                              paste(dropped, collapse = ", ")))
        ri <- ri + 1L
        result_rows[[ri]] <- .cox_rows(dfg, source, cohort, gene, endpoint,
                                       as.list(terms), "cox_adjusted", min_n, 3L)
      }
      ci_ <- ci_ + 1L
      curve_rows[[ci_]] <- .km_curve_rows(dfg, source, cohort, gene, endpoint,
                                          split_cfg$method, attr(dfg, "cutpoint"))
    }
  }

  results <- do.call(rbind, result_rows)
  curves <- do.call(rbind, curve_rows)

  # ---- multiple-testing correction ----------------------------------------------
  adj <- toupper(as.character(get_field(sur_cfg, "adjust_method", "BH")))
  fam <- paste(results$dataset, results$cohort, results$target_gene, results$model, sep = "\r")
  results$p_adjust <- NA_real_
  tested <- !is.na(results$p_value)
  for (f in unique(fam[tested])) {
    idx <- which(tested & fam == f)
    results$p_adjust[idx] <- if (adj == "NONE") results$p_value[idx] else
      stats::p.adjust(results$p_value[idx], method = switch(adj, HOLM = "holm", "BH"))
  }

  write_tsv_atomic(results, file.path(tables_dir, "survival_results.tsv"))
  write_tsv_atomic(curves, file.path(tables_dir, "km_curves.tsv"))
  save_config_used(cfg_raw, file.path(meta_dir, "config_used.yaml"))

  pkg_ver <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  write_yaml_atomic(list(
    pipeline_version = PDA_PIPELINE_VERSION,
    git_commit = .git_commit_if_available(),
    analysis = list(name = analysis_name, type = "survival", endpoint = endpoint,
                    split_method = split_cfg$method, cutoff = split_cfg$cutoff,
                    covariates = covariates, models = models_wanted),
    datasets = lapply(get_field(cfg_raw, "datasets"), function(d) d),
    input_files = list(analysis_config = cfg_path),
    output_files = list(survival_results = file.path(tables_dir, "survival_results.tsv"),
                        km_curves = file.path(tables_dir, "km_curves.tsv")),
    counts = list(statistic_rows = nrow(results),
                  ok_rows = sum(results$status == "ok"),
                  km_step_rows = nrow(curves)),
    software = list(R = R.version.string, survival = pkg_ver("survival"),
                    yaml = pkg_ver("yaml"), ggplot2 = pkg_ver("ggplot2")),
    timestamps = list(started_at_utc = started_at,
                      finished_at_utc = pda_timestamp_utc(),
                      duration_s = round(as.numeric(Sys.time() - t_start), 2))
  ), file.path(meta_dir, "analysis_metadata.yaml"))
  writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

  log_outputs(list(survival_results = file.path(tables_dir, "survival_results.tsv"),
                   km_curves = file.path(tables_dir, "km_curves.tsv")))
  cat(sprintf("\n%d statistic rows, %d KM step rows -> %s\n",
              nrow(results), nrow(curves), tables_dir))
}

main <- function() {
  spec <- list(
    config        = list(help = "path to analysis YAML", takes_value = TRUE,
                         required = TRUE, default = NULL),
    `output-root` = list(help = "override output root", takes_value = TRUE,
                         required = FALSE, default = NULL),
    force         = list(help = "accepted for CLI symmetry", takes_value = FALSE,
                         required = FALSE, default = FALSE)
  )
  args <- cli_parse("Run survival analysis on standardized processed data.", spec)
  .main_impl(args)
}

run_main(main)
