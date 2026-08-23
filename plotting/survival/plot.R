#!/usr/bin/env Rscript
# ============================================================================
# Survival plotting -- result tables -> figures (NO statistics recomputation)
#
# Usage:
#   Rscript plotting/survival/plot.R \
#       --results results/<name>/tables/survival_results.tsv \
#       --curves  results/<name>/tables/km_curves.tsv \
#       [--plot-config config/plotting/default.yaml] \
#       [--out-dir results/<name>/figures] [--types km,forest]
#
# Outputs: survival_km.<pdf|png>   Kaplan-Meier curves per block
#          survival_forest.<pdf|png> Cox hazard-ratio forest plot
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

require_pkg("ggplot2", "figure rendering")

main <- function() {
  spec <- list(
    results       = list(help = "survival_results.tsv from analyses/survival",
                         takes_value = TRUE, required = TRUE, default = NULL),
    curves        = list(help = "km_curves.tsv from analyses/survival",
                         takes_value = TRUE, required = TRUE, default = NULL),
    `plot-config` = list(help = "plotting YAML", takes_value = TRUE,
                         required = FALSE, default = "config/plotting/default.yaml"),
    `out-dir`     = list(help = "figure output directory", takes_value = TRUE,
                         required = FALSE, default = NULL),
    types         = list(help = "comma list among km,forest", takes_value = TRUE,
                         required = FALSE, default = "km,forest"),
    prefix        = list(help = "filename prefix base", takes_value = TRUE,
                         required = FALSE, default = "survival")
  )
  args <- cli_parse("Plot KM curves and a Cox forest plot from survival tables.", spec)

  res_path <- resolve_path(args$results)
  cur_path <- resolve_path(args$curves)
  out_dir <- if (!is.null(args[["out-dir"]])) resolve_path(args[["out-dir"]])
             else file.path(dirname(dirname(res_path)), "figures")
  ensure_dir(out_dir)
  log_init(file.path(out_dir, "plot_survival.log"))

  pcfg <- read_yaml_file(resolve_path(args[["plot-config"]]), "plot config")

  log_section("Survival plots")
  log_inputs(list(results = res_path, curves = cur_path, plot_config = args[["plot-config"]]))
  log_params(list(types = args$types, out_dir = out_dir))

  res <- read_tsv(res_path, label = "survival results")
  cur <- read_tsv(curves, label = "KM curves")
  need_res <- c("dataset", "cohort", "target_gene", "endpoint", "model",
                "term", "estimate_type", "estimate", "ci_low", "ci_high",
                "p_value", "p_adjust", "n", "n_events")
  missing_cols <- setdiff(need_res, colnames(res))
  if (length(missing_cols)) fail(paste0("survival_results.tsv lacks: ",
                                        paste(missing_cols, collapse = ", ")))
  need_cur <- c("dataset", "cohort", "target_gene", "endpoint", "group",
                "time", "surv", "n_risk", "n_group", "events_group", "cutpoint")
  missing_cols <- setdiff(need_cur, colnames(cur))
  if (length(missing_cols)) fail(paste0("km_curves.tsv lacks: ",
                                        paste(missing_cols, collapse = ", ")))

  for (col in c("time", "surv")) cur[[col]] <- suppressWarnings(as.numeric(cur[[col]]))
  for (col in c("estimate", "ci_low", "ci_high", "p_value", "p_adjust"))
    res[[col]] <- suppressWarnings(as.numeric(res[[col]]))

  res$block <- paste(res$dataset, res$cohort, res$endpoint, res$target_gene, sep = "\n")
  cur$block <- paste(cur$dataset, cur$cohort, cur$endpoint, cur$target_gene, sep = "\n")
  cur$group <- factor(cur$group, levels = c("Low", "High"))

  base_size <- as.numeric(get_field(pcfg, "heatmap.font_size", 10))
  formats <- tolower(as.character(get_field(pcfg, "figure.formats", c("pdf", "png"))))
  dpi <- as.integer(get_field(pcfg, "figure.dpi", 600))
  width_in <- as.numeric(get_field(pcfg, "figure.width_in", 9))
  height_in <- as.numeric(get_field(pcfg, "figure.height_in", 5))

  save_all <- function(plot_obj, prefix) {
    written <- character(0)
    for (fmt in formats) {
      path <- file.path(out_dir, paste0(prefix, ".", fmt))
      ok <- tryCatch({ ggplot2::ggsave(path, plot_obj, width = width_in,
                                       height = height_in, dpi = dpi); TRUE },
                     error = function(e) { log_warn(paste0("failed saving ", fmt,
                                                           ": ", conditionMessage(e))); FALSE })
      if (ok) written <- c(written, path)
    }
    written
  }
  outputs <- character(0)
  types <- strsplit(args$types, ",", fixed = TRUE)[[1]]
  types <- trimws(tolower(types))

  if ("km" %in% types) {
    # annotation: log-rank p per block (read from results table; no recomputation)
    lr <- res[res$model == "logrank", , drop = FALSE]
    blocks <- unique(cur$block)
    caps <- vapply(blocks, function(b) {
      r <- lr[lr$block == b & !is.na(lr$p_value), , drop = FALSE]
      if (!nrow(r)) return(paste0(b, ": no valid log-rank result"))
      q <- if (is.na(r$p_adjust[[1]])) NA else r$p_adjust[[1]]
      sprintf("%s: log-rank p=%.2g%s", b, r$p_value[[1]],
              ifelse(is.na(q), "", sprintf(", q=%.2g", q)))
    }, character(1))

    group_meta <- unique(cur[, c("block", "group", "n_group", "events_group")])
    cur$legend <- vapply(seq_len(nrow(cur)), function(i) {
      g <- group_meta[group_meta$block == cur$block[i] &
                        group_meta$group == cur$group[i], ][1, ]
      if (is.na(g$group)) return(as.character(cur$group[i]))
      sprintf("%s (n=%d, events=%d)", cur$group[i], g$n_group, g$events_group)
    }, character(1))

    legend_vals <- unique(cur[, c("group", "legend")])
    pal <- c(Low = "#2166AC", High = "#B2182B")
    names(pal) <- legend_vals$legend[match(names(pal), legend_vals$group)]
    pal <- pal[!is.na(names(pal))]

    p_km <- ggplot2::ggplot(cur, ggplot2::aes(x = time, y = surv, color = legend)) +
      ggplot2::geom_step(linewidth = 0.8) +
      ggplot2::facet_wrap(~block) +
      ggplot2::scale_color_manual(values = pal, name = NULL) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"),
                     legend.position = "bottom") +
      ggplot2::labs(x = "Time (source units)", y = "Survival probability",
                    title = "Kaplan-Meier by target-gene expression group",
                    caption = paste(caps, collapse = "\n"))
    outputs <- c(outputs, save_all(p_km, paste0(args$prefix, "_km")))
  }

  if ("forest" %in% types) {
    fx <- res[res$model %in% c("cox_univariate", "cox_adjusted") & res$status == "ok", , drop = FALSE]
    fx$term_display <- ifelse(fx$term == ".z",
                              paste0(fx$target_gene, " (per SD)"), fx$term)
    fx$ann <- sprintf("HR=%.2f [%.2f-%.2f], q=%.2g",
                      fx$estimate, fx$ci_low, fx$ci_high, fx$p_adjust)
    p_f <- ggplot2::ggplot(fx, ggplot2::aes(x = estimate, y = term_display)) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low, xmax = ci_high),
                              height = 0.15) +
      ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
      ggplot2::facet_grid(model + endpoint ~ block, scales = "free_y") +
      ggplot2::scale_x_log10() +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) +
      ggplot2::labs(x = "Hazard ratio (log scale)", y = NULL,
                    title = "Cox proportional-hazards models")
    outputs <- c(outputs, save_all(p_f, paste0(args$prefix, "_forest")))
  }

  if (!length(outputs)) fail("no figure could be produced (check --types and inputs)")
  log_outputs(setNames(as.list(outputs), paste0("figure_", seq_along(outputs))))
}

run_main(main)
