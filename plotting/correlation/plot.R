#!/usr/bin/env Rscript
# ============================================================================
# Correlation plotting -- result table -> figures (NO statistics recomputation)
#
# Usage:
#   Rscript plotting/correlation/plot.R \
#       --input  results/<name>/tables/correlation_results.tsv \
#       [--plot-config config/plotting/default.yaml] \
#       [--out-dir results/<name>/figures] [--prefix correlation_heatmap]
#
# Reads the long-format correlation_results.tsv produced by analyses/correlation.
# Changing plot parameters never re-triggers download/prepare/analysis.
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
    input        = list(help = "correlation_results.tsv from the analysis module",
                        takes_value = TRUE, required = TRUE, default = NULL),
    `plot-config` = list(help = "plotting YAML", takes_value = TRUE,
                         required = FALSE, default = "config/plotting/default.yaml"),
    `out-dir`    = list(help = "figure output directory",
                        takes_value = TRUE, required = FALSE, default = NULL),
    prefix       = list(help = "output filename prefix", takes_value = TRUE,
                        required = FALSE, default = "correlation_heatmap")
  )
  args <- cli_parse("Plot correlation heatmaps from a results table.", spec)

  require_pkg("ggplot2", "figure rendering")

  input_path <- resolve_path(args$input)
  out_dir <- if (!is.null(args[["out-dir"]])) resolve_path(args[["out-dir"]])
             else file.path(dirname(dirname(input_path)), "figures")
  ensure_dir(out_dir)
  log_init(file.path(out_dir, "plot.log"))

  pcfg_path <- resolve_path(args[["plot-config"]])
  pcfg <- read_yaml_file(pcfg_path, "plot config")

  log_section("Correlation plot")
  log_inputs(list(results_table = input_path, plot_config = pcfg_path))
  log_params(list(prefix = args$prefix,
                  formats = as.character(get_field(pcfg, "figure.formats", c("pdf", "png"))),
                  out_dir = out_dir))

  res <- read_tsv(input_path, label = "correlation results")
  need_cols <- c("dataset", "cohort", "target_gene", "comparison_gene",
                 "correlation", "p_adjust")
  missing_cols <- setdiff(need_cols, colnames(res))
  if (length(missing_cols)) {
    fail(paste0("results table lacks column(s): ", paste(missing_cols, collapse = ", ")))
  }
  res$correlation <- as.numeric(res$correlation)
  res$p_adjust <- suppressWarnings(as.numeric(res$p_adjust))
  res$block <- paste(res$dataset, res$cohort, sep = " / ")

  # significance marks purely for display; statistics are NOT recomputed
  marks_cfg <- get_field(pcfg, "significance_marks", list())
  mark_of <- function(p) {
    if (is.na(p)) return("")
    for (m in marks_cfg) {
      if (!is.null(m$threshold) && !is.na(p) && p <= m$threshold) return(as.character(m$mark))
    }
    ""
  }
  res$mark <- vapply(res$p_adjust, mark_of, character(1))

  digits <- as.integer(get_field(pcfg, "heatmap.value_digits", 2))
  annotate <- isTRUE(get_field(pcfg, "heatmap.annotate_values", TRUE))
  res$label <- if (annotate) {
    ifelse(is.na(res$correlation), "NA",
           paste0(formatC(res$correlation, format = "f", digits = digits), res$mark))
  } else res$mark

  limits <- as.numeric(get_field(pcfg, "heatmap.limits", c(-1, 1)))
  midpoint <- as.numeric(get_field(pcfg, "heatmap.midpoint", 0))
  base_size <- as.numeric(get_field(pcfg, "heatmap.font_size", 10))

  p <- ggplot2::ggplot(res, ggplot2::aes(x = comparison_gene, y = target_gene,
                                         fill = correlation)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = base_size * 0.30,
                       na.rm = TRUE) +
    ggplot2::facet_wrap(~block, scales = "free") +
    ggplot2::scale_fill_gradient2(low = get_field(pcfg, "heatmap.color_low", "#2166AC"),
                                  mid = get_field(pcfg, "heatmap.color_mid", "#F7F7F7"),
                                  high = get_field(pcfg, "heatmap.color_high", "#B2182B"),
                                  na.value = get_field(pcfg, "heatmap.na_color", "#BEBEBE"),
                                  midpoint = midpoint,
                                  limits = limits,
                                  name = paste0(method_label(res), "\nrho")) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   strip.text = ggplot2::element_text(face = "bold"),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Target vs comparison gene correlation",
                  subtitle = sprintf("marks: %s",
                                     paste(vapply(marks_cfg, function(m)
                                       paste0(m$mark, " p<=", m$threshold), character(1)),
                                       collapse = ", ")))

  formats <- tolower(as.character(get_field(pcfg, "figure.formats", c("pdf", "png"))))
  dpi <- as.integer(get_field(pcfg, "figure.dpi", 600))
  width_in <- as.numeric(get_field(pcfg, "figure.width_in", 9))
  height_in <- as.numeric(get_field(pcfg, "figure.height_in", 5))

  outputs <- c()
  for (fmt in formats) {
    if (!fmt %in% c("pdf", "png", "svg")) fail(paste0("unsupported figure format: ", fmt))
    path <- file.path(out_dir, paste0(args$prefix, ".", fmt))
    ok <- tryCatch({ ggplot2::ggsave(path, p, width = width_in, height = height_in,
                                     dpi = dpi); TRUE },
                   error = function(e) { log_warn(paste0("failed saving ", fmt, ": ",
                                                         conditionMessage(e))); FALSE })
    if (ok) outputs[[length(outputs) + 1L]] <- path
  }
  if (!length(outputs)) fail("no figure could be written")
  log_outputs(setNames(as.list(outputs), paste0("figure_", seq_along(outputs))))
}

method_label <- function(df) {
  m <- unique(df$method)
  paste(unique(m[!is.na(m)]), collapse = "/")
}

run_main(main)
