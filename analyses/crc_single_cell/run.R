#!/usr/bin/env Rscript
# Workflow entry point for the validated CRC single-cell panels.

PROJECT_ROOT <- local({
  fa <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  d <- dirname(normalizePath(sub("^--file=", "", fa[[1]]), winslash = "/"))
  repeat { if (file.exists(file.path(d, "utils", "common.R"))) break; p <- dirname(d); if (identical(p,d)) stop("project root not found"); d <- p }
  d
})
source(file.path(PROJECT_ROOT, "utils", "common.R"))

copy_tree <- function(src, dst) {
  if (!dir.exists(src)) fail(paste0("missing panel path: ", src))
  ensure_dir(dst)
  files <- list.files(src, recursive=TRUE, full.names=TRUE, include.dirs=FALSE)
  for (f in files) {
    rel <- substring(f, nchar(src)+2L); out <- file.path(dst, rel)
    ensure_dir(dirname(out)); file.copy(f, out, overwrite=TRUE)
  }
}

main <- function() {
  args <- cli_parse("Run CRC single-cell panels through public_data_analysis_workflow.", list(
    config=list(help="CRC analysis YAML",takes_value=TRUE,required=TRUE,default=NULL),
    `output-root`=list(help="override PDA_RESULTS_ROOT",takes_value=TRUE,required=FALSE,default=NULL),
    `panel-root`=list(help="override CRC_PANEL_ROOT",takes_value=TRUE,required=FALSE,default=NULL),
    `skip-processing`=list(help="reuse existing panel source files",takes_value=FALSE,required=FALSE,default=FALSE),
    `skip-plotting`=list(help="stage existing figures without redrawing",takes_value=FALSE,required=FALSE,default=FALSE)
  ))
  cfg <- read_yaml_file(resolve_path(args$config), "CRC analysis config")
  panel_root <- if (!is.null(args[["panel-root"]])) resolve_path(args[["panel-root"]]) else {
    x <- Sys.getenv("CRC_PANEL_ROOT", ""); if (!nzchar(x)) fail("set CRC_PANEL_ROOT or pass --panel-root"); resolve_path(x)
  }
  out_root <- if (!is.null(args[["output-root"]])) resolve_path(args[["output-root"]]) else resolve_path(Sys.getenv("PDA_RESULTS_ROOT", "results"))
  out <- file.path(out_root, get_field(cfg, "analysis.name", "crc_single_cell_panels")); ensure_dir(out)
  ensure_dir(file.path(out,"logs")); ensure_dir(file.path(out,"metadata"))
  log_init(file.path(out, "logs", "run.log")); log_section("CRC single-cell panels")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  panel_scripts <- file.path(panel_root, c("01_MHC_I_correlation/script", "02_T_cell_state_correlation/script", "03_tumor_treatment_trajectories/script"))
  for (s in panel_scripts) {
    proc <- file.path(s, "data_processing.R"); plot <- file.path(s, "plot.R")
    if (isTRUE(get_field(cfg,"execution.run_processing",TRUE)) && !isTRUE(args[["skip-processing"]])) { status <- system2(rscript, proc); if (status != 0) fail(paste0("processing failed: ", proc)) }
    if (isTRUE(get_field(cfg,"execution.run_plotting",TRUE)) && !isTRUE(args[["skip-plotting"]])) { status <- system2(rscript, plot); if (status != 0) fail(paste0("plotting failed: ", plot)) }
    key <- basename(dirname(dirname(s))); copy_tree(dirname(s), file.path(out,"panels",key))
  }
  save_config_used(cfg, file.path(out,"metadata/config_used.yaml"))
  writeLines(c("STATUS: SUCCESS", paste0("panel_root: ",panel_root), paste0("output: ",out)), file.path(out,"metadata/provenance.txt"))
  writeLines(capture.output(sessionInfo()), file.path(out,"sessionInfo.txt"))
  cat("CRC single-cell panels completed: ", out, "\n", sep="")
}
run_main(main)
