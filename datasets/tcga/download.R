#!/usr/bin/env Rscript
# ============================================================================
# TCGA download -- fetch raw public files once, reuse forever, never modify.
#
# Usage:
#   Rscript datasets/tcga/download.R --config config/datasets/tcga.yaml [--force]
#
# Behaviour:
#   * skips files that already exist (unless --force)
#   * verifies existing files against manifest checksums; mismatch => loud error
#   * writes/updates raw_root/manifest.yaml with source URL, date, bytes, sha256
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
    config = list(help = "path to dataset YAML", takes_value = TRUE,
                  required = TRUE, default = NULL),
    force  = list(help = "re-download even if the file already exists",
                  takes_value = FALSE, required = FALSE, default = FALSE)
  )
  args <- cli_parse("Download raw TCGA (UCSC Xena Toil) files.", spec)

  cfg <- read_yaml_file(resolve_path(args$config), "dataset")
  require_fields(cfg, c("paths.raw_root", "download.files"), context = "tcga.yaml")

  raw_root <- resolve_path(cfg$paths$raw_root)
  ensure_dir(raw_root, "raw_root")
  log_init(file.path(raw_root, "logs", "download.log"))

  log_section("TCGA download")
  log_inputs(list(config = normalizePath(args$config, winslash = "/"),
                  raw_root = raw_root,
                  force = isTRUE(args$force)))
  log_params(list(files = names(cfg$download$files),
                  http_timeout_s = get_field(cfg, "download.http_timeout_seconds", 3600)))

  manifest_path <- file.path(raw_root, "manifest.yaml")
  old_manifest <- list()
  if (file.exists(manifest_path)) {
    old_manifest <- tryCatch(read_yaml_file(manifest_path, "existing manifest"),
                             error = function(e) list())
    old_entries <- get_field(old_manifest, "files", list())
  } else old_entries <- list()

  entries <- list()
  for (key in names(cfg$download$files)) {
    f <- cfg$download$files[[key]]
    dest <- file.path(raw_root, f$filename)
    log_section(paste0("file: ", key))

    entry <- list(url = f$url, filename = f$filename,
                  description = get_field(f, "description", ""),
                  downloaded_at = NULL, bytes = NA_real_,
                  sha256 = NA_character_, sha_algorithm = NA_character_, status = NA_character_)

    declared <- get_field(f, "approx_bytes", NULL)

    if (file.exists(dest) && !isTRUE(args$force)) {
      size <- file_size_bytes(dest)
      chk <- file_sha256(dest)
      prev <- old_entries[[key]]
      if (!is.null(prev) && !is.null(prev$sha256) &&
          identical(tolower(prev$sha256), chk$value) && !identical(chk$value, "NA")) {
        log_info(sprintf("reuse (checksum verified): %s (%d bytes)", dest, size))
      } else if (!is.null(prev) && !is.null(prev$sha256)) {
        fail(paste0("existing file checksum differs from manifest record for ", key,
                    " (", basename(dest), "). The upstream file may have changed or the ",
                    "local copy was altered. Re-run with --force to re-download deliberately."))
      } else {
        log_warn(sprintf("file exists but has no prior manifest record; adopting it: %s", dest))
      }
      entry$status <- "already_present"
      entry$bytes <- size
      entry$sha256 <- tolower(chk$value)
      entry$sha_algorithm <- chk$algorithm
    } else {
      if (isTRUE(args$force)) log_info("--force given; re-downloading")
      timeout <- as.numeric(get_field(cfg, "download.http_timeout_seconds", 3600))
      download_file(f$url, dest, timeout_s = timeout)
      if (declared && file_size_bytes(dest) < declared * 0.5) {
        fail(paste0("downloaded file is much smaller than expected for ", key,
                    " (got ", file_size_bytes(dest), " bytes, expected ~", declared, ")"))
      }
      chk <- file_sha256(dest)
      entry$status <- "downloaded"
      entry$downloaded_at <- Sys.Date()
      entry$bytes <- file_size_bytes(dest)
      entry$sha256 <- tolower(chk$value)
      entry$sha_algorithm <- chk$algorithm
    }
    if (declared) entry$declared_bytes <- declared
    entries[[key]] <- entry
  }

  manifest <- list(
    dataset = get_field(cfg, "dataset.source", "UNKNOWN"),
    pipeline_version = PDA_PIPELINE_VERSION,
    updated_at_utc = pda_timestamp_utc(),
    files = entries
  )
  write_yaml_atomic(manifest, manifest_path)

  log_outputs(list(raw_root = raw_root, manifest = manifest_path))
  for (key in names(entries)) {
    log_info(sprintf("  %-16s %-14s %12s bytes", key, entries[[key]]$status,
                     format(entries[[key]]$bytes, big.mark = ",")))
  }
}

run_main(main)
