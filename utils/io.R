# ============================================================================
# io.R -- atomic writes, checksums, readers, downloads
#
# Rules enforced here:
#   * all writes go through a temp sibling file + rename (never partial files)
#   * raw data is written once and never modified in place
#   * checksums are computed with whatever tool the system provides
# ============================================================================

ensure_dir <- function(path, label = NULL) {
  if (!dir.exists(path)) {
    ok <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!ok && !dir.exists(path)) {
      fail(paste0("cannot create directory", (if (!is.null(label)) paste0(" (", label, ")") else ""), ": ", path))
    }
  }
  invisible(normalizePath(path, winslash = "/"))
}

file_size_bytes <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  as.numeric(file.info(path)$size)
}

require_file <- function(path, label = "input") {
  if (!file.exists(path)) fail(paste0(label, " file does not exist: ", path))
  if (file_size_bytes(path) == 0) fail(paste0(label, " file is empty: ", path))
  invisible(path)
}

# ---- checksums -------------------------------------------------------------

file_sha256 <- function(path) {
  require_file(path, "checksum target")
  # 1) R package `digest` if available
  if (requireNamespace("digest", quietly = TRUE)) {
    h <- digest::digest(file = path, algo = "sha256")
    return(list(algorithm = "sha256", value = h))
  }
  # 2) system sha256sum / shasum / openssl
  candidates <- list(
    c("sha256sum", path),
    c("shasum", "-a", "256", path),
    c("openssl", "dgst", "-sha256", path)
  )
  for (cmd in candidates) {
    if (nzchar(Sys.which(cmd[[1]]))) {
      out <- suppressWarnings(
        tryCatch(system2(cmd[[1]], shQuote(cmd[-1]), stdout = TRUE, stderr = FALSE),
                 error = function(e) NULL))
      if (length(out)) {
        m <- regmatches(paste(out, collapse = " "),
                        regexpr("[0-9a-fA-F]{64}", paste(out, collapse = " ")))
        if (length(m) == 1L && nzchar(m)) {
          return(list(algorithm = "sha256", value = tolower(m)))
        }
      }
    }
  }
  log_warn("no sha256 tool available; falling back to md5 (tools::md5sum)")
  list(algorithm = "md5", value = unname(tools::md5sum(path)))
}

# ---- atomic writers ---------------------------------------------------------

.atomic_write <- function(writer, path) {
  dir_path <- dirname(path)
  ensure_dir(dir_path)
  tmp <- tempfile(pattern = paste0(".", basename(path), "."),
                  tmpdir = dir_path, fileext = ".tmp")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  writer(tmp)
  if (file.exists(path)) {
    ok <- suppressWarnings(file.rename(tmp, path))
    if (!ok) { file.remove(path); ok <- file.rename(tmp, path) }
  } else {
    ok <- file.rename(tmp, path)
  }
  if (!ok) fail(paste0("atomic rename failed for: ", path))
  invisible(normalizePath(path, winslash = "/"))
}

write_tsv_atomic <- function(df, path, na = "NA") {
  writer <- function(tmp) write.table(df, file = tmp, sep = "\t", quote = FALSE,
                                      row.names = FALSE, col.names = TRUE,
                                      na = na, eol = "\n")
  .atomic_write(writer, path)
}

write_rds_atomic <- function(obj, path) {
  writer <- function(tmp) saveRDS(obj, file = tmp, version = 3)
  .atomic_write(writer, path)
}

write_yaml_atomic <- function(obj, path) {
  require_pkg("yaml", "YAML config support")
  txt <- yaml::as.yaml(obj, indent.mapping.sequence = TRUE)
  writer <- function(tmp) writeLines(enc2utf8(txt), con = tmp, useBytes = TRUE)
  .atomic_write(writer, path)
}

write_text_atomic <- function(txt, path) {
  writer <- function(tmp) writeLines(enc2utf8(txt), con = tmp, useBytes = TRUE)
  .atomic_write(writer, path)
}

# ---- readers ----------------------------------------------------------------

read_rds_input <- function(path, label = "RDS") {
  require_file(path, label)
  tryCatch(readRDS(path), error = function(e)
    fail(paste0("failed to read ", label, " file: ", path, " (", conditionMessage(e), ")")))
}

read_tsv <- function(path, header = TRUE, label = "TSV") {
  require_file(path, label)
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(path, sep = "\t", header = header,
                            stringsAsFactors = FALSE, data.table = FALSE,
                            showProgress = FALSE, check.names = FALSE,
                            na.strings = c("NA", "", "#N/A"))
    return(dt)
  }
  read.delim(path, sep = "\t", header = header, check.names = FALSE,
             stringsAsFactors = FALSE, quote = "", comment.char = "",
             na.strings = c("NA", "", "#N/A"))
}

# Read a gzipped TSV keeping only selected columns (by name).
# Uses data.table when present; otherwise base R with colClasses pruning.
read_gz_tsv_columns <- function(path, keep_cols, first_col_name = NULL) {
  require_file(path, "expression matrix")
  if (is.null(first_col_name)) {
    con <- gzfile(path, "rt"); on.exit(close(con))
    hdr_line <- readLines(con, n = 1L)
    first_col_name <- strsplit(hdr_line, "\t", fixed = TRUE)[[1]][[1]]
  }
  select_cols <- unique(c(first_col_name, keep_cols))
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(path, sep = "\t", header = TRUE,
                            select = select_cols, stringsAsFactors = FALSE,
                            data.table = FALSE, showProgress = FALSE,
                            check.names = FALSE,
                            na.strings = c("NA", "", "NaN"))
    return(dt)
  }
  log_warn("data.table not installed; reading large TSV with base R (slow)")
  con <- gzfile(path, "rt"); hdr <- strsplit(readLines(con, n = 1L), "\t", fixed = TRUE)[[1]]; close(con)
  cc <- rep("NULL", length(hdr)); cc[hdr %in% select_cols] <- NA
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             colClasses = cc, quote = "", comment.char = "",
             stringsAsFactors = FALSE, na.strings = c("NA", "", "NaN"))
}

# ---- download ---------------------------------------------------------------

download_file <- function(url, dest, timeout_s = 3600, attempts = 3) {
  ensure_dir(dirname(dest))
  tmp <- paste0(dest, ".part")
  old_timeout <- getOption("timeout")
  options(timeout = timeout_s)
  on.exit(options(timeout = old_timeout), add = TRUE)

  methods <- if (capabilities("libcurl")) c(getOption("download.file.method"), "curl") else c("curl", "wget")

  for (attempt in seq_len(attempts)) {
    method <- methods[[min(attempt, length(methods))]]
    log_info(sprintf("downloading (attempt %d/%d, method=%s): %s",
                     attempt, attempts, method, url))
    unlink(tmp)
    status <- tryCatch({
      if (method %in% c("curl", "wget")) {
        args <- if (method == "curl") c("-L", "--fail", "-o", shQuote(tmp), shQuote(url)) else c("-qO", shQuote(tmp), shQuote(url))
        system_result <- suppressWarnings(tryCatch(system2(method, args, stdout = FALSE, stderr = FALSE),
                                                   error = function(e) 1L))
        system_result
      } else {
        suppressWarnings(download.file(url, destfile = tmp, mode = "wb", quiet = TRUE))
        0L
      }
    }, error = function(e) 1L)
    if (identical(status, 0L) && file.exists(tmp) && file_size_bytes(tmp) > 0) {
      if (file.exists(dest)) file.remove(dest)
      ok <- file.rename(tmp, dest)
      if (!ok) { file.copy(tmp, dest, overwrite = TRUE); unlink(tmp) }
      log_info(sprintf("downloaded %s -> %s (%d bytes)", url, basename(dest),
                       file_size_bytes(dest)))
      return(invisible(file_size_bytes(dest)))
    }
    log_warn(sprintf("attempt %d failed for %s", attempt, url))
    Sys.sleep(min(30, 2^attempt))
  }
  unlink(tmp)
  fail(paste0("download failed after ", attempts, " attempts: ", url))
}
