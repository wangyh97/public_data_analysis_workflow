# ============================================================================
# logging.R -- structured run logging with loud failures
#
# Every script logs: start time, inputs, config summary, major parameters,
# sample/feature counts, warnings, outputs, and a final completion status.
# Severe errors must stop execution (fail loudly), never continue silently.
# ============================================================================

.pda_log_state <- new.env(parent = emptyenv())
.pda_log_state$file <- NULL
.pda_log_state$n_warnings <- 0L
.pda_log_state$start <- Sys.time()

pda_timestamp_utc <- function() {
  format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ")
}

log_init <- function(log_file = NULL) {
  .pda_log_state$start <- Sys.time()
  .pda_log_state$n_warnings <- 0L
  if (!is.null(log_file)) {
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    .pda_log_state$file <- log_file
    cat(paste0("# log started: ", pda_timestamp_utc(), "\n"),
        file = log_file, append = TRUE)
  }
  log_line("INFO", paste("R version:", R.version.string))
}

log_line <- function(level, msg) {
  line <- paste0("[", pda_timestamp_utc(), "] [", level, "] ", msg)
  message(line)
  if (!is.null(.pda_log_state$file)) {
    try(cat(line, "\n", sep = "", file = .pda_log_state$file, append = TRUE),
        silent = TRUE)
  }
  invisible(line)
}

log_info <- function(...) log_line("INFO", paste0(...))
log_warn <- function(...) {
  .pda_log_state$n_warnings <- .pda_log_state$n_warnings + 1L
  log_line("WARN", paste0(...))
}
log_error <- function(...) log_line("ERROR", paste0(...))

# Fail loudly: logs the error and raises; the top-level wrapper turns it into exit code 1.
fail <- function(msg) {
  log_error(msg)
  stop(msg, call. = FALSE)
}

log_section <- function(title) {
  log_line("INFO", "----------------------------------------------------------------")
  log_line("INFO", title)
}

log_kv <- function(items) {
  for (nm in names(items)) {
    v <- items[[nm]]
    if (is.null(v)) {
      txt <- "<null>"
    } else if (length(v) > 8) {
      txt <- paste0(paste(format(head(v, 8)), collapse = ", "),
                    ", ... (n=", length(v), ")")
    } else {
      txt <- paste(format(v), collapse = ", ")
    }
    log_line("INFO", sprintf("  %-28s %s", paste0(nm, ":"), txt))
  }
}

log_inputs <- function(inputs) {
  log_section("Inputs")
  log_kv(inputs)
}

log_params <- function(params) {
  log_section("Major parameters")
  log_kv(params)
}

log_outputs <- function(outputs) {
  log_section("Outputs")
  log_kv(outputs)
}

log_finalize <- function(status) {
  elapsed <- as.numeric(difftime(Sys.time(), .pda_log_state$start, units = "secs"))
  log_line("INFO", sprintf("warnings emitted: %d | elapsed: %.1f s",
                           .pda_log_state$n_warnings, elapsed))
  log_line("INFO", paste("STATUS:", status))
}

# Top-level wrapper for every script's main(). Converts errors into exit codes.
run_main <- function(fn) {
  tryCatch({
    fn()
    log_finalize("SUCCESS")
    .pda_exit(0)
  }, error = function(e) {
    log_line("ERROR", conditionMessage(e))
    log_finalize("FAILED")
    .pda_exit(1)
  }, warning = function(w) {
    # warnings raised as conditions via `stop(warningCondition(...))` land here
    log_line("ERROR", conditionMessage(w))
    log_finalize("FAILED")
    .pda_exit(1)
  })
}
