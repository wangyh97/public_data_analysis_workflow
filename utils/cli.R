# ============================================================================
# cli.R -- minimal dependency-free command line argument parser
#
# Every executable script in this project is a CLI tool. This parser supports:
#   --key value   --key=value   --flag
# and auto-adds --help / -h.
#
# Exit codes (project-wide convention):
#   0  success
#   1  runtime error
#   2  usage / configuration error
# ============================================================================

.pda_exit <- function(status) {
  quit(save = "no", status = status)
}

cli_fail <- function(msg) {
  message("USAGE ERROR: ", msg)
  message("Run with --help for usage information.")
  .pda_exit(2)
}

.cli_print_help <- function(description, spec) {
  cat("================================================================\n")
  if (!is.null(description)) cat(description, "\n\n", sep = "")
  cat("Options:\n")
  keys <- names(spec)
  for (k in keys) {
    s <- spec[[k]]
    takes <- if (isTRUE(s$takes_value)) " VALUE" else ""
    req <- if (isTRUE(s$required)) " [required]" else ""
    def <- if (!is.null(s$default)) paste0(" [default: ", as.character(s$default), "]") else ""
    cat(sprintf("  --%s%s%s%s\n      %s\n", k, takes, req, def,
                ifelse(is.null(s$help), "", s$help)))
  }
}

cli_parse <- function(description, spec) {
  spec[["help"]] <- list(help = "show this help and exit",
                         takes_value = FALSE, required = FALSE, default = NULL)
  raw <- commandArgs(trailingOnly = TRUE)

  if ("--help" %in% raw || "-h" %in% raw || length(raw) == 0 && FALSE) {
    .cli_print_help(description, spec)
    .pda_exit(0)
  }
  # tolerate bare invocation with no args: show help instead of failing cryptically
  if (length(raw) == 0) {
    .cli_print_help(description, spec)
    .pda_exit(2)
  }

  values <- list()
  i <- 1L
  while (i <= length(raw)) {
    tok <- raw[[i]]
    if (!grepl("^--[A-Za-z][A-Za-z0-9_-]*$", tok)) {
      cli_fail(paste0("unrecognized argument token: ", shQuote(tok)))
    }
    key <- substring(tok, 3L)
    if (!key %in% names(spec)) {
      cli_fail(paste0("unknown option --", key))
    }
    s <- spec[[key]]
    if (isTRUE(s$takes_value)) {
      if (grepl("=", tok, fixed = TRUE)) {
        values[[key]] <- sub("^[^=]*=", "", tok)
        i <- i + 1L
      } else {
        if (i + 1L > length(raw)) {
          cli_fail(paste0("option --", key, " expects a value"))
        }
        values[[key]] <- raw[[i + 1L]]
        i <- i + 2L
      }
    } else {
      if (grepl("=", tok, fixed = TRUE)) {
        cli_fail(paste0("option --", key, " is a flag and takes no value"))
      }
      values[[key]] <- TRUE
      i <- i + 1L
    }
  }

  for (key in names(spec)) {
    s <- spec[[key]]
    if (is.null(values[[key]]) && !is.null(s$default)) values[[key]] <- s$default
    if (isTRUE(s$required) && is.null(values[[key]])) {
      cli_fail(paste0("missing required option --", key))
    }
    if (identical(s$type, "integer") && !is.null(values[[key]])) {
      v <- suppressWarnings(as.integer(values[[key]]))
      if (is.na(v)) cli_fail(paste0("--", key, " must be an integer"))
      values[[key]] <- v
    }
  }
  values[["help"]] <- NULL
  values
}
