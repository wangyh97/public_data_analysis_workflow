# ============================================================================
# config.R -- YAML config loading, env-var path expansion, registry lookup,
#             gene-set resolution, centralized defaults
#
# Path resolution order (highest priority first):
#   CLI flag  >  config file value  >  environment variable
# Config strings may reference environment variables as ${VAR_NAME}.
# ============================================================================

require_pkg <- function(pkg, purpose = "") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    fail(paste0("required R package '", pkg, "' is not installed",
                (if (nzchar(purpose)) paste0(" (", purpose, ")") else ""),
                ". Install with: install.packages('", pkg, "')"))
  }
}

read_yaml_file <- function(path, label = "config") {
  require_pkg("yaml", "reading YAML configuration")
  require_file(path, label)
  cfg <- tryCatch(yaml::read_yaml(path), error = function(e)
    fail(paste0("invalid YAML in ", label, " file ", path, ": ", conditionMessage(e))))
  if (!is.list(cfg)) fail(paste0(label, " must deserialize to a mapping: ", path))
  cfg
}

# ---- ${ENV_VAR} expansion ---------------------------------------------------

.expand_env_string <- function(s) {
  matches <- gregexpr("\\$\\{[A-Za-z_][A-Za-z0-9_]*\\}", s)[[1]]
  if (matches[[1]] == -1L) return(s)
  tokens <- regmatches(s, matches)[[1]]
  for (tok in tokens) {
    var <- substr(tok, 3L, nchar(tok) - 1L)
    val <- Sys.getenv(var)
    if (!nzchar(val)) {
      fail(paste0("environment variable '", var, "' is referenced but not set",
                  " (string: ", s, "). Set it or pass the path explicitly."))
    }
    s <- gsub(tok, val, s, fixed = TRUE)
  }
  s
}

expand_env_recursive <- function(x) {
  if (is.character(x)) return(vapply(x, .expand_env_string, character(1)))
  if (is.list(x)) {
    for (nm in names(x)) x[[nm]] <- expand_env_recursive(x[[nm]])
    return(x)
  }
  x
}

resolve_path <- function(path, base = PROJECT_ROOT) {
  p <- .expand_env_string(path)
  p <- path.expand(p)
  if (!grepl("^([A-Za-z]:[/\\\\]|/|\\\\)", p)) p <- file.path(base, p)
  normalizePath(p, winslash = "/", mustWork = FALSE)
}

env_or_fail <- function(var, hint = NULL) {
  v <- Sys.getenv(var)
  if (!nzchar(v)) {
    fail(paste0("environment variable '", var, "' is not set.",
                (if (!is.null(hint)) paste0(" ", hint) else "")))
  }
  v
}

data_root <- function() env_or_fail(
  "PDA_DATA_ROOT",
  hint = "It should point to the shared public-data directory, e.g. /data/public_data"
)

results_root <- function() env_or_fail(
  "PDA_RESULTS_ROOT",
  hint = "It should point to the analysis output directory, e.g. /project/results"
)

# ---- deep merge / field access ----------------------------------------------

deep_merge <- function(defaults, override) {
  if (is.null(defaults)) return(override)
  if (is.null(override)) return(defaults)
  out <- defaults
  for (nm in names(override)) {
    d <- defaults[[nm]]; o <- override[[nm]]
    if (is.list(d) && is.list(o) && !is.data.frame(d) && !is.data.frame(o)) {
      out[[nm]] <- deep_merge(d, o)
    } else {
      out[[nm]] <- o
    }
  }
  out
}

get_field <- function(cfg, dotted_path, default = NULL) {
  parts <- strsplit(dotted_path, ".", fixed = TRUE)[[1]]
  cur <- cfg
  for (p in parts) {
    if (!is.list(cur) || !p %in% names(cur)) return(default)
    cur <- cur[[p]]
  }
  if (is.null(cur)) default else cur
}

require_fields <- function(cfg, dotted_paths, context = "config") {
  for (p in dotted_paths) {
    if (is.null(get_field(cfg, p))) {
      fail(paste0(context, ": missing required field '", p, "'"))
    }
  }
  invisible(TRUE)
}

# ---- dataset registry --------------------------------------------------------

registry_path <- function() file.path(PROJECT_ROOT, "config", "dataset_registry.yaml")

load_registry <- function() {
  reg <- read_yaml_file(registry_path(), "dataset registry")
  if (is.null(get_field(reg, "datasets"))) fail("dataset registry has no 'datasets' section")
  reg
}

lookup_dataset <- function(registry, source) {
  entry <- get_field(registry, paste0("datasets.", source))
  if (is.null(entry)) {
    fail(paste0("dataset '", source, "' is not registered in ",
                registry_path(), ". Register it before use."))
  }
  status <- tolower(as.character(get_field(entry, "status", "active")))
  if (status %in% c("planned", "disabled")) {
    fail(paste0("dataset adapter '", source, "' exists in the registry but its status is '",
                status, "'; it is not usable yet"))
  }
  entry
}

processed_dir_for <- function(source, cohort, registry, data_root_value = NULL) {
  dr <- if (is.null(data_root_value)) data_root() else data_root_value
  layout <- get_field(registry, paste0("datasets.", source, ".processed_layout"),
                      get_field(registry, "defaults.processed_layout",
                                "{processed_root}/{source}/{cohort}"))
  out <- gsub("{processed_root}", file.path(dr, "processed"), layout, fixed = TRUE)
  out <- gsub("{source}", source, out, fixed = TRUE)
  out <- gsub("{cohort}", cohort, out, fixed = TRUE)
  normalizePath(out, winslash = "/", mustWork = FALSE)
}

# ---- gene lists / gene sets ---------------------------------------------------

resolve_gene_spec <- function(spec, base_dir = PROJECT_ROOT, role = "genes") {
  # spec may be: list(genes=c(...)) and/or list(geneset_file="path"); both combinable.
  genes <- character(0)
  inline <- spec$genes
  if (!is.null(inline)) genes <- c(genes, as.character(inline))
  gs_file <- spec$geneset_file
  if (!is.null(gs_file)) {
    path <- resolve_path(gs_file, base_dir)
    require_file(path, paste0(role, " gene set"))
    lines <- readLines(path, warn = FALSE)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    genes <- c(genes, lines)
    log_info(sprintf("%s: loaded %d genes from gene set %s", role, length(lines), basename(path)))
  }
  genes <- unique(genes[nzchar(genes)])
  if (!length(genes)) fail(paste0(role, " specification resolved to an empty gene list"))
  genes
}

# Save the fully-resolved effective config next to results (provenance).
save_config_used <- function(cfg, path) write_yaml_atomic(cfg, path)
