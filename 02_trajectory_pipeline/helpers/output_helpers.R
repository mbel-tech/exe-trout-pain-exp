# GOT pipeline — shared output-management helpers
# Sourced by 00_master_pipeline_GOT.R and defensively by each step script.

# =============================================================================
# resolve_dataset_path()
#
# Resolves a (possibly machine-specific) file path with a fallback to
# data/ (relative to the repository root).
#
# Lookup order:
#   1. original path exists on disk → return it unchanged
#   2. search fallback_dir for a file with the same basename
#      - exactly one match → return it (with an informative message)
#      - no match          → stop with a clear error listing both locations
#      - multiple matches  → stop listing all candidates (user must disambiguate)
# =============================================================================
resolve_dataset_path <- function(path,
                                 fallback_dir = file.path(get0(".repo_root", ifnotfound = ".."), "data")) {
  if (is.null(path) || is.na(path) || !nzchar(trimws(path))) return(path)
  if (file.exists(path)) return(path)

  bname <- basename(path)

  if (!dir.exists(fallback_dir)) {
    stop(
      "File not found: ", path,
      "\nFallback directory does not exist: ", fallback_dir,
      "\nCreate the directory or fix the path in the pipeline config.",
      call. = FALSE
    )
  }

  candidates <- list.files(fallback_dir, full.names = TRUE, recursive = FALSE)
  matches    <- candidates[basename(candidates) == bname]

  if (length(matches) == 1L) {
    message("resolve_dataset_path: '", bname,
            "' not at original location — using fallback:\n  ", matches[1L])
    return(matches[1L])
  }

  if (length(matches) == 0L) {
    stop(
      "File not found: ", path,
      "\nAlso searched fallback dir: ", fallback_dir,
      "\nNo file named '", bname, "' found in either location.",
      "\nCopy the file to '", fallback_dir, "' or fix the path in the pipeline config.",
      call. = FALSE
    )
  }

  stop(
    "File not found at: ", path,
    "\nMultiple files named '", bname, "' in ", fallback_dir, ":\n",
    paste0("  ", matches, collapse = "\n"),
    "\nRemove duplicates or set the correct path explicitly.",
    call. = FALSE
  )
}

.got_make_output_dir <- function(step_name) {
  parent  <- file.path(getwd(), step_name)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  run_dir <- file.path(parent, paste0(step_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  run_dir
}

.got_find_latest_output <- function(step_name, filename = NULL) {
  parent  <- file.path(getwd(), step_name)
  if (!dir.exists(parent)) return(NULL)
  subdirs <- list.dirs(parent, full.names = TRUE, recursive = FALSE)
  pattern <- paste0("^", step_name, "_\\d{8}_\\d{6}$")
  subdirs <- subdirs[grepl(pattern, basename(subdirs))]
  if (length(subdirs) == 0L) return(NULL)
  latest  <- subdirs[which.max(file.mtime(subdirs))]
  if (is.null(filename)) return(latest)
  candidate <- file.path(latest, filename)
  if (file.exists(candidate)) candidate else NULL
}

.got_write_session_info <- function(dir) {
  tryCatch({
    f <- file.path(dir, "sessionInfo.txt")
    sink(f); on.exit(sink(), add = TRUE)
    cat("sessionInfo written:", format(Sys.time()), "\n\n")
    print(sessionInfo())
  }, error = function(e) warning("sessionInfo write failed: ", conditionMessage(e), call. = FALSE))
}

.got_write_run_settings <- function(dir, ...) {
  settings <- c(list(...),
                list(run_timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
                     r_version     = paste(R.version$major, R.version$minor, sep = ".")))
  df <- data.frame(
    parameter = names(settings),
    value     = vapply(settings, function(x) paste(x, collapse = ","), character(1L)),
    stringsAsFactors = FALSE
  )
  write.csv(df, file.path(dir, "run_settings.csv"), row.names = FALSE)
}

.got_append_pipeline_log <- function(step_number, step_name, out_dir, elapsed_s) {
  log_file <- file.path(getwd(), "GOT_pipeline_runs.csv")
  entry <- data.frame(
    timestamp  = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    step       = step_number,
    step_name  = step_name,
    output_dir = out_dir,
    elapsed_s  = round(elapsed_s, 1),
    r_version  = paste(R.version$major, R.version$minor, sep = "."),
    stringsAsFactors = FALSE
  )
  if (file.exists(log_file)) {
    existing <- tryCatch(read.csv(log_file, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(existing)) entry <- rbind(existing, entry)
  }
  write.csv(entry, log_file, row.names = FALSE)
}
