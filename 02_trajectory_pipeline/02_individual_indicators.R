# =============================================================================
# MASTER OUTPUTS PIPELINE (GOT) — ROBUST VERSION FOR "GROUP_* + MASTER_FISH_BY_FRAME"
#
# Supports by_frame subfolder layout like:
#   group_group/ , group_ind/ , group_pairwise/ , occupancy/ , position/ , speed/
#   + master_fish_by_frame.(csv|xlsx)
#
# Changes vs previous:
# - Canonical fish-by-frame source is now ONLY: master_fish_by_frame
#   Columns expected:
#     session_id, trial_id, segment_id, fish_id, fish_in_session, frame, time,
#     x_interp, y_interp, notes, zone, speed_cm_s, heading_rad
# - position_by_frame / speed_by_frame / zone_by_frame are DERIVED from master_fish_by_frame
# - group_by_frame is loaded from master_group_dynamics_by_frame.csv found under ANY "group_*" folder
#   (instead of master_groupmetrics_by_frame.csv).
# - Folder detection is robust: clipboard path may point to master_outputs root OR directly to by_frame.
# - Adds framewise speed in body lengths per second:
#     speed_BL_s = speed_cm_s / body_length_cm
#   and appends it to master_fish_by_frame in .GlobalEnv and on disk.
# =============================================================================

# =============================================================================
# ==== 0) PACKAGES =============================================================
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(ggplot2)
  library(patchwork)
  library(officer)
  library(rvg)
  library(data.table)
  library(tibble)
  library(purrr)
  library(readr)
})

# ---- output helpers (standalone guard) ----
if (!exists(".got_make_output_dir", mode = "function")) {
  tryCatch(
    source(file.path(getwd(), "02_trajectory_pipeline", "helpers", "output_helpers.R"), local = FALSE),
    error = function(e) tryCatch(
      source(file.path(getwd(), "helpers/output_helpers.R"), local = FALSE),
      error = function(e2) stop(
        "Cannot find 00_output_helpers.R. Run via master pipeline or set working directory appropriately.",
        call. = FALSE
      )
    )
  )
}


# =============================================================================
# ==== 1) GLOBAL CONSTANTS (NON-PLOT) =========================================
# =============================================================================

# ---- Input metadata ----
TRIAL_SUMMARY_PATH <- if (exists("TRIAL_SUMMARY_XLSX", envir = .GlobalEnv, inherits = FALSE)) {
  get("TRIAL_SUMMARY_XLSX", envir = .GlobalEnv)
} else {
  "data/trials summary GOT.xlsx"
}

# ---- Canonical load policy ----
FORCE_RELOAD_FROM_DISK <- FALSE

# ---- Calibration ----
LENGTH_UNIT_PX_PER_CM <- 8.94681

# ---- Activity budget settings ----
MIN_BOUT_DURATION <- 300
THRESHOLD_MULTIPLIERS <- c(`0.5` = 0.5, `1` = 1, `1.5` = 1.5, `2` = 2)

# ---- Common join keys ----
JOIN_KEY_FISH_SEGMENT <- c("session_id","trial_id","timepoint","fish_in_session")
JOIN_KEY_FISH_FRAME   <- c("session_id","trial_id","timepoint","fish_in_session","fish_id")


# =============================================================================
# ==== 2) UTILITIES ============================================================
# =============================================================================

ts_msg <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), " | ", ..., "\n", sep = ""); flush.connection(stdout()) }

# ---- Safe console printer for tibbles/data.frames ----------------------------
print_df_safe <- function(x, title = NULL) {
  if (!is.null(title)) {
    cat("\n", title, "\n", sep = "")
  }
  
  if (is.null(x)) {
    cat("<NULL>\n")
    return(invisible(NULL))
  }
  
  x_tbl <- tibble::as_tibble(x)
  print(x_tbl, n = Inf, width = Inf)
  invisible(x_tbl)
}

# ---- Validate canonical master_fish_by_frame structure ----
is_valid_master_fish_by_frame <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0) return(FALSE)
  
  required_base <- c(
    "session_id", "trial_id", "fish_id", "fish_in_session", "frame", "time",
    "x_interp", "y_interp", "notes", "speed_cm_s", "heading_rad"
  )
  
  has_time_col <- ("segment_id" %in% names(x)) || ("timepoint" %in% names(x))
  
  all(required_base %in% names(x)) && has_time_col
}

require_cols <- function(df, cols, df_name = deparse(substitute(df))) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) stop(df_name, " missing: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

# ---- Windows clipboard path -> normalized folder ----
clipboard_dir <- function() {
  if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
    p <- normalizePath(get("PIPELINE_DATA_DIR", envir = .GlobalEnv), winslash = "/", mustWork = FALSE)
    mo <- file.path(p, "master_outputs")
    if (dir.exists(mo)) return(mo)
    return(p)
  }
  default_dir <- "data/tracking_sessions"
  if (dir.exists(default_dir)) {
    mo <- file.path(default_dir, "master_outputs")
    if (dir.exists(mo)) {
      message("clipboard_dir: using default master_outputs directory: ", mo)
      return(normalizePath(mo, winslash = "/", mustWork = FALSE))
    }
    message("clipboard_dir: using default data directory: ", default_dir)
    return(normalizePath(default_dir, winslash = "/", mustWork = FALSE))
  }
  d <- trimws(readLines("clipboard", warn = FALSE))
  d <- d[nzchar(d)][1]
  if (is.na(d) || !nzchar(d)) stop("Clipboard is empty. Copy folder path first.", call. = FALSE)

  normalizePath(
    gsub("\\\\", "/", gsub('^"|"$', "", gsub("\r", "", d))),
    winslash = "/",
    mustWork = FALSE
  )
}

# ---- Fast CSV reader ----
read_csv_safe <- function(path) {
  tryCatch({
    as_tibble(data.table::fread(path, showProgress = FALSE))
  }, error = function(e1) {
    tryCatch(
      readr::read_csv(path, show_col_types = FALSE),
      error = function(e2) {
        warning("Failed to read CSV: ", path, "\n  ", conditionMessage(e2), call. = FALSE)
        NULL
      }
    )
  })
}

# ---- Safe console printer for tibbles/data.frames ----------------------------
print_df_safe <- function(x, title = NULL) {
  if (!is.null(title)) {
    cat("\n", title, "\n", sep = "")
  }
  
  if (is.null(x)) {
    cat("<NULL>\n")
    return(invisible(NULL))
  }
  
  x_tbl <- tibble::as_tibble(x)
  print(x_tbl, n = Inf, width = Inf)
  invisible(x_tbl)
}

# ---- Read "table" from CSV or Excel ----
read_table_safe <- function(path) {
  if (is.null(path) || !is.character(path) || length(path) != 1) return(NULL)
  if (!file.exists(path)) return(NULL)
  
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("csv", "txt")) return(read_csv_safe(path))
  if (ext %in% c("xlsx", "xls")) {
    return(tryCatch(
      as_tibble(readxl::read_excel(path)),
      error = function(e) {
        warning("Failed to read Excel: ", path, "\n  ", conditionMessage(e), call. = FALSE)
        NULL
      }
    ))
  }
  
  warning("Unsupported file extension: ", path, call. = FALSE)
  NULL
}

# ---- Timestamp helper ----
timestamp_tag <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

# ---- Detect by_frame directory robustly ----
detect_by_frame_dir <- function(root_dir) {
  stopifnot(dir.exists(root_dir))
  
  cand1 <- file.path(root_dir, "by_frame")
  
  subdirs <- list.dirs(root_dir, full.names = FALSE, recursive = FALSE)
  has_group <- length(subdirs[grepl("^group_", subdirs)]) > 0
  has_pos   <- dir.exists(file.path(root_dir, "position"))
  has_speed <- dir.exists(file.path(root_dir, "speed"))
  has_master_fish <- length(list.files(
    root_dir,
    pattern = "^master_fish_by_frame\\.(csv|xlsx|xls)$",
    ignore.case = TRUE,
    full.names = TRUE
  )) > 0
  
  if (dir.exists(cand1)) return(cand1)
  if (has_group || has_pos || has_speed || has_master_fish) return(root_dir)
  
  stop(
    "Could not locate by_frame. Clipboard path must be either:\n",
    " - the master_outputs root containing a 'by_frame' folder, OR\n",
    " - the by_frame folder itself (with group_*/position/speed/master_fish_by_frame.*).",
    call. = FALSE
  )
}

# ---- Find a file by basename, recursively ----
find_file_recursive <- function(root_dir, basename_regex, must_contain_regex = NULL) {
  stopifnot(dir.exists(root_dir))
  
  all <- list.files(
    root_dir,
    pattern = basename_regex,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (!is.null(must_contain_regex)) {
    all <- all[grepl(must_contain_regex, all, ignore.case = TRUE)]
  }
  
  if (length(all) == 0) return(NULL)
  
  all_norm <- normalizePath(all, winslash = "/", mustWork = FALSE)
  depth <- stringr::str_count(all_norm, "/")
  
  fi <- file.info(all)
  mtime_num <- as.numeric(fi$mtime)
  mtime_num[is.na(mtime_num)] <- -Inf
  
  ord <- order(depth, -mtime_num, na.last = TRUE)
  all[ord][1]
}

# ---- Prefer object from .GlobalEnv unless forced reload ----
get_global_or_load <- function(object_name, loader_expr, validate_fn = NULL, force_reload = FALSE) {
  if (!isTRUE(force_reload) && exists(object_name, envir = .GlobalEnv, inherits = FALSE)) {
    obj <- get(object_name, envir = .GlobalEnv, inherits = FALSE)
    ts_msg("Using existing object from .GlobalEnv: ", object_name)
    
    if (!is.null(validate_fn)) {
      ok <- tryCatch(validate_fn(obj), error = function(e) FALSE)
      if (!isTRUE(ok)) {
        warning(
          "Object '", object_name, "' found in .GlobalEnv but failed validation; reloading from disk.",
          call. = FALSE
        )
      } else {
        return(obj)
      }
    } else {
      return(obj)
    }
  }
  
  if (isTRUE(force_reload)) {
    ts_msg("FORCE_RELOAD_FROM_DISK = TRUE -> reloading object from disk/build step: ", object_name)
  } else {
    ts_msg("Object not found in .GlobalEnv (or failed validation): ", object_name, " -> loading from disk/build step")
  }
  
  loader_expr
}

# ---- Resolve by_frame path on demand from clipboard ----
resolve_by_frame_from_clipboard <- function() {
  root_dir <- clipboard_dir()
  ts_msg("Clipboard path: ", root_dir)
  stopifnot(dir.exists(root_dir))
  
  by_frame_dir <- detect_by_frame_dir(root_dir)
  ts_msg("Detected by_frame_dir: ", by_frame_dir)
  
  master_root_guess <- if (basename(by_frame_dir) == "by_frame") dirname(by_frame_dir) else dirname(by_frame_dir)
  morph_dir <- file.path(master_root_guess, "morph")
  
  list(
    root_dir = root_dir,
    by_frame_dir = by_frame_dir,
    morph_dir = morph_dir
  )
}


# =============================================================================
# ==== 2b) STEP OUTPUT DIRECTORY ==============================================
# =============================================================================

STEP_OUTPUT_DIR <- .got_make_output_dir("STEP2_output")

ts_msg("STEP output directory: ", normalizePath(STEP_OUTPUT_DIR, winslash = "/", mustWork = FALSE))


# =============================================================================
# ==== 3) LOADERS + ANNOTATION =================================================
# =============================================================================

# ---- Trial meta builder ----
build_trial_meta <- function(trial_summary,
                             meta_trial_col = "Trial",
                             meta_time_col  = "timepoint",
                             meta_cols      = c("Group", "Treatment", "Date")) {
  required_meta <- c(meta_trial_col, meta_time_col, meta_cols)
  missing_meta <- setdiff(required_meta, names(trial_summary))
  if (length(missing_meta) > 0) {
    stop("trial_summary missing required cols: ", paste(missing_meta, collapse = ", "), call. = FALSE)
  }
  
  trial_summary %>%
    select(all_of(c(meta_trial_col, meta_time_col, meta_cols))) %>%
    mutate(
      .trial_key = suppressWarnings(as.integer(.data[[meta_trial_col]])),
      .time_key  = suppressWarnings(as.integer(.data[[meta_time_col]]))
    ) %>%
    filter(!is.na(.trial_key), !is.na(.time_key)) %>%
    distinct(.trial_key, .time_key, .keep_all = TRUE) %>%
    select(-all_of(c(meta_trial_col, meta_time_col)))
}

# ---- Annotate any by-frame-ish table that has trial_id + (segment_id|timepoint) ----
annotate_by_frame_df <- function(df,
                                 meta_tbl,
                                 df_trial_col = "trial_id",
                                 df_time_col_candidates = c("segment_id", "timepoint"),
                                 rename_time_to = "timepoint",
                                 rename_date_to = "trial_date",
                                 on_missing_keys = c("warn_keep", "keep_silent", "null")) {
  on_missing_keys <- match.arg(on_missing_keys)
  if (is.null(df)) return(NULL)
  
  df_time_col <- df_time_col_candidates[df_time_col_candidates %in% names(df)][1]
  if (!(df_trial_col %in% names(df)) || is.na(df_time_col)) {
    msg <- paste0("Skipping annotation: df missing join keys. Needed '", df_trial_col,
                  "' and one of {", paste(df_time_col_candidates, collapse = ", "), "}.")
    if (on_missing_keys == "null") return(NULL)
    if (on_missing_keys == "warn_keep") warning(msg, call. = FALSE)
    return(df)
  }
  
  out <- df %>%
    mutate(
      .trial_key = suppressWarnings(as.integer(.data[[df_trial_col]])),
      .time_key  = suppressWarnings(as.integer(.data[[df_time_col]]))
    )
  
  if (df_time_col != rename_time_to && !(rename_time_to %in% names(out))) {
    out <- rename(out, !!rename_time_to := !!df_time_col)
  }
  
  out <- out %>%
    left_join(meta_tbl, by = c(".trial_key" = ".trial_key", ".time_key" = ".time_key")) %>%
    select(-.trial_key, -.time_key)
  
  if ("Date" %in% names(out) && !(rename_date_to %in% names(out))) {
    out <- rename(out, !!rename_date_to := Date)
  }
  
  out
}



# =============================================================================
# ==== 4) LOAD + ANNOTATE (ROBUST, STRICTER) ==================================
# =============================================================================

ts_msg("Checking .GlobalEnv first for canonical objects, then clipboard path if needed...")

# ---- Trial summary / meta ----
trial_summary <- get_global_or_load(
  "trial_summary",
  {
    ts_msg("Reading trial summary Excel: ", TRIAL_SUMMARY_PATH)
    readxl::read_excel(TRIAL_SUMMARY_PATH)
  },
  validate_fn = function(x) is.data.frame(x) && nrow(x) > 0,
  force_reload = FORCE_RELOAD_FROM_DISK
)

meta_tbl <- get_global_or_load(
  "meta_tbl",
  {
    build_trial_meta(trial_summary)
  },
  validate_fn = function(x) is.data.frame(x) && all(c(".trial_key", ".time_key") %in% names(x)),
  force_reload = FORCE_RELOAD_FROM_DISK
)

# ---- Resolve clipboard path only if canonical master object is absent or reload forced ----
have_master_global <- FALSE

if (!isTRUE(FORCE_RELOAD_FROM_DISK) &&
    exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
  master_candidate <- get("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)
  have_master_global <- is_valid_master_fish_by_frame(master_candidate)
  
  if (!have_master_global) {
    warning(
      "master_fish_by_frame exists in .GlobalEnv but is not valid for reuse; will fall back to disk load.",
      call. = FALSE
    )
  }
}

if (isTRUE(FORCE_RELOAD_FROM_DISK) || !have_master_global) {
  ts_msg("Need disk access for canonical master_fish_by_frame -> reading master_outputs_dir (or by_frame_dir) from clipboard...")
  
  path_info <- resolve_by_frame_from_clipboard()
  root_dir     <- path_info$root_dir
  by_frame_dir <- path_info$by_frame_dir
  morph_dir    <- path_info$morph_dir
  
} else {
  ts_msg("Using valid master_fish_by_frame from .GlobalEnv")
  root_dir     <- NULL
  by_frame_dir <- NULL
  morph_dir    <- NULL
}

# ---- Load / reuse canonical master_fish_by_frame only ----
master_fish_by_frame <- get_global_or_load(
  "master_fish_by_frame",
  {
    local_by_frame_dir <- by_frame_dir
    
    if (is.null(local_by_frame_dir)) {
      ts_msg("master_fish_by_frame requires disk load -> resolving by_frame_dir from clipboard on demand...")
      path_info <- resolve_by_frame_from_clipboard()
      root_dir     <- path_info$root_dir
      by_frame_dir <- path_info$by_frame_dir
      morph_dir    <- path_info$morph_dir
      local_by_frame_dir <- by_frame_dir
    }
    
    ts_msg("Locating master_fish_by_frame (csv/xlsx)...")
    fish_master_path <- find_file_recursive(
      root_dir = local_by_frame_dir,
      basename_regex = "^master_fish_by_frame\\.(csv|xlsx|xls)$"
    )
    if (is.null(fish_master_path)) {
      stop("Could not find master_fish_by_frame.(csv|xlsx|xls) under: ", local_by_frame_dir, call. = FALSE)
    }
    ts_msg("Found: ", fish_master_path)
    
    out <- read_table_safe(fish_master_path)
    if (is.null(out) || nrow(out) == 0) {
      stop("master_fish_by_frame loaded as NULL/empty: ", fish_master_path, call. = FALSE)
    }
    out
  },
  validate_fn = is_valid_master_fish_by_frame,
  force_reload = FORCE_RELOAD_FROM_DISK
)

fish_cols_needed <- c(
  "session_id","trial_id","fish_id","fish_in_session","frame","time",
  "x_interp","y_interp","notes","speed_cm_s","heading_rad"
)
require_cols(master_fish_by_frame, fish_cols_needed, "master_fish_by_frame")

if (!any(c("segment_id", "timepoint") %in% names(master_fish_by_frame))) {
  stop("master_fish_by_frame must contain either 'segment_id' or 'timepoint'.", call. = FALSE)
}

# ---- Normalize fish keys including frame -----------------------------------
norm_fish_keys_with_frame <- function(df) {
  df %>%
    mutate(
      session_id = trimws(as.character(session_id)),
      trial_id   = suppressWarnings(as.integer(trial_id)),
      timepoint  = suppressWarnings(as.integer(timepoint)),
      fish_in_session = suppressWarnings(as.integer(fish_in_session)),
      fish_id    = trimws(as.character(fish_id)),
      frame      = suppressWarnings(as.integer(frame))
    )
}


# ---- Normalize fish keys excluding frame -----------------------------------
norm_fish_keys_no_frame <- function(df) {
  df %>%
    mutate(
      session_id = trimws(as.character(session_id)),
      trial_id   = suppressWarnings(as.integer(trial_id)),
      timepoint  = suppressWarnings(as.integer(timepoint)),
      fish_in_session = suppressWarnings(as.integer(fish_in_session)),
      fish_id    = trimws(as.character(fish_id))
    )
}

# ---- Annotate canonical master_fish_by_frame if metadata columns are absent ----
if (!all(c("Group", "Treatment", "trial_date", "timepoint") %in% names(master_fish_by_frame))) {
  ts_msg("Annotating master_fish_by_frame with metadata...")
  master_fish_by_frame <- annotate_by_frame_df(
    master_fish_by_frame,
    meta_tbl = meta_tbl,
    df_trial_col = "trial_id",
    df_time_col_candidates = c("segment_id", "timepoint"),
    rename_time_to = "timepoint",
    on_missing_keys = "warn_keep"
  )
} else {
  ts_msg("master_fish_by_frame already contains metadata columns; keeping existing version.")
}

assign("trial_summary", trial_summary, envir = .GlobalEnv)
assign("meta_tbl", meta_tbl, envir = .GlobalEnv)
assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)

# ---- OPT-0 downstream row filter: apply KEEP_TP1 / DROP_PHARMA to master data ----
{
  .n0 <- nrow(master_fish_by_frame)
  if (exists("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE) &&
      !isTRUE(get("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE))) {
    master_fish_by_frame <- master_fish_by_frame %>%
      dplyr::filter(as.integer(timepoint) != 1L)
    ts_msg("OPT-0 downstream: removed ", .n0 - nrow(master_fish_by_frame),
           " TP1 rows (KEEP_TP1=FALSE); ", nrow(master_fish_by_frame), " rows remain")
    .n0 <- nrow(master_fish_by_frame)
  }
  if (exists("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE) &&
      isTRUE(get("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE)) &&
      exists("PHARMA_GROUPS", envir = .GlobalEnv, inherits = FALSE)) {
    .pg <- get("PHARMA_GROUPS", envir = .GlobalEnv, inherits = FALSE)
    master_fish_by_frame <- master_fish_by_frame %>%
      dplyr::filter(!Group %in% .pg)
    ts_msg("OPT-0 downstream: removed ", .n0 - nrow(master_fish_by_frame),
           " pharma rows (Group %in% ", paste(.pg, collapse = ","),
           "); ", nrow(master_fish_by_frame), " rows remain")
  }
  assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)
  rm(.n0)
}

# ---- Always standardize + derive downstream canonical tables from master_fish_by_frame ----
ts_msg("Standardizing master_fish_by_frame and re-deriving position_by_frame / speed_by_frame / zone_by_frame...")

fish_bf <- master_fish_by_frame %>%
  mutate(
    session_id = trimws(as.character(session_id)),
    trial_id   = suppressWarnings(as.integer(trial_id)),
    timepoint  = suppressWarnings(as.integer(timepoint)),
    fish_id    = trimws(as.character(fish_id)),
    fish_in_session = suppressWarnings(as.integer(fish_in_session)),
    frame      = suppressWarnings(as.integer(frame)),
    time       = suppressWarnings(as.numeric(time)),
    x_interp   = suppressWarnings(as.numeric(x_interp)),
    y_interp   = suppressWarnings(as.numeric(y_interp)),
    speed_cm_s = suppressWarnings(as.numeric(speed_cm_s)),
    heading_rad = suppressWarnings(as.numeric(heading_rad))
  )

position_by_frame <- fish_bf %>%
  select(
    session_id, trial_id, timepoint, fish_id, fish_in_session, frame, time,
    x_interp, y_interp, notes,
    Group, Treatment, trial_date
  )

speed_by_frame <- fish_bf %>%
  transmute(
    individual = fish_id,
    frame, time,
    x = x_interp, y = y_interp,
    session_id, trial_id, timepoint, fish_in_session, fish_id,
    x_interp, y_interp, notes,
    speed_cm_s,
    Group, Treatment, trial_date
  )

zone_by_frame <- fish_bf %>%
  transmute(
    session_id, trial_id, timepoint, fish_id, fish_in_session,
    frame, time,
    zone,
    Group, Treatment, trial_date
  )

assign("position_by_frame", position_by_frame, envir = .GlobalEnv)
assign("speed_by_frame",    speed_by_frame,    envir = .GlobalEnv)
assign("zone_by_frame",     zone_by_frame,     envir = .GlobalEnv)

# ---- Reuse or load group_by_frame separately ----
group_by_frame <- get_global_or_load(
  "group_by_frame",
  {
    local_by_frame_dir <- by_frame_dir
    local_root_dir     <- root_dir
    
    if (is.null(local_by_frame_dir)) {
      ts_msg("group_by_frame missing in .GlobalEnv -> resolving by_frame_dir from clipboard on demand...")
      path_info <- resolve_by_frame_from_clipboard()
      local_root_dir     <- path_info$root_dir
      local_by_frame_dir <- path_info$by_frame_dir
      
      if (is.null(root_dir))     root_dir     <- path_info$root_dir
      if (is.null(by_frame_dir)) by_frame_dir <- path_info$by_frame_dir
      if (is.null(morph_dir))    morph_dir    <- path_info$morph_dir
    }
    
    ts_msg("Locating master_group_dynamics_by_frame.csv under group_* folders...")
    group_dyn_path <- find_file_recursive(
      root_dir = local_by_frame_dir,
      basename_regex = "^master_group_dynamics_by_frame\\.csv$",
      must_contain_regex = "/group_"
    )
    
    if (is.null(group_dyn_path)) {
      warning(
        "Could not find master_group_dynamics_by_frame.csv under any group_* folder.\n",
        "Group-level outputs will be unavailable (group_by_frame not created).",
        call. = FALSE
      )
      NULL
    } else {
      ts_msg("Found: ", group_dyn_path)
      out <- read_csv_safe(group_dyn_path)
      out <- annotate_by_frame_df(
        out,
        meta_tbl = meta_tbl,
        df_trial_col = "trial_id",
        df_time_col_candidates = c("segment_id", "timepoint"),
        rename_time_to = "timepoint",
        on_missing_keys = "warn_keep"
      )
      out
    }
  },
  validate_fn = function(x) is.null(x) || is.data.frame(x),
  force_reload = FORCE_RELOAD_FROM_DISK
)

if (!is.null(group_by_frame)) {
  assign("group_by_frame", group_by_frame, envir = .GlobalEnv)
}

# ---- Morph objects: reuse if present; otherwise load if path available ----
have_bodylength_by_fish <- !isTRUE(FORCE_RELOAD_FROM_DISK) &&
  exists("bodylength_by_fish", envir = .GlobalEnv, inherits = FALSE)

have_bodylength_summary <- !isTRUE(FORCE_RELOAD_FROM_DISK) &&
  exists("bodylength_summary_by_session", envir = .GlobalEnv, inherits = FALSE)

if (isTRUE(FORCE_RELOAD_FROM_DISK) || !have_bodylength_by_fish || !have_bodylength_summary) {
  
  if (is.null(morph_dir)) {
    ts_msg("Morph objects missing in .GlobalEnv -> resolving morph_dir from clipboard on demand...")
    path_info <- resolve_by_frame_from_clipboard()
    morph_dir <- path_info$morph_dir
    
    if (is.null(by_frame_dir)) by_frame_dir <- path_info$by_frame_dir
    if (is.null(root_dir))     root_dir     <- path_info$root_dir
  }
  
  if (!is.null(morph_dir) && dir.exists(morph_dir)) {
    ts_msg("Loading morph outputs from: ", morph_dir)
    
    morph_bodylen_path <- find_file_recursive(morph_dir, "^master_bodylength_by_fish\\.csv$")
    morph_summary_path <- find_file_recursive(morph_dir, "^master_bodylength_summary_by_session\\.csv$")
    
    bodylength_by_fish <- get_global_or_load(
      "bodylength_by_fish",
      {
        out <- read_csv_safe(morph_bodylen_path)
        if (!is.null(out)) {
          if (!("timepoint" %in% names(out)) && ("segment_id" %in% names(out))) {
            out <- dplyr::rename(out, timepoint = segment_id)
          }
          out <- annotate_by_frame_df(
            out,
            meta_tbl = meta_tbl,
            df_trial_col = "trial_id",
            df_time_col_candidates = c("timepoint", "segment_id"),
            rename_time_to = "timepoint",
            on_missing_keys = "warn_keep"
          )
        }
        out
      },
      validate_fn = function(x) is.null(x) || is.data.frame(x),
      force_reload = FORCE_RELOAD_FROM_DISK
    )
    
    bodylength_summary_by_session <- get_global_or_load(
      "bodylength_summary_by_session",
      {
        out <- read_csv_safe(morph_summary_path)
        if (!is.null(out)) {
          if (!("timepoint" %in% names(out)) && ("segment_id" %in% names(out))) {
            out <- dplyr::rename(out, timepoint = segment_id)
          }
          out <- annotate_by_frame_df(
            out,
            meta_tbl = meta_tbl,
            df_trial_col = "trial_id",
            df_time_col_candidates = c("timepoint", "segment_id"),
            rename_time_to = "timepoint",
            on_missing_keys = "warn_keep"
          )
        }
        out
      },
      validate_fn = function(x) is.null(x) || is.data.frame(x),
      force_reload = FORCE_RELOAD_FROM_DISK
    )
    
    if (!is.null(bodylength_by_fish)) {
      assign("bodylength_by_fish", bodylength_by_fish, envir = .GlobalEnv)
    }
    if (!is.null(bodylength_summary_by_session)) {
      assign("bodylength_summary_by_session", bodylength_summary_by_session, envir = .GlobalEnv)
    }
    
  } else {
    warning("morph directory not found (skipping): ", morph_dir, call. = FALSE)
  }
} else {
  ts_msg("Using existing morph objects from .GlobalEnv: bodylength_by_fish, bodylength_summary_by_session")
}

ts_msg(
  "Canonical objects ready: master_fish_by_frame, position_by_frame, speed_by_frame, zone_by_frame",
  if (!is.null(group_by_frame)) ", group_by_frame" else ""
)
# =============================================================================
# ==== 5) ACTIVITY TABLES + APPEND TO master_fish_by_frame (CLEAN BLOCK) ======
# =============================================================================

ts_msg("Building activity tables and updating master_fish_by_frame (clean rerunnable block)...")

# ---- Required objects --------------------------------------------------------
stopifnot(exists("speed_by_frame", envir = .GlobalEnv, inherits = FALSE))
stopifnot(exists("bodylength_by_fish", envir = .GlobalEnv, inherits = FALSE))
stopifnot(exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE))

speed_by_frame     <- get("speed_by_frame", envir = .GlobalEnv)
bodylength_by_fish <- get("bodylength_by_fish", envir = .GlobalEnv)
mf                 <- get("master_fish_by_frame", envir = .GlobalEnv)

# ── Resolve idmatcher_id_map: GlobalEnv first, then newest STEP1_5 run on disk ─
# STEP 1.5 assigns idmatcher_id_map into .GlobalEnv, but only runs when the
# pipeline starts at STEP 1. Resuming at STEP 2 (run_pipeline_from_step2.R) left
# it absent, and the two consumers below then silently fell back to per-segment
# body length — producing NA fish_id_matched_trial and the wrong
# speed_threshold_TS* for every matched fish. That is how the 2026-05-24 vintage
# came to have segment-BL thresholds for all 96 fish (see LOGS.md 2026-08-21).
# The disk fallback makes a STEP 2 resume equivalent to a full run; a genuinely
# absent map now warns instead of passing by in a message.
if (!exists("idmatcher_id_map", envir = .GlobalEnv, inherits = FALSE)) {
  .idmap_dir <- if (exists(".got_find_latest_output", mode = "function")) {
    .got_find_latest_output("STEP1_5")
  } else NULL
  .idmap_csv <- if (!is.null(.idmap_dir)) {
    sort(list.files(.idmap_dir, pattern = "^idmatcher_id_map_.*\\.csv$",
                    full.names = TRUE))
  } else character(0)
  if (length(.idmap_csv) > 0L) {
    .idmap_csv <- .idmap_csv[length(.idmap_csv)]
    assign("idmatcher_id_map",
           suppressMessages(readr::read_csv(.idmap_csv, show_col_types = FALSE)),
           envir = .GlobalEnv)
    ts_msg("STEP2: idmatcher_id_map not in GlobalEnv — loaded from disk: ", .idmap_csv)
  }
}

# ── idmatcher: join fish_id_matched_trial + body_length_mean_cm onto mf ──────
# body_length_mean_cm is the mean of seg2 + seg3 BL for matched fish.
# Speed thresholds (speed_threshold_TS*) are recomputed here for matched fish
# so make_activity_alias() picks up the canonical values automatically.
if (exists("idmatcher_id_map", envir = .GlobalEnv, inherits = FALSE)) {
  .idmap_thr <- idmatcher_id_map |>
    dplyr::select(tracker_fish_id, fish_id_matched_trial, body_length_mean_cm) |>
    dplyr::filter(!is.na(fish_id_matched_trial))
  mf <- mf |>
    dplyr::left_join(.idmap_thr, by = c("fish_id" = "tracker_fish_id"))
  for (.ts_mult in c(0.5, 1.0, 1.5, 2.0)) {
    .thr_col <- paste0("speed_threshold_TS", .ts_mult)
    .ls_col  <- paste0("locomotor_status_TS", .ts_mult)
    if (.thr_col %in% names(mf)) {
      mf <- mf |>
        dplyr::mutate(
          !!.thr_col := dplyr::if_else(
            !is.na(body_length_mean_cm),
            body_length_mean_cm * .ts_mult,
            .data[[.thr_col]]
          ),
          !!.ls_col := dplyr::if_else(
            is.finite(speed_cm_s) & is.finite(.data[[.thr_col]]),
            speed_cm_s >= .data[[.thr_col]],
            NA
          )
        )
    }
  }
  ts_msg("idmatcher: speed_threshold_TS* recomputed for matched fish using mean BL.")
} else {
  mf <- mf |>
    dplyr::mutate(fish_id_matched_trial = NA_character_,
                  body_length_mean_cm   = NA_real_)
  warning("STEP2: idmatcher_id_map not found in GlobalEnv or under STEP1_5/ — ",
          "fish_id_matched_trial will be NA for ALL fish and speed_threshold_TS* ",
          "will fall back to per-segment body length. Run STEP 1.5 first.",
          call. = FALSE)
}

# ---- Join keys ---------------------------------------------------------------
mf_join    <- c("session_id","trial_id","timepoint","fish_in_session","fish_id","frame")
mf_bl_join <- c("session_id","trial_id","timepoint","fish_in_session","fish_id")

# ---- Required columns --------------------------------------------------------
speed_cols_needed <- c(
  "session_id", "trial_id", "timepoint", "fish_in_session", "fish_id",
  "frame", "time", "x_interp", "y_interp", "notes",
  "speed_cm_s", "Group", "Treatment", "trial_date"
)

bodylen_cols_needed <- c(
  "session_id", "trial_id", "timepoint", "fish_in_session",
  "fish_id", "body_length_cm"
)

mf_cols_needed <- c(
  "session_id","trial_id","timepoint","fish_in_session","fish_id","frame","speed_cm_s"
)

require_cols(speed_by_frame, speed_cols_needed, "speed_by_frame")
require_cols(bodylength_by_fish, bodylen_cols_needed, "bodylength_by_fish")
require_cols(mf, mf_cols_needed, "master_fish_by_frame")

# ---- Standardize framewise speed table --------------------------------------
activity_base <- speed_by_frame %>%
  select(any_of(c(
    "session_id", "trial_id", "timepoint", "fish_in_session", "fish_id",
    "frame", "time", "x_interp", "y_interp", "notes",
    "speed_cm_s", "Group", "Treatment", "trial_date"
  ))) %>%
  norm_fish_keys_with_frame() %>%
  mutate(
    time        = suppressWarnings(as.numeric(time)),
    x_interp    = suppressWarnings(as.numeric(x_interp)),
    y_interp    = suppressWarnings(as.numeric(y_interp)),
    speed_cm_s  = suppressWarnings(as.numeric(speed_cm_s)),
    Group       = as.character(Group),
    Treatment   = as.character(Treatment)
  ) %>%
  filter(!is.na(speed_cm_s))

# ---- Standardize body-length lookup -----------------------------------------
bodylen_lookup <- bodylength_by_fish %>%
  select(any_of(c(
    "session_id", "trial_id", "timepoint", "fish_in_session", "fish_id",
    "body_length_cm", "body_length_sd_cm"
  ))) %>%
  norm_fish_keys_no_frame() %>%
  mutate(
    body_length_cm = suppressWarnings(as.numeric(body_length_cm)),
    body_length_sd_cm = if ("body_length_sd_cm" %in% names(.)) {
      suppressWarnings(as.numeric(body_length_sd_cm))
    } else {
      NA_real_
    }
  ) %>%
  filter(is.finite(body_length_cm), body_length_cm > 0) %>%
  group_by(across(all_of(mf_bl_join))) %>%
  summarise(
    body_length_cm    = first(body_length_cm),
    body_length_sd_cm = first(body_length_sd_cm),
    .groups = "drop"
  )

# ---- Join body length onto framewise speed, add idmatcher canonical BL ------
activity_base <- activity_base %>%
  left_join(bodylen_lookup, by = mf_bl_join)

# ── idmatcher: add fish_id_matched_trial + body_length_canonical_cm ───────────
if (exists("idmatcher_id_map", envir = .GlobalEnv, inherits = FALSE)) {
  .idmap_act <- idmatcher_id_map |>
    dplyr::select(tracker_fish_id, fish_id_matched_trial, body_length_mean_cm)
  activity_base <- activity_base |>
    dplyr::left_join(.idmap_act, by = c("fish_id" = "tracker_fish_id")) |>
    dplyr::mutate(
      body_length_canonical_cm = dplyr::coalesce(body_length_mean_cm, body_length_cm),
      body_length_source = dplyr::if_else(!is.na(body_length_mean_cm),
                                          "seg2_seg3_mean", "per_segment_fallback")
    )
} else {
  activity_base <- activity_base |>
    dplyr::mutate(
      fish_id_matched_trial    = NA_character_,
      body_length_mean_cm      = NA_real_,
      body_length_canonical_cm = body_length_cm,
      body_length_source       = "per_segment_fallback"
    )
}

activity_base <- activity_base %>%
  mutate(
    speed_BL_s = dplyr::if_else(
      is.finite(speed_cm_s) & is.finite(body_length_canonical_cm) & body_length_canonical_cm > 0,
      speed_cm_s / body_length_canonical_cm,
      NA_real_
    )
  )

# ---- locomotor_status_TS* comes from STEP1 (idtracker_pipelineSTEP1GOT) -----
# build_activity_table() was removed: STEP1 now computes speed_threshold_TS* and
# locomotor_status_TS* per frame and writes them into master_fish_by_frame.
# STEP2 reads those columns directly from mf (loaded above); the activity_TS*
# table aliases below preserve backward-compatibility for STEP2's join loop.

ts_msg("Aliasing activity_TS* tables from master_fish_by_frame (locomotor_status already in STEP1 output)...")

make_activity_alias <- function(mf, thr_label) {
  ls_col <- paste0("locomotor_status_", thr_label)
  st_col <- paste0("speed_threshold_",  thr_label)
  if (!ls_col %in% names(mf)) {
    warning("Column '", ls_col, "' not found in master_fish_by_frame. ",
            "Ensure STEP1 ran with locomotor_status computation enabled.",
            call. = FALSE)
    return(NULL)
  }
  result <- mf %>%
    dplyr::select(
      dplyr::all_of(c("session_id", "trial_id", "fish_id", "fish_in_session", "frame")),
      dplyr::any_of(c(st_col, ls_col, "speed_cm_s", "body_length_cm",
                      "Group", "Treatment", "trial_date", "segment_id", "timepoint", "time"))
    )

  # Rename the activity columns if they exist
  if (st_col %in% names(result)) {
    names(result)[names(result) == st_col] <- "speed_threshold"
  }
  if (ls_col %in% names(result)) {
    names(result)[names(result) == ls_col] <- "locomotor_status"
  }

  result
}

activity_TS0.5 <- make_activity_alias(mf, "TS0.5")
activity_TS1   <- make_activity_alias(mf, "TS1")
activity_TS1.5 <- make_activity_alias(mf, "TS1.5")
activity_TS2   <- make_activity_alias(mf, "TS2")

assign("activity_TS0.5", activity_TS0.5, envir = .GlobalEnv)
assign("activity_TS1",   activity_TS1,   envir = .GlobalEnv)
assign("activity_TS1.5", activity_TS1.5, envir = .GlobalEnv)
assign("activity_TS2",   activity_TS2,   envir = .GlobalEnv)

ts_msg("Activity tables ready: activity_TS0.5, activity_TS1, activity_TS1.5, activity_TS2")

# ---- Clean master_fish_by_frame before appending derived columns ------------
mf <- mf %>%
  norm_fish_keys_with_frame() %>%
  mutate(
    speed_cm_s = suppressWarnings(as.numeric(speed_cm_s))
  ) %>%
  select(
    -any_of(c(
      "body_length_cm", "body_length_sd_cm", "speed_BL_s",
      "speed_threshold_TS0.5", "locomotor_status_TS0.5",
      "speed_threshold_TS1",   "locomotor_status_TS1",
      "speed_threshold_TS1.5", "locomotor_status_TS1.5",
      "speed_threshold_TS2",   "locomotor_status_TS2"
    )),
    -matches("^body_length_cm\\.(x|y)$"),
    -matches("^body_length_sd_cm\\.(x|y)$"),
    -matches("^speed_threshold_TS(0\\.5|1|1\\.5|2)\\.(x|y)$"),
    -matches("^locomotor_status_TS(0\\.5|1|1\\.5|2)\\.(x|y)$")
  )

# ---- Append activity columns one table at a time ----------------------------
ts_datasets <- c("activity_TS0.5", "activity_TS1", "activity_TS1.5", "activity_TS2")

for (nm in ts_datasets) {
  df_act <- get(nm, envir = .GlobalEnv)
  
  require_cols(df_act, c(mf_join, "speed_threshold", "locomotor_status"), nm)
  
  suffix <- sub("^activity_", "", nm)
  st_col <- paste0("speed_threshold_", suffix)
  ls_col <- paste0("locomotor_status_", suffix)
  
  df_act2 <- df_act %>%
    norm_fish_keys_with_frame() %>%
    transmute(
      session_id,
      trial_id,
      timepoint,
      fish_in_session,
      fish_id,
      frame,
      !!st_col := suppressWarnings(as.numeric(speed_threshold)),
      !!ls_col := as.character(locomotor_status)
    ) %>%
    group_by(across(all_of(mf_join))) %>%
    summarise(
      !!st_col := first(.data[[st_col]]),
      !!ls_col := first(.data[[ls_col]]),
      .groups = "drop"
    )
  
  # ---- Diagnostics -----------------------------------------------------------
  mf_keys <- mf %>% distinct(across(all_of(mf_join)))
  act_keys <- df_act2 %>% distinct(across(all_of(mf_join)))
  
  missing_in_activity <- mf_keys %>%
    anti_join(act_keys, by = mf_join) %>%
    arrange(trial_id, timepoint, fish_in_session, fish_id, frame)
  
  activity_not_in_mf <- act_keys %>%
    anti_join(mf_keys, by = mf_join) %>%
    arrange(trial_id, timepoint, fish_in_session, fish_id, frame)
  
  cat("\n================ ACTIVITY JOIN DIAGNOSTIC:", suffix, "=================\n")
  
  print_df_safe(
    missing_in_activity %>%
      count(session_id, trial_id, timepoint, fish_in_session, fish_id, name = "n_frames_missing") %>%
      arrange(trial_id, timepoint, fish_in_session, fish_id),
    paste0("-- Keys present in master_fish_by_frame but missing in ", nm, " --")
  )
  
  print_df_safe(
    activity_not_in_mf %>%
      count(session_id, trial_id, timepoint, fish_in_session, fish_id, name = "n_frames_extra") %>%
      arrange(trial_id, timepoint, fish_in_session, fish_id),
    paste0("-- Keys present in ", nm, " but not in master_fish_by_frame --")
  )
  
  print_df_safe(
    missing_in_activity %>%
      count(trial_id, timepoint, name = "n_frames_missing") %>%
      arrange(trial_id, timepoint),
    "-- Summary of missing master_fish_by_frame -> activity join keys by trial/timepoint --"
  )
  
  cat("================================================================\n\n")
  
  readr::write_csv(
    missing_in_activity,
    file.path(STEP_OUTPUT_DIR, paste0("diagnostic_activity_join_missing_in_", suffix, ".csv"))
  )
  
  readr::write_csv(
    activity_not_in_mf,
    file.path(STEP_OUTPUT_DIR, paste0("diagnostic_activity_join_extra_in_", suffix, ".csv"))
  )
  
  readr::write_csv(
    missing_in_activity %>%
      count(session_id, trial_id, timepoint, fish_in_session, fish_id, name = "n_frames_missing") %>%
      arrange(trial_id, timepoint, fish_in_session, fish_id),
    file.path(STEP_OUTPUT_DIR, paste0("diagnostic_activity_join_missing_summary_by_fish_", suffix, ".csv"))
  )
  
  readr::write_csv(
    missing_in_activity %>%
      count(trial_id, timepoint, name = "n_frames_missing") %>%
      arrange(trial_id, timepoint),
    file.path(STEP_OUTPUT_DIR, paste0("diagnostic_activity_join_missing_summary_by_trial_tp_", suffix, ".csv"))
  )
  
  before_n <- nrow(mf)
  mf <- mf %>% left_join(df_act2, by = mf_join)
  stopifnot(nrow(mf) == before_n)
  
  ts_msg(" - Added ", st_col, " + ", ls_col, " | missing rows after join: ", sum(is.na(mf[[ls_col]])))
}

# ---- Append body_length_cm + body_length_sd_cm + canonical speed_BL_s -------
ts_msg("Appending body_length_cm + body_length_sd_cm + speed_BL_s (canonical) to master_fish_by_frame...")

before_n <- nrow(mf)

mf <- mf %>%
  left_join(bodylen_lookup, by = mf_bl_join) %>%
  mutate(
    body_length_cm    = suppressWarnings(as.numeric(body_length_cm)),
    body_length_sd_cm = suppressWarnings(as.numeric(body_length_sd_cm)),
    speed_cm_s        = suppressWarnings(as.numeric(speed_cm_s)),
    # body_length_mean_cm was joined earlier (idmatcher block at top of STEP2);
    # canonical = mean BL for matched fish, per-segment fallback for unmatched
    body_length_canonical_cm = dplyr::coalesce(body_length_mean_cm, body_length_cm),
    body_length_source = dplyr::if_else(!is.na(body_length_mean_cm),
                                        "seg2_seg3_mean", "per_segment_fallback"),
    speed_BL_s = dplyr::if_else(
      is.finite(speed_cm_s) & is.finite(body_length_canonical_cm) & body_length_canonical_cm > 0,
      speed_cm_s / body_length_canonical_cm,
      NA_real_
    )
  )

stopifnot(nrow(mf) == before_n)

# ---- Final guard against join artifacts -------------------------------------
bad_dup_cols <- names(mf)[grepl("\\.(x|y)$", names(mf))]
if (length(bad_dup_cols) > 0) {
  stop(
    "Duplicate join-artifact columns remain in master_fish_by_frame: ",
    paste(bad_dup_cols, collapse = ", "),
    call. = FALSE
  )
}

ts_msg(" - body_length_cm missing rows: ", sum(is.na(mf$body_length_cm)))
ts_msg(" - speed_BL_s missing rows: ", sum(is.na(mf$speed_BL_s)))

assign("master_fish_by_frame", mf, envir = .GlobalEnv)
ts_msg("Done. master_fish_by_frame updated cleanly with activity/body-length columns.")

mf_activity_path <- file.path(
  STEP_OUTPUT_DIR,
  paste0("master_fish_by_frame_with_activity_", timestamp_tag(), ".csv")
)
readr::write_csv(mf, mf_activity_path)
ts_msg("Saved: ", mf_activity_path)
# =============================================================================
# ==== 5d) REFRESH speed_by_frame WITH speed_BL_s =============================
# =============================================================================

ts_msg("Refreshing speed_by_frame from updated master_fish_by_frame...")

require_cols(
  mf,
  c(
    "session_id","trial_id","timepoint","fish_in_session","fish_id","frame","time",
    "x_interp","y_interp","notes","speed_cm_s","speed_BL_s","Group","Treatment","trial_date"
  ),
  "master_fish_by_frame"
)

speed_by_frame <- mf %>%
  transmute(
    individual = fish_id,
    frame, time,
    x = x_interp, y = y_interp,
    session_id, trial_id, timepoint, fish_in_session, fish_id,
    x_interp, y_interp, notes,
    speed_cm_s,
    speed_BL_s,
    Group, Treatment, trial_date
  )

assign("speed_by_frame", speed_by_frame, envir = .GlobalEnv)
ts_msg("Done. speed_by_frame refreshed with speed_BL_s.")

# =============================================================================
# ==== 6) DISTANCE PER FISH (meters + body lengths) ===========================
# =============================================================================
ts_msg("Computing total distance per fish from position_by_frame...")

position_by_frame  <- get("position_by_frame", envir = .GlobalEnv)
bodylength_by_fish <- get("bodylength_by_fish", envir = .GlobalEnv)

req_cols_pos <- c(
  "session_id","trial_id","timepoint","fish_id","fish_in_session",
  "Group","Treatment","trial_date","frame","x_interp","y_interp"
)
require_cols(position_by_frame, req_cols_pos, "position_by_frame")

req_cols_bl <- c(
  "session_id","trial_id","timepoint","fish_id","fish_in_session","body_length_cm"
)
require_cols(bodylength_by_fish, req_cols_bl, "bodylength_by_fish")

bl_dist <- bodylength_by_fish %>%
  mutate(
    session_id = trimws(as.character(session_id)),
    trial_id   = suppressWarnings(as.integer(trial_id)),
    timepoint  = suppressWarnings(as.integer(timepoint)),
    fish_id    = trimws(as.character(fish_id)),
    fish_in_session = suppressWarnings(as.integer(fish_in_session)),
    body_length_cm = suppressWarnings(as.numeric(body_length_cm))
  ) %>%
  filter(!is.na(body_length_cm), body_length_cm > 0) %>%
  group_by(session_id, trial_id, timepoint, fish_id, fish_in_session) %>%
  summarise(body_length_cm = first(body_length_cm), .groups = "drop")

distance_by_fish <- position_by_frame %>%
  mutate(
    session_id = trimws(as.character(session_id)),
    trial_id   = suppressWarnings(as.integer(trial_id)),
    timepoint  = suppressWarnings(as.integer(timepoint)),
    fish_id    = trimws(as.character(fish_id)),
    fish_in_session = suppressWarnings(as.integer(fish_in_session)),
    frame      = suppressWarnings(as.integer(frame)),
    x          = suppressWarnings(as.numeric(x_interp)),
    y          = suppressWarnings(as.numeric(y_interp))
  ) %>%
  group_by(session_id, trial_id, timepoint, fish_id, fish_in_session) %>%
  arrange(frame, .by_group = TRUE) %>%
  mutate(
    dx_px   = x - lag(x),
    dy_px   = y - lag(y),
    step_ok = is.finite(dx_px) & is.finite(dy_px),
    step_cm = if_else(step_ok, sqrt(dx_px^2 + dy_px^2) / LENGTH_UNIT_PX_PER_CM, NA_real_)
  ) %>%
  summarise(
    Group             = first(Group),
    Treatment         = first(Treatment),
    trial_date        = first(trial_date),
    total_distance_cm = sum(step_cm, na.rm = TRUE),
    distance_m        = total_distance_cm / 100,
    n_steps_used      = sum(!is.na(step_cm)),
    .groups = "drop"
  ) %>%
  left_join(
    bl_dist,
    by = c("session_id","trial_id","timepoint","fish_id","fish_in_session")
  ) %>%
  mutate(
    distance_bl = if_else(
      is.finite(total_distance_cm) & is.finite(body_length_cm) & body_length_cm > 0,
      total_distance_cm / body_length_cm,
      NA_real_
    )
  ) %>%
  select(
    session_id, trial_id, timepoint, fish_id, fish_in_session,
    Group, Treatment, trial_date,
    body_length_cm,
    total_distance_cm, distance_m, distance_bl,
    n_steps_used
  )

distance_csv_path <- file.path(STEP_OUTPUT_DIR, "distance_by_fish.csv")
write.csv(distance_by_fish, distance_csv_path, row.names = FALSE)
ts_msg("Saved: ", distance_csv_path)

assign("distance_by_fish", distance_by_fish, envir = .GlobalEnv)


# =============================================================================
# ==== 6b) DISTANCE PLOTS (meters + body lengths) -> PPT ======================
# =============================================================================
ts_msg("Exporting distance plots (meters + body lengths) to PPT...")

distance_by_fish <- get("distance_by_fish", envir = .GlobalEnv)

dist_cols_needed <- c(
  "fish_id","trial_id","timepoint","Group","Treatment","trial_date",
  "distance_m","distance_bl"
)
require_cols(distance_by_fish, dist_cols_needed, "distance_by_fish")

distance_ppt_file <- file.path(STEP_OUTPUT_DIR, paste0("distance_by_fish_", timestamp_tag(), ".pptx"))

make_distance_graph_A <- function(df_fish, y_col, y_lab, title_text, y_max, y_breaks) {
  df2 <- df_fish %>%
    mutate(
      Group     = factor(Group),
      timepoint = factor(timepoint, levels = sort(unique(timepoint))),
      Treatment = factor(Treatment),
      x_key = interaction(Group, timepoint, Treatment, sep = " | ", drop = TRUE),
      x_lab = paste0(Treatment, "\nTP ", timepoint, "\n", Group)
    )
  
  level_tbl <- df2 %>%
    distinct(Group, timepoint, Treatment, x_key, x_lab) %>%
    arrange(Group, timepoint, Treatment)
  
  df2 <- df2 %>% mutate(x_key = factor(x_key, levels = level_tbl$x_key))
  lab_map <- setNames(level_tbl$x_lab, level_tbl$x_key)
  
  ggplot(df2, aes(x = x_key, y = .data[[y_col]])) +
    geom_point(
      aes(color = Treatment),
      position = position_jitter(width = 0.12, height = 0),
      alpha = 0.7, size = 2
    ) +
    stat_summary(aes(fill = Treatment), fun = mean, geom = "col", width = 0.55, alpha = 0.35) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18, linewidth = 0.7) +
    scale_x_discrete(labels = lab_map) +
    scale_y_continuous(limits = c(0, y_max), breaks = y_breaks) +
    labs(x = NULL, y = y_lab, title = title_text) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(size = 10, lineheight = 0.95))
}

make_distance_graph_B <- function(df_fish, y_col, y_lab, title_text, y_max, y_breaks) {
  pd_bar <- position_dodge(width = 0.65)
  
  ggplot(df_fish, aes(x = Group, y = .data[[y_col]], fill = Treatment)) +
    geom_point(
      aes(color = Treatment),
      position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0, dodge.width = 0.65),
      alpha = 0.7, size = 2
    ) +
    stat_summary(fun = mean, geom = "col", width = 0.55, alpha = 0.35, position = pd_bar) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18, linewidth = 0.7, position = pd_bar) +
    facet_wrap(~ timepoint, nrow = 1) +
    scale_y_continuous(limits = c(0, y_max), breaks = y_breaks) +
    labs(x = "Group", y = y_lab, title = title_text) +
    theme_minimal(base_size = 13)
}

plot_distance_metric <- function(df, y_col, y_lab, label_prefix, ppt) {
  d <- df %>%
    mutate(
      Group = factor(Group),
      Treatment = factor(Treatment),
      timepoint = factor(as.character(timepoint), levels = sort(unique(as.character(timepoint))))
    ) %>%
    filter(is.finite(.data[[y_col]]))
  
  if (nrow(d) == 0) return(ppt)
  
  y_max_raw <- max(d[[y_col]], na.rm = TRUE)
  y_max <- if (is.finite(y_max_raw) && y_max_raw > 0) y_max_raw * 1.10 else 1
  y_breaks <- pretty(c(0, y_max), n = 10)
  
  pA <- make_distance_graph_A(d, y_col, y_lab, paste0(label_prefix, " — Graph A"), y_max, y_breaks)
  pB <- make_distance_graph_B(d, y_col, y_lab, paste0(label_prefix, " — Graph B"), y_max, y_breaks)
  
  ppt <- add_plot_slide(ppt, pA, paste0(label_prefix, " — Graph A"))
  ppt <- add_plot_slide(ppt, pB, paste0(label_prefix, " — Graph B"))
  ppt
}

add_plot_slide <- function(ppt, plot_obj, slide_title) {
  ppt %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = slide_title, location = ph_location_type(type = "title")) %>%
    ph_with(value = rvg::dml(ggobj = plot_obj), location = ph_location_type(type = "body"))
}

ppt_dist <- read_pptx()

ppt_dist <- plot_distance_metric(
  distance_by_fish,
  y_col = "distance_m",
  y_lab = "Total distance travelled (m)",
  label_prefix = "Distance travelled in meters",
  ppt = ppt_dist
)

ppt_dist <- plot_distance_metric(
  distance_by_fish,
  y_col = "distance_bl",
  y_lab = "Total distance travelled (body lengths)",
  label_prefix = "Distance travelled in body lengths",
  ppt = ppt_dist
)

print(ppt_dist, target = distance_ppt_file)
ts_msg("Saved PPT: ", distance_ppt_file)


# =============================================================================
# ==== 7) ACTIVITY BUDGET (MOVING SECONDS) -> PPT =============================
# =============================================================================

ts_msg("Exporting moving-seconds activity budget to PPT...")

dataset_names <- c("activity_TS0.5", "activity_TS1", "activity_TS1.5", "activity_TS2")
p_lo <- 0.01
p_hi <- 0.99
timepoint_to_drop <- "1"
ppt_file <- file.path(STEP_OUTPUT_DIR, paste0("activity_budget_GOT_", timestamp_tag(), ".pptx"))

assert_has_cols <- function(df, cols, df_name = "data") {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) stop(df_name, " missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

compute_moving_seconds_by_fish <- function(df, df_name = "data") {
  required <- c("fish_id","trial_id","timepoint","Group","Treatment","trial_date","time","locomotor_status")
  assert_has_cols(df, required, df_name)
  
  df %>%
    mutate(
      fish_id = as.character(fish_id),
      locomotor_status = as.character(locomotor_status),
      time = as.numeric(time),
      timepoint = as.character(timepoint),
      Group = as.character(Group),
      Treatment = as.character(Treatment),
      trial_date = as.character(trial_date)
    ) %>%
    group_by(fish_id) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(dt = time - lag(time), dt_ok = is.finite(dt) & dt > 0) %>%
    ungroup() %>%
    filter(dt_ok) %>%
    group_by(fish_id) %>%
    summarise(
      trial_id   = first(trial_id),
      timepoint  = first(timepoint),
      Group      = first(Group),
      Treatment  = first(Treatment),
      trial_date = first(trial_date),
      moving_seconds = sum(dt[locomotor_status == "moving"], na.rm = TRUE),
      n_dt_steps = n(),
      .groups = "drop"
    )
}

trim_outliers_percentile <- function(df, y_col, p_lo = 0.01, p_hi = 0.99) {
  y <- df[[y_col]]
  y_finite <- y[is.finite(y)]
  if (length(y_finite) == 0) {
    return(list(df_trim = df[0, ], removed_n = nrow(df), bounds = c(NA, NA)))
  }
  lo <- as.numeric(quantile(y_finite, probs = p_lo, na.rm = TRUE, names = FALSE))
  hi <- as.numeric(quantile(y_finite, probs = p_hi, na.rm = TRUE, names = FALSE))
  keep <- is.finite(y) & y >= lo & y <= hi
  list(df_trim = df[keep, , drop = FALSE], removed_n = sum(!keep), bounds = c(lo, hi))
}

make_graph_A <- function(df_fish, title_text, y_max, y_breaks) {
  df2 <- df_fish %>%
    mutate(
      Group     = factor(Group),
      timepoint = factor(timepoint, levels = sort(unique(timepoint))),
      Treatment = factor(Treatment),
      x_key = interaction(Group, timepoint, Treatment, sep = " | ", drop = TRUE),
      x_lab = paste0(Treatment, "\nTP ", timepoint, "\n", Group)
    )
  
  level_tbl <- df2 %>%
    distinct(Group, timepoint, Treatment, x_key, x_lab) %>%
    arrange(Group, timepoint, Treatment)
  
  df2 <- df2 %>% mutate(x_key = factor(x_key, levels = level_tbl$x_key))
  lab_map <- setNames(level_tbl$x_lab, level_tbl$x_key)
  
  ggplot(df2, aes(x = x_key, y = moving_seconds)) +
    geom_point(
      aes(color = Treatment),
      position = position_jitter(width = 0.12, height = 0),
      alpha = 0.7, size = 2
    ) +
    stat_summary(aes(fill = Treatment), fun = mean, geom = "col", width = 0.55, alpha = 0.35) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18, linewidth = 0.7) +
    scale_x_discrete(labels = lab_map) +
    scale_y_continuous(limits = c(0, y_max), breaks = y_breaks) +
    labs(x = NULL, y = "Moving time (s)", title = title_text) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(size = 10, lineheight = 0.95))
}

make_graph_B <- function(df_fish, title_text, y_max, y_breaks) {
  pd_bar <- position_dodge(width = 0.65)
  
  ggplot(df_fish, aes(x = Group, y = moving_seconds, fill = Treatment)) +
    geom_point(
      aes(color = Treatment),
      position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0, dodge.width = 0.65),
      alpha = 0.7, size = 2
    ) +
    stat_summary(fun = mean, geom = "col", width = 0.55, alpha = 0.35, position = pd_bar) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18, linewidth = 0.7, position = pd_bar) +
    facet_wrap(~ timepoint, nrow = 1) +
    scale_y_continuous(limits = c(0, y_max), breaks = y_breaks) +
    labs(x = "Group", y = "Moving time (s)", title = title_text) +
    theme_minimal(base_size = 13)
}

ppt <- read_pptx()
fish_summaries <- list()

for (nm in dataset_names) {
  if (!exists(nm, inherits = TRUE)) {
    warning("Dataset not found: ", nm, call. = FALSE)
    next
  }
  
  df_raw <- get(nm, inherits = TRUE)
  df_fish <- compute_moving_seconds_by_fish(df_raw, df_name = nm)
  
  df_fish_filt <- df_fish %>%
    mutate(timepoint_chr = as.character(timepoint)) %>%
    filter(timepoint_chr != timepoint_to_drop) %>%
    mutate(
      Group = factor(Group),
      Treatment = factor(Treatment),
      timepoint = factor(timepoint_chr, levels = sort(unique(timepoint_chr)))
    ) %>%
    select(-timepoint_chr)
  
  if (nrow(df_fish_filt) == 0) next
  
  trim_res <- trim_outliers_percentile(df_fish_filt, "moving_seconds", p_lo, p_hi)
  df_trim <- trim_res$df_trim
  removed_n <- trim_res$removed_n
  
  if (nrow(df_trim) == 0) next
  
  y_max_raw <- max(df_trim$moving_seconds, na.rm = TRUE)
  y_max <- if (is.finite(y_max_raw) && y_max_raw > 0) y_max_raw * 1.10 else 1
  y_breaks <- pretty(c(0, y_max), n = 10)
  
  fish_summaries[[nm]] <- df_trim %>%
    select(fish_id, moving_seconds, trial_id, timepoint, Group, Treatment, trial_date)
  
  title_A <- paste0(nm, " — Graph A (moving seconds)\nFiltered: timepoint != 1 | Outliers removed: ", removed_n)
  title_B <- paste0(nm, " — Graph B (moving seconds)\nFiltered: timepoint != 1 | Outliers removed: ", removed_n)
  
  pA <- make_graph_A(df_trim, title_A, y_max, y_breaks)
  pB <- make_graph_B(df_trim, title_B, y_max, y_breaks)
  
  ppt <- add_plot_slide(ppt, pA, paste0(nm, " — Graph A"))
  ppt <- add_plot_slide(ppt, pB, paste0(nm, " — Graph B"))
}

print(ppt, target = ppt_file)
ts_msg("Saved PPT: ", ppt_file)


# =============================================================================
# ==== 8) ZONE OCCUPANCY (SECONDS) -> PPT =====================================
# =============================================================================

ts_msg("Computing zone occupancy (seconds) and exporting to PPT...")

zone_by_frame <- get("zone_by_frame", envir = .GlobalEnv)

zone_cols_needed <- c(
  "session_id","trial_id","timepoint","fish_id","fish_in_session",
  "frame","time","zone","Group","Treatment","trial_date"
)
require_cols(zone_by_frame, zone_cols_needed, "zone_by_frame")

Y_MAX   <- 600
Y_BREAK <- 60
zone_ppt_file <- file.path(STEP_OUTPUT_DIR, paste0("occupancy_by_zone_time_", timestamp_tag(), ".pptx"))

z <- zone_by_frame %>%
  mutate(
    session_id = trimws(as.character(session_id)),
    trial_id   = suppressWarnings(as.integer(trial_id)),
    timepoint  = suppressWarnings(as.integer(timepoint)),
    fish_in_session = suppressWarnings(as.integer(fish_in_session)),
    fish_id    = trimws(as.character(fish_id)),
    frame      = suppressWarnings(as.integer(frame)),
    time       = suppressWarnings(as.numeric(time)),
    zone_chr   = tolower(trimws(as.character(zone))),
    zone_norm  = case_when(
      zone_chr == "flow" ~ "flow",
      zone_chr == "calm" ~ "calm",
      TRUE ~ zone_chr
    )
  ) %>%
  filter(zone_norm %in% c("flow","calm")) %>%
  arrange(session_id, trial_id, timepoint, fish_in_session, fish_id, time, frame)

if (nrow(z) == 0) stop("After filtering to (flow/calm), 0 rows remain.")

z_dt <- z %>%
  group_by(session_id, trial_id, timepoint, fish_in_session, fish_id) %>%
  arrange(time, frame, .by_group = TRUE) %>%
  mutate(
    dt = lead(time) - time,
    dt = ifelse(is.na(dt) | dt <= 0, median(dt[is.finite(dt) & dt > 0], na.rm = TRUE), dt),
    dt = ifelse(!is.finite(dt) | dt <= 0, NA_real_, dt)
  ) %>%
  ungroup()

zone_time_long <- z_dt %>%
  group_by(
    session_id, trial_id, timepoint, fish_in_session, fish_id,
    Group, Treatment, trial_date, zone_norm
  ) %>%
  summarise(time_seconds = sum(dt, na.rm = TRUE), .groups = "drop") %>%
  rename(zone = zone_norm)

zone_time_by_fish <- zone_time_long %>%
  pivot_wider(names_from = zone, values_from = time_seconds, values_fill = 0)

make_stacked_x <- function(df) {
  df %>%
    mutate(
      Group     = as.character(Group),
      Treatment = as.character(Treatment),
      timepoint = as.character(timepoint),
      x_stack   = factor(paste(Group, timepoint, Treatment, sep = "\n"))
    )
}

plot_zone_seconds <- function(df_long, zone_name, letter_tag) {
  d <- df_long %>%
    filter(zone == zone_name) %>%
    make_stacked_x() %>%
    mutate(time_seconds = pmin(time_seconds, Y_MAX))
  
  ggplot(d, aes(x = x_stack, y = time_seconds, color = Treatment)) +
    geom_point(position = position_jitter(width = 0.14, height = 0), size = 2.2, alpha = 0.85) +
    stat_summary(fun = mean, geom = "col", alpha = 0.25, width = 0.60) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18, linewidth = 0.7) +
    scale_y_continuous(
      limits = c(0, Y_MAX),
      breaks = seq(0, Y_MAX, by = Y_BREAK),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(x = NULL, y = "Time in zone (seconds)", title = str_to_title(zone_name)) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(size = 11, face = "bold", lineheight = 0.95),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      legend.position = "bottom",
      plot.margin = margin(14, 14, 14, 14)
    ) +
    annotate("text", x = -Inf, y = Inf, label = letter_tag,
             hjust = -0.35, vjust = 1.25, size = 6, fontface = "bold")
}

pA_flow <- plot_zone_seconds(zone_time_long, "flow", "A")
pB_calm <- plot_zone_seconds(zone_time_long, "calm", "B")
p_side  <- pA_flow + pB_calm + plot_layout(ncol = 2)

ppt2 <- read_pptx() %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with("Occupancy by zone (time in seconds)", location = ph_location_type(type = "title")) %>%
  ph_with(dml(ggobj = p_side), location = ph_location_type(type = "body"))

print(ppt2, target = zone_ppt_file)
ts_msg("Saved PPT: ", zone_ppt_file)

# =============================================================================
# ==== 9) QUICK INSPECTION (OPTIONAL) =========================================
# =============================================================================
# colnames(master_fish_by_frame)
# colnames(position_by_frame)
# colnames(speed_by_frame)
# if (exists("group_by_frame")) colnames(group_by_frame)
# colnames(zone_by_frame)


# =============================================================================
# ==== 10) EXPORT CORE TABLES TO STEP OUTPUT ==================================
# =============================================================================

ts_msg("Exporting core tables to STEP output directory...")

write_csv_if_exists <- function(obj_name, file_name) {
  if (!exists(obj_name, envir = .GlobalEnv)) {
    ts_msg("Skipping export (object not found): ", obj_name)
    return(invisible(NULL))
  }
  
  obj <- get(obj_name, envir = .GlobalEnv)
  if (!is.data.frame(obj)) {
    ts_msg("Skipping export (not a data.frame): ", obj_name)
    return(invisible(NULL))
  }
  
  out_path <- file.path(STEP_OUTPUT_DIR, file_name)
  readr::write_csv(obj, out_path)
  ts_msg("Saved: ", out_path)
}

write_csv_if_exists("master_fish_by_frame", "master_fish_by_frame.csv")
write_csv_if_exists("position_by_frame", "position_by_frame.csv")
write_csv_if_exists("speed_by_frame", "speed_by_frame.csv")
write_csv_if_exists("zone_by_frame", "zone_by_frame.csv")
write_csv_if_exists("group_by_frame", "group_by_frame.csv")
write_csv_if_exists("bodylength_by_fish", "bodylength_by_fish.csv")
write_csv_if_exists("bodylength_summary_by_session", "bodylength_summary_by_session.csv")
write_csv_if_exists("activity_TS0.5", "activity_TS0.5.csv")
write_csv_if_exists("activity_TS1", "activity_TS1.csv")
write_csv_if_exists("activity_TS1.5", "activity_TS1.5.csv")
write_csv_if_exists("activity_TS2", "activity_TS2.csv")
write_csv_if_exists("distance_by_fish", "distance_by_fish.csv")

# ---- R2: re-save enriched master_fish_by_frame to canonical store -----------
# The version of master_fish_by_frame in .GlobalEnv at this point has been
# enriched with locomotor-status columns + body-length speed (speed_BL_s).
# Writing it back to master_outputs/by_frame/master_fish_by_frame.csv ensures
# the canonical persistent store is the enriched version, so that re-running
# from STEP 4 in a fresh R session reads the right thing (instead of the
# pre-enriched STEP 1 output).
#
# Trade-off: overwrites the pristine STEP 1 output. STEP 1 output is
# reproducible by re-running STEP 1, so this is acceptable. Adds ~10-30 s of
# disk-write time per STEP 2 run.
.canonical_by_frame_dir  <- file.path(PIPELINE_DATA_DIR, "master_outputs", "by_frame")
.canonical_by_frame_path <- file.path(.canonical_by_frame_dir,
                                       "master_fish_by_frame.csv")
if (exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE) &&
    is.data.frame(master_fish_by_frame) && nrow(master_fish_by_frame) > 0L) {
  dir.create(.canonical_by_frame_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(master_fish_by_frame, .canonical_by_frame_path)
  ts_msg("R2: re-saved enriched master_fish_by_frame to canonical: ",
         .canonical_by_frame_path)
} else {
  ts_msg("R2: master_fish_by_frame not present in .GlobalEnv; canonical re-save skipped.")
}
rm(.canonical_by_frame_dir, .canonical_by_frame_path)

.got_write_run_settings(
  STEP_OUTPUT_DIR,
  TRIAL_SUMMARY_PATH    = TRIAL_SUMMARY_PATH,
  FORCE_RELOAD_FROM_DISK = FORCE_RELOAD_FROM_DISK,
  LENGTH_UNIT_PX_PER_CM  = LENGTH_UNIT_PX_PER_CM
)
.got_write_session_info(STEP_OUTPUT_DIR)

ts_msg("All exports complete.")