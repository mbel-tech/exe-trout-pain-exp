###############################################################################
##  1) LIBRARIES                                                             ##
###############################################################################
library(dplyr)
library(tidyr)
library(readr)
library(jsonlite)
library(stringr)
library(purrr)
library(tibble)
library(data.table)
library(swaRmverse)
library(stats)

###############################################################################
##  1b) CONFIG                                                               ##
###############################################################################
# ---- Conda / Python bboxes generator ----
# extract_bounding_boxes.py reads each session's preprocessing/list_of_blobs.pickle and
# writes <session_id>_bboxes.csv into trajectories/trajectories_csv/. STEP 08
# derives body_length_cm from it, so this is the origin of every body-length
# number in the pipeline. It is normally invisible because
# ensure_bboxes_for_session() short-circuits whenever a cached CSV is present,
# and all 70 sessions are currently cached — a broken BBOX_SCRIPT path therefore
# stayed silent for a long time (it was corrected here in 2026-08-20).
#
# Verified working end-to-end 2026-08-21 on R 4.6.0: conda env "idtrackerai"
# (python 3.10.19, numpy 2.2.6, idtrackerai.list_of_blobs importable), invoked
# through the same system2() path ensure_bboxes_for_session() builds. Regenerated
# output is md5-identical to the cached CSVs, so the cache and the generator
# agree. To force regeneration for one session, call
# ensure_bboxes_for_session(..., force = TRUE).
CONDA_EXE   <- Sys.getenv("CONDA_EXE", "conda")
CONDA_ENV   <- "idtrackerai"
BBOX_SCRIPT <- "01_tracking/extract_bounding_boxes.py"

DEFAULT_N_FISH_MAX      <- 4L
DEFAULT_MAX_GAP_SECONDS <- 1
DEFAULT_SPEED_PROB_MIN  <- NULL
DEFAULT_USE_TIME_DT     <- TRUE
DEFAULT_PROB_THRESHOLD  <- 0.80
DEFAULT_Y0_PX           <- mean(c(483, 477))
DEFAULT_FLOW_HIGHER_Y   <- TRUE

## [C3] Verbosity toggle — set FALSE to suppress per-session messages
VERBOSE            <- TRUE
## [A1] Timeout for conda bbox generation (seconds); requires R.utils if >0
BBOX_TIMEOUT_S     <- 600L

# ---- Bbox identity correction policy ----
ALLOW_DROP_SPURIOUS_BBOX_ID_ZERO <- TRUE

# ---- Jump detection (RETIRED 2026-08-20 -- see STEP 3.5 / OUTLIER_REJECT_*) --
# This mean+SD threshold is pooled across ALL tracks, so a genuine tracker
# runaway (thousands of px/frame) inflates the SD enough to raise the bar
# above itself -- self-defeating on exactly the frames it exists to catch.
# It also ran ~700 lines after speed_cm_s was computed, on the assembled
# master_fish_by_frame, so any correction it did make never reached the
# kinematics columns (verified: 47/47 corrected frames kept stale speeds up
# to 1,008,298 cm/s in run STEP2_output_20260524_160853). Superseded by
# STEP 3.5 (dt_reject_outliers_inplace), which runs per-session BEFORE
# interpolation and BEFORE speed is computed, using a robust (MAD-based)
# threshold that outliers cannot inflate. Kept here, disabled, for reference.
JUMP_DETECTION_ENABLED <- FALSE         # set FALSE to disable entirely
JUMP_THRESHOLD_METHOD  <- "sd_multiple" # "sd_multiple" | "percentile"
JUMP_SD_MULT           <- 10.0          # multiplier (idtracker.ai GUI default)
JUMP_PCT               <- 99L           # percentile for "percentile" method
JUMP_PCT_MULT          <- 2.0           # multiplier for "percentile" method
JUMP_INTERP            <- TRUE          # TRUE: interpolate over jumps; FALSE: leave as NA

# ---- STEP 3.5: outlier rejection (replaces jump detection, see above) -----
# Runs on raw x/y, per session, BEFORE dt_interpolate_short_gaps_inplace and
# BEFORE speed is computed -- so rejected points are genuinely removed, not
# just flagged after the fact, and the resulting gaps get filled by the
# existing (already-correct) short-gap interpolator rather than by ad hoc
# extrapolation. Two independent tests, OR'd together; either one nulls a
# point:
#   (a) out-of-frame: x/y outside [0-PAD, width/height+PAD] -- a physical
#       constraint from session.json/attributes.json, not a statistical one.
#   (b) implausible displacement: frame-to-frame speed (px/s, over
#       consecutive VALID points) exceeds max(median + K*MAD, FLOOR_BL_S *
#       median_body_length_px), pooled across the session's fish. MAD (not
#       SD) so a handful of extreme frames cannot inflate their own escape
#       threshold. Pooled per session (not per fish) because a per-fish MAD
#       collapses to ~0 for near-stationary fish and over-flags routine
#       jitter; the BL/s floor is what keeps pooling safe for quiet sessions.
# Iterates (OUTLIER_MAX_PASSES) because nulling the worst point in a jump can
# expose the next-worst as the new max -- a single pass cannot unwind a
# multi-frame excursion (see the "constant-speed ramp" defect in LOGS.md).
OUTLIER_REJECT_ENABLED <- TRUE
OUTLIER_FRAME_PAD_PX   <- 5    # tolerance outside the video frame, in px
OUTLIER_MAD_MULT       <- 30   # K; ~133-421 cm/s (6-19 BL/s) across this dataset's sessions
OUTLIER_FLOOR_BL_S     <- 8    # floor, body lengths/second
OUTLIER_MAX_PASSES     <- 5L

###############################################################################
##  2) UTILITIES                                                             ##
###############################################################################

## [C3] ts_msg now gated by VERBOSE global
ts_msg <- function(...) {
  if (!isTRUE(VERBOSE)) return(invisible(NULL))
  msg <- paste0(format(Sys.time(), "%H:%M:%S"), " | ", paste0(..., collapse = ""))
  message(msg)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

norm_path <- function(p, mustWork = FALSE) {
  normalizePath(p, winslash = "/", mustWork = mustWork)
}

is_valid_positive <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
}

clipboard_dir <- function() {
  if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
    p <- get("PIPELINE_DATA_DIR", envir = .GlobalEnv)
    return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  raw <- readLines("clipboard", warn = FALSE)
  raw <- trimws(raw)
  raw <- raw[nzchar(raw)]
  if (length(raw) == 0L) stop("Clipboard is empty or has no non-empty lines.")
  p <- raw[1L]
  p <- gsub('^["\']|["\']$', "", p)
  p <- gsub("\\\\", "/", p)
  if (!grepl("^[A-Za-z]:/|^//", p))
    stop("Clipboard path does not look like a drive path or UNC: ", p)
  p <- normalizePath(p, winslash = "/", mustWork = TRUE)
  if (!dir.exists(p)) stop("Clipboard path is not a directory: ", p)
  p
}

## [audit C1] n_max parameter was unused — now documented as ignored (header is always line 1)
read_csv_header_only <- function(path, n_max = 5) {
  hdr <- readLines(path, n = 1L, warn = FALSE)
  strsplit(trimws(hdr), ",")[[1L]]
}

internal_na_mask <- function(is_na_vec) {
  n <- length(is_na_vec)
  if (n == 0L) return(logical(0))
  mask <- rep(FALSE, n)
  r <- rle(is_na_vec)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  for (i in seq_along(r$lengths)) {
    if (r$values[i] && starts[i] > 1L && ends[i] < n) {
      mask[starts[i]:ends[i]] <- TRUE
    }
  }
  mask
}

## [C2] Portable na.locf replacement — eliminates zoo dependency
locf_fill <- function(x) {
  idx <- cummax(seq_along(x) * (!is.na(x)))
  idx[idx == 0L] <- NA_integer_
  x[idx]
}

# -----------------------------------------------------------------------------
# Bbox identity diagnostics + optional correction
# -----------------------------------------------------------------------------

read_bbox_identity_info <- function(bboxes_path) {
  if (is.na(bboxes_path) || !file.exists(bboxes_path)) {
    return(list(
      ok_read = FALSE,
      reason = "bboxes_missing",
      ids = integer(0),
      n_ids = NA_integer_
    ))
  }
  
  bb <- tryCatch(data.table::fread(bboxes_path, nrows = Inf), error = function(e) NULL)
  if (is.null(bb)) {
    return(list(
      ok_read = FALSE,
      reason = "bboxes_read_error",
      ids = integer(0),
      n_ids = NA_integer_
    ))
  }
  
  names(bb) <- tolower(names(bb))
  
  if ("identity" %in% names(bb) && !"individual" %in% names(bb)) {
    data.table::setnames(bb, "identity", "individual")
  }
  
  if (!"individual" %in% names(bb)) {
    return(list(
      ok_read = FALSE,
      reason = "bboxes_no_identity_col",
      ids = integer(0),
      n_ids = NA_integer_
    ))
  }
  
  ids <- sort(unique(as.integer(bb$individual)))
  ids <- ids[!is.na(ids)]
  
  list(
    ok_read = TRUE,
    reason = "ok",
    ids = ids,
    n_ids = length(ids)
  )
}

diagnose_bbox_id_set <- function(ids, n_fish, allow_drop_zero = TRUE) {
  ids <- sort(unique(as.integer(ids)))
  ids <- ids[!is.na(ids)]
  expected <- seq_len(n_fish)
  
  # perfect match
  if (identical(ids, expected)) {
    return(list(
      mismatch = FALSE,
      reason = "ok",
      ids = ids,
      ids_used = ids,
      drop_ids = integer(0),
      autocorrected = FALSE
    ))
  }
  
  # allow extra zero only
  if (isTRUE(allow_drop_zero)) {
    ids_wo_zero <- ids[ids != 0L]
    if (length(ids) == n_fish + 1L &&
        0L %in% ids &&
        identical(ids_wo_zero, expected)) {
      return(list(
        mismatch = FALSE,
        reason = sprintf(
          "autocorrect_drop_zero(expected=%s; raw_ids=%s)",
          paste(expected, collapse = ","),
          paste(ids, collapse = ",")
        ),
        ids = ids,
        ids_used = ids_wo_zero,
        drop_ids = 0L,
        autocorrected = TRUE
      ))
    }
  }
  
  # count mismatch
  if (length(ids) != n_fish) {
    return(list(
      mismatch = TRUE,
      reason = sprintf(
        "count_mismatch(n_fish=%d, n_ids=%d, ids=%s)",
        n_fish, length(ids), paste(ids, collapse = ",")
      ),
      ids = ids,
      ids_used = ids,
      drop_ids = integer(0),
      autocorrected = FALSE
    ))
  }
  
  # same count, wrong labels
  return(list(
    mismatch = TRUE,
    reason = sprintf(
      "label_mismatch(expected=%s; got=%s)",
      paste(expected, collapse = ","),
      paste(ids, collapse = ",")
    ),
    ids = ids,
    ids_used = ids,
    drop_ids = integer(0),
    autocorrected = FALSE
  ))
}

###############################################################################
## UPDATED FUNCTIONS: add support for session.json (idtracker.ai) + robust file
## discovery that only triggers fallback search if the canonical path fails.
##
## Drop these into your script (replacing the existing versions).
###############################################################################

## [audit C2] Fallback for headerless CSVs: count columns / 2
detect_n_fish_from_trajectories <- function(trajectories_path) {
  hdr <- read_csv_header_only(trajectories_path)
  x_cols <- grep("^x\\d+$", hdr, value = TRUE)
  if (length(x_cols) > 0L) {
    nums <- as.integer(sub("^x", "", x_cols))
    return(max(nums, na.rm = TRUE))
  }
  ## Fallback: headerless numeric CSV. Count columns.
  ## All tokens are numeric -> no named headers
  all_numeric <- all(grepl("^-?[0-9]", trimws(hdr)))
  if (all_numeric) {
    nc <- length(hdr)
    ## Odd column count -> first col is time
    n_data <- if (nc %% 2L == 1L) nc - 1L else nc
    n_fish <- n_data %/% 2L
    if (n_fish >= 1L) return(as.integer(n_fish))
  }
  NA_integer_
}

# -----------------------------------------------------------------------------
# 1) Robust file resolver (canonical first, recursive fallback second)
# -----------------------------------------------------------------------------

.choose_best_hit <- function(paths, prefer_dir_regex = "(?i)trajectories_csv") {
  if (length(paths) == 0L) return(NA_character_)
  paths <- norm_path(paths, mustWork = FALSE)
  
  ## Prefer hits under trajectories_csv, then shorter paths
  score <- rep(0L, length(paths))
  score <- score + ifelse(grepl(prefer_dir_regex, paths, perl = TRUE), 100L, 0L)
  score <- score - nchar(paths)
  paths[which.max(score)]
}

.find_in_session_recursive <- function(session_dir, filename_pattern, max_hits = 200L) {
  if (!dir.exists(session_dir)) return(character())
  hits <- list.files(
    session_dir,
    pattern = filename_pattern,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  if (length(hits) > max_hits) hits <- hits[seq_len(max_hits)]
  hits[file.exists(hits)]
}

resolve_session_file <- function(session_dir,
                                 canonical_rel = NULL,
                                 fallback_pattern = NULL,
                                 prefer_dir_regex = "(?i)trajectories_csv") {
  if (!dir.exists(session_dir)) return(NA_character_)
  
  ## --- canonical fast path ---
  if (!is.null(canonical_rel)) {
    canon_full <- file.path(session_dir, canonical_rel)
    if (file.exists(canon_full)) return(norm_path(canon_full))
  }
  
  ## --- fallback recursive search (ONLY if canonical missing) ---
  if (!is.null(fallback_pattern)) {
    hits <- .find_in_session_recursive(session_dir, fallback_pattern)
    best <- .choose_best_hit(hits, prefer_dir_regex = prefer_dir_regex)
    if (!is.na(best) && file.exists(best)) return(norm_path(best))
  }
  
  NA_character_
}


# -----------------------------------------------------------------------------
# 2) Index sessions (updated: includes session_json_path and robust discovery)
# -----------------------------------------------------------------------------

index_session_files_smart <- function(base_dir, session_dir_regex = "(?i)^session_") {
  all_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  keep <- all_dirs[grepl(session_dir_regex, basename(all_dirs), perl = TRUE)]
  if (length(keep) == 0L) return(tibble(
    session_id = character(), session_dir = character(),
    trajectories_path = character(), attributes_path = character(),
    session_json_path = character(),
    areas_path = character(), bboxes_path = character()
  ))
  
  out <- tibble(
    session_id  = basename(keep),
    session_dir = norm_path(keep)
  ) %>%
    rowwise() %>%
    mutate(
      ## Canonical: trajectories/trajectories_csv/trajectories.csv
      trajectories_path = resolve_session_file(
        session_dir,
        canonical_rel = "trajectories/trajectories_csv/trajectories.csv",
        fallback_pattern = "(?i)^trajectories(.*)?\\.csv$",
        prefer_dir_regex = "(?i)trajectories_csv"
      ),
      
      ## Canonical: trajectories/trajectories_csv/attributes.json
      attributes_path = resolve_session_file(
        session_dir,
        canonical_rel = "trajectories/trajectories_csv/attributes.json",
        fallback_pattern = "(?i)^attributes(.*)?\\.json$",
        prefer_dir_regex = "(?i)trajectories_csv"
      ),
      
      ## NEW: session.json at session root (idtracker.ai), but allow fallback
      session_json_path = resolve_session_file(
        session_dir,
        canonical_rel = "session.json",
        fallback_pattern = "(?i)^session\\.json$",
        prefer_dir_regex = "(?i)$"  # no strong preference needed
      ),
      
      ## Optional: areas.csv (keep canonical, fallback if you want)
      areas_path = resolve_session_file(
        session_dir,
        canonical_rel = "trajectories/trajectories_csv/areas.csv",
        fallback_pattern = "(?i)^areas(.*)?\\.csv$",
        prefer_dir_regex = "(?i)trajectories_csv"
      ),
      
      ## Canonical expected bboxes name based on session_id
      bboxes_path = resolve_session_file(
        session_dir,
        canonical_rel = file.path("trajectories/trajectories_csv", paste0(session_id, "_bboxes.csv")),
        fallback_pattern = "(?i)_bboxes\\.csv$",
        prefer_dir_regex = "(?i)trajectories_csv"
      )
    ) %>%
    ungroup() %>%
    arrange(session_id)
  
  out
}


# -----------------------------------------------------------------------------
# 3) Bboxes discovery (updated: canonical first; fallback only if missing)
# -----------------------------------------------------------------------------

find_bboxes_in_session <- function(session_dir) {
  sid <- basename(session_dir)

  p <- resolve_session_file(
    session_dir,
    canonical_rel = file.path("trajectories", "trajectories_csv", paste0(sid, "_bboxes.csv")),
    fallback_pattern = "(?i)_bboxes\\.csv$",
    prefer_dir_regex = "(?i)trajectories_csv"
  )

  if (!is.na(p) && file.exists(p)) return(p)
  NA_character_
}

ensure_bboxes_for_session <- function(session_dir, conda_exe, conda_env, script_path,
                                      log_dir = NULL, timeout_s = BBOX_TIMEOUT_S, force = FALSE) {
  sid <- basename(session_dir)
  out_dir_bb <- file.path(session_dir, "trajectories", "trajectories_csv")
  ensure_dir(out_dir_bb)
  expected <- file.path(out_dir_bb, paste0(sid, "_bboxes.csv"))

  # [OPT-7] Check if output already exists; skip extraction if found
  if (file.exists(expected) && !isTRUE(force)) {
    cat(sprintf("  [cached] bboxes already exist: %s\n", basename(expected)))
    return(norm_path(expected))
  }

  ## force=TRUE must actually force: only accept an existing cached file
  ## found by fallback search when we are NOT forcing regeneration.
  if (!isTRUE(force)) {
    found <- find_bboxes_in_session(session_dir)
    if (!is.na(found) && file.exists(found)) return(found)
  }

  ts_msg("Running conda bbox script for: ", sid)

  ## [A1] Build command depending on .bat vs .exe
  sd_q  <- shQuote(norm_path(session_dir))
  od_q  <- shQuote(norm_path(out_dir_bb))
  sp_q  <- shQuote(norm_path(script_path))
  conda_args <- c("run", "-n", conda_env, "python", sp_q,
                   "--session_dir", sd_q, "--out_dir", od_q)

  is_bat <- grepl("\\.bat$", conda_exe, ignore.case = TRUE)
  if (is_bat) {
    exe  <- "cmd.exe"
    args <- c("/c", shQuote(norm_path(conda_exe)), conda_args)
  } else {
    exe  <- norm_path(conda_exe)
    args <- conda_args
  }

  ## [A2] Log to per-session file
  log_file <- NULL
  if (!is.null(log_dir) && dir.exists(log_dir)) {
    log_file <- file.path(log_dir, paste0("bboxes_", sid, ".log"))
  }

  run_system2 <- function() {
    if (!is.null(log_file)) {
      system2(exe, args = args, stdout = log_file, stderr = log_file)
    } else {
      system2(exe, args = args, stdout = TRUE, stderr = TRUE)
    }
  }

  res <- NULL
  timed_out <- FALSE

  if (requireNamespace("R.utils", quietly = TRUE) && is.numeric(timeout_s) && timeout_s > 0) {
    res <- tryCatch(
      R.utils::withTimeout(run_system2(), timeout = timeout_s),
      TimeoutException = function(e) { timed_out <<- TRUE; NULL },
      error = function(e) { warning("system2 failed: ", e$message); NULL }
    )
  } else {
    res <- tryCatch(run_system2(),
                    error = function(e) { warning("system2 failed: ", e$message); NULL })
  }

  if (timed_out) {
    warning("Bbox script timed out (", timeout_s, "s) for ", sid)
    return(NA_character_)
  }

  if (is.integer(res) && length(res) == 1L) {
    status <- res
  } else {
    status <- attr(res, "status")
  }

  if (!is.null(status) && status != 0L) {
    if (!is.null(log_file) && file.exists(log_file)) {
      log_lines <- readLines(log_file, warn = FALSE)
      tail_lines <- utils::tail(log_lines, 50)
      warning("Bbox script failed for ", sid, " (see ", log_file, "):\n",
              paste(tail_lines, collapse = "\n"))
    } else if (is.character(res)) {
      warning("Bbox script failed for ", sid, ":\n",
              paste(utils::tail(res, 50), collapse = "\n"))
    }
    return(NA_character_)
  }

  if (!file.exists(expected)) {
    warning("Bboxes not found after running script for ", sid)
    NA_character_
  } else {
    norm_path(expected)
  }
}


# -----------------------------------------------------------------------------
# 4) Read attributes with session.json fallback
#    - Primary: attributes.json (trajectories/trajectories_csv/attributes.json)
#    - Fallback: session.json (session root)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# session.json SAFE reader (idtracker.ai): handles Infinity/NaN tokens
# -----------------------------------------------------------------------------
read_session_json_safe <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  
  txt <- tryCatch(paste(readLines(path, warn = FALSE), collapse = "\n"),
                  error = function(e) NULL)
  if (is.null(txt) || !nzchar(txt)) return(NULL)
  
  # Replace invalid JSON tokens used by idtracker.ai (Infinity / -Infinity / NaN)
  # Use `null` so jsonlite parses and you still have the keys present.
  txt2 <- txt
  txt2 <- gsub("(?<![A-Za-z0-9_\"'])-?Infinity(?![A-Za-z0-9_\"'])", "null", txt2, perl = TRUE)
  txt2 <- gsub("(?<![A-Za-z0-9_\"'])NaN(?![A-Za-z0-9_\"'])", "null", txt2, perl = TRUE)
  
  obj <- tryCatch(
    jsonlite::fromJSON(txt2, simplifyVector = TRUE),
    error = function(e) NULL
  )
  
  obj
}

read_session_attributes <- function(attributes_path = NA_character_,
                                    session_json_path = NA_character_) {
  raw_attr <- NULL
  raw_sess <- NULL
  
  # Try attributes.json first (normal route)
  if (!is.na(attributes_path) && file.exists(attributes_path)) {
    raw_attr <- tryCatch(jsonlite::fromJSON(attributes_path, simplifyVector = TRUE),
                         error = function(e) NULL)
  }
  
  # session.json (fallback or extra metadata)
  if (!is.na(session_json_path) && file.exists(session_json_path)) {
    raw_sess <- read_session_json_safe(session_json_path)
  }
  
  raw <- if (!is.null(raw_attr)) raw_attr else raw_sess
  if (is.null(raw)) {
    stop("Could not read attributes.json or session.json for this session.")
  }
  
  # fps: in both attributes.json and session.json (your example uses frames_per_second)
  fps <- raw$frames_per_second
  if (is.null(fps) && !is.null(raw$frames_per_second)) fps <- raw$frames_per_second
  if (is.null(fps)) fps <- NA_real_
  
  # length_unit: usually ONLY in attributes.json
  lu <- raw$length_unit
  if (is.null(lu)) lu <- NA_real_

  # ---- Frame geometry + median body length, for the STEP 3.5 outlier gate ----
  # width/height appear under the same key in both attributes.json and
  # session.json. Body length differs: attributes.json calls it "body_length",
  # session.json calls it "median_body_length" — try both, preferring whichever
  # `raw` resolved to (raw_attr if present, else raw_sess), then fall back to
  # explicitly checking the other source so a partial attributes.json doesn't
  # lose geometry that session.json has.
  pick_num <- function(...) {
    for (v in list(...)) {
      if (!is.null(v) && is.numeric(v) && length(v) == 1L && is.finite(v)) return(as.numeric(v))
    }
    NA_real_
  }
  frame_width  <- pick_num(raw$width,  raw_attr$width,  raw_sess$width)
  frame_height <- pick_num(raw$height, raw_attr$height, raw_sess$height)
  median_bl_px <- pick_num(raw$body_length, raw$median_body_length,
                           raw_attr$body_length, raw_sess$median_body_length)

  list(
    fps                 = as.numeric(fps),
    length_unit         = as.numeric(lu),
    frame_width         = frame_width,
    frame_height        = frame_height,
    median_body_length_px = median_bl_px,
    raw_attributes      = raw_attr,
    raw_session         = raw_sess,
    raw                 = raw
  )
}
# -----------------------------------------------------------------------------
# 5) Read session wide -> long (updated: passes session_json_path into attribute reader)
# -----------------------------------------------------------------------------

read_session_wide <- function(session_id, files_tbl, length_unit_fallback = NA) {
  ## Single-row lookup
  row <- files_tbl[files_tbl$session_id == session_id, ]
  stopifnot(nrow(row) >= 1L)
  if (nrow(row) > 1L) row <- row[1L, ]
  
  traj_path <- row$trajectories_path
  attr_path <- row$attributes_path
  sess_json <- if ("session_json_path" %in% names(row)) row$session_json_path else NA_character_
  bbox_path <- row$bboxes_path
  
  if (is.na(traj_path) || !file.exists(traj_path))
    stop("Trajectories file missing for ", session_id)
  
  ## Allow attribute source to be either attributes.json or session.json
  if ((is.na(attr_path) || !file.exists(attr_path)) && (is.na(sess_json) || !file.exists(sess_json))) {
    stop("Neither attributes.json nor session.json found for ", session_id)
  }
  
  att <- read_session_attributes(attributes_path = attr_path, session_json_path = sess_json)
  fps <- att$fps
  if (!is_valid_positive(fps)) stop("Invalid fps for ", session_id)
  
  length_unit <- att$length_unit
  fallback_used <- FALSE
  if (!is_valid_positive(length_unit)) {
    if (is_valid_positive(length_unit_fallback)) {
      length_unit <- length_unit_fallback
      fallback_used <- TRUE
      ts_msg("Using fallback length_unit=", round(length_unit, 4), " for ", session_id)
    } else {
      stop("Invalid length_unit and no valid fallback for ", session_id,
           " (note: session.json usually does not include length_unit).")
    }
  }
  
  ## Detect x/y/probability columns from header
  hdr <- read_csv_header_only(traj_path)
  x_cols <- grep("^x\\d+$", hdr, value = TRUE)
  
  na_strings <- c("", "NA", "NaN", "nan", "None", "null", "-")
  
  wide <- tryCatch({
    data.table::fread(traj_path, na.strings = na_strings)
  }, error = function(e) NULL)
  
  if (is.null(wide)) {
    ts_msg("fread failed for ", session_id, ", falling back to readr::read_csv")
    wide <- readr::read_csv(
      traj_path,
      col_types = cols(.default = col_double()),
      na = na_strings,
      show_col_types = FALSE
    )
    wide <- as.data.table(wide)
  }
  
  cn <- names(wide)
  
  if (length(x_cols) > 0) {
    n_fish <- length(x_cols)
    nums <- sort(as.integer(sub("^x", "", x_cols)))
    
    has_time <- "time" %in% cn
    sel_cols <- c()
    if (has_time) sel_cols <- "time"
    
    for (i in nums) {
      sel_cols <- c(sel_cols, paste0("x", i), paste0("y", i))
      
      ## optional probability columns
      p_i <- c(
        paste0("p", i),
        paste0("prob", i),
        paste0("probability", i)
      )
      sel_cols <- c(sel_cols, p_i[p_i %in% cn])
    }
    
    wide <- wide[, .SD, .SDcols = intersect(sel_cols, cn)]
    n_frames <- nrow(wide)
    
    if (has_time) {
      time_vec <- as.numeric(wide$time)
    } else {
      time_vec <- seq(0, (n_frames - 1)) / fps
    }
    
    for (cc in names(wide)) {
      if (is.numeric(wide[[cc]])) {
        vals <- wide[[cc]]
        vals[is.nan(vals) | is.infinite(vals)] <- NA_real_
        data.table::set(wide, j = cc, value = vals)
      } else {
        data.table::set(wide, j = cc, value = as.numeric(wide[[cc]]))
      }
    }
    
    rows_list <- vector("list", n_fish)
    for (i in seq_len(n_fish)) {
      fi <- nums[i]
      xc <- paste0("x", fi)
      yc <- paste0("y", fi)
      
      pc <- c(paste0("p", fi), paste0("prob", fi), paste0("probability", fi))
      pc <- pc[pc %in% names(wide)]
      prob_vec <- if (length(pc) > 0) wide[[pc[1]]] else rep(NA_real_, n_frames)
      
      rows_list[[i]] <- data.table(
        frame       = 0L:(n_frames - 1L),
        individual  = as.integer(i),
        time        = time_vec,
        x           = wide[[xc]],
        y           = wide[[yc]],
        probability = prob_vec
      )
    }
    long <- rbindlist(rows_list)
    
  } else {
    nc <- ncol(wide)
    first_col <- as.numeric(wide[[1]])
    
    starts_near_zero <- !is.na(first_col[1]) && first_col[1] < 1.0
    is_monotonic <- !any(is.na(first_col)) && all(diff(first_col) >= 0, na.rm = TRUE)
    has_time <- (nc %% 2L == 1L) && is_monotonic && starts_near_zero
    
    if (has_time) {
      time_vec <- first_col
      data_cols <- 2:nc
    } else {
      data_cols <- 1:nc
      time_vec <- seq(0, nrow(wide) - 1) / fps
    }
    
    n_data_cols <- length(data_cols)
    n_fish <- n_data_cols %/% 2L
    if (n_fish < 1L) stop("Cannot detect fish columns for ", session_id)
    
    n_frames <- nrow(wide)
    
    mat <- as.matrix(wide[, .SD, .SDcols = data_cols])
    storage.mode(mat) <- "double"
    mat[is.nan(mat) | is.infinite(mat)] <- NA_real_
    
    rows_list <- vector("list", n_fish)
    for (i in seq_len(n_fish)) {
      ci_x <- (i - 1) * 2 + 1
      ci_y <- (i - 1) * 2 + 2
      rows_list[[i]] <- data.table(
        frame       = 0L:(n_frames - 1L),
        individual  = as.integer(i),
        time        = time_vec,
        x           = mat[, ci_x],
        y           = mat[, ci_y],
        probability = NA_real_
      )
    }
    long <- rbindlist(rows_list)
  }
  
  list(
    fps                        = fps,
    length_unit                = length_unit,
    length_unit_fallback_used  = fallback_used,
    n_fish                     = n_fish,
    trajectories_long          = long,
    bboxes_path                = bbox_path,
    trajectories_path          = traj_path,
    attributes_path            = attr_path,
    session_json_path          = sess_json,
    ## Frame geometry + median body length for the STEP 3.5 outlier gate.
    ## NA when neither attributes.json nor session.json carried them for this
    ## session; STEP 3.5 must degrade gracefully (speed rule only) in that case.
    frame_width                = att$frame_width,
    frame_height               = att$frame_height,
    median_body_length_px      = att$median_body_length_px
  )
}

######

base_dir <- clipboard_dir()
files_index <- index_session_files_smart(base_dir)

dplyr::glimpse(files_index)

# See which sessions are missing attributes.json but have session.json
files_index %>%
  dplyr::transmute(
    session_id,
    has_traj   = !is.na(trajectories_path),
    has_attr   = !is.na(attributes_path),
    has_sjson  = !is.na(session_json_path)
  ) %>%
  dplyr::count(has_traj, has_attr, has_sjson)


library(dplyr)
library(purrr)
library(jsonlite)

quick_sessionjson_audit <- function(files_index, n_show = 20) 
  {
  res <- files_index %>%
    mutate(
      sess = purrr::map(session_json_path, ~ read_session_json_safe(.x)),
      fps_sjson = purrr::map_dbl(sess, ~ if (is.null(.x)) NA_real_ else as.numeric(.x$frames_per_second)),
      n_animals_sjson = purrr::map_int(sess, ~ if (is.null(.x)) NA_integer_ else as.integer(.x$number_of_animals)),
      n_frames_sjson  = purrr::map_int(sess, ~ if (is.null(.x)) NA_integer_ else as.integer(.x$number_of_frames))
    ) %>%
    select(-sess)
  
  # detect n_fish from trajectories header (fast)
  res <- res %>%
    mutate(
      n_fish_detected = purrr::map_int(trajectories_path, ~{
        if (is.na(.x) || !file.exists(.x)) return(NA_integer_)
        detect_n_fish_from_trajectories(.x)
      })
    )
  
  # flags
  res <- res %>%
    mutate(
      flag_animals_mismatch = !is.na(n_animals_sjson) & !is.na(n_fish_detected) & (n_animals_sjson != n_fish_detected),
      flag_fps_missing      = is.na(fps_sjson),
      flag_frames_missing   = is.na(n_frames_sjson)
    )
  
  bad <- res %>% filter(flag_animals_mismatch | flag_fps_missing | flag_frames_missing)
  
  message("Sessions with any session.json flags: ", nrow(bad))
  
  if (nrow(bad) > 0) {
    print(
      bad %>%
        select(session_id, n_animals_sjson, n_fish_detected, fps_sjson, n_frames_sjson,
               flag_animals_mismatch, flag_fps_missing, flag_frames_missing) %>%
        head(n_show),
      n = n_show
    )
  }
  
  invisible(res)
}

audit <- quick_sessionjson_audit(files_index)

###############################################################################
##  5) SANITY CHECK: BBOX IDs VS N FISH                                     ##
###############################################################################

sanity_check_bboxes_vs_trajectories <- function(
    bboxes_path,
    trajectories_path,
    allow_drop_zero = ALLOW_DROP_SPURIOUS_BBOX_ID_ZERO
) {
  if (is.na(bboxes_path) || !file.exists(bboxes_path)) {
    return(list(
      mismatch = NA,
      reason = "bboxes_missing",
      n_fish = NA_integer_,
      n_ids = NA_integer_,
      ids = integer(0),
      ids_used = integer(0),
      drop_ids = integer(0),
      autocorrected = FALSE
    ))
  }
  
  if (is.na(trajectories_path) || !file.exists(trajectories_path)) {
    return(list(
      mismatch = NA,
      reason = "trajectories_missing",
      n_fish = NA_integer_,
      n_ids = NA_integer_,
      ids = integer(0),
      ids_used = integer(0),
      drop_ids = integer(0),
      autocorrected = FALSE
    ))
  }
  
  n_fish <- detect_n_fish_from_trajectories(trajectories_path)
  if (is.na(n_fish)) {
    return(list(
      mismatch = TRUE,
      reason = "trajectories_cannot_detect_n_fish",
      n_fish = NA_integer_,
      n_ids = NA_integer_,
      ids = integer(0),
      ids_used = integer(0),
      drop_ids = integer(0),
      autocorrected = FALSE
    ))
  }
  
  id_info <- read_bbox_identity_info(bboxes_path)
  if (!isTRUE(id_info$ok_read)) {
    return(list(
      mismatch = TRUE,
      reason = id_info$reason,
      n_fish = n_fish,
      n_ids = id_info$n_ids,
      ids = id_info$ids,
      ids_used = integer(0),
      drop_ids = integer(0),
      autocorrected = FALSE
    ))
  }
  
  diag <- diagnose_bbox_id_set(
    ids = id_info$ids,
    n_fish = n_fish,
    allow_drop_zero = allow_drop_zero
  )
  
  list(
    mismatch = diag$mismatch,
    reason = diag$reason,
    n_fish = n_fish,
    n_ids = id_info$n_ids,
    ids = diag$ids,
    ids_used = diag$ids_used,
    drop_ids = diag$drop_ids,
    autocorrected = diag$autocorrected
  )
}
###############################################################################
##  6) BODY LENGTH UTILITIES                                                 ##
###############################################################################

read_bboxes_safe <- function(path, drop_ids = integer(0)) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  
  df <- tryCatch(data.table::fread(path), error = function(e) NULL)
  if (is.null(df)) return(NULL)
  
  names(df) <- tolower(names(df))
  
  if ("identity" %in% names(df) && !"individual" %in% names(df)) {
    setnames(df, "identity", "individual")
  }
  
  required <- c("frame", "individual", "body_length_px")
  if (!all(required %in% names(df))) return(NULL)
  
  df[, frame := as.integer(frame)]
  df[, individual := as.integer(individual)]
  df[, body_length_px := as.numeric(body_length_px)]
  
  # optional drop of known spurious IDs
  if (length(drop_ids) > 0L) {
    df <- df[!(individual %in% drop_ids)]
  }
  
  df <- df[is.finite(body_length_px) & body_length_px > 0]
  if (nrow(df) == 0L) return(NULL)
  
  df <- df[, .SD, .SDcols = intersect(
    names(df),
    c("frame", "individual", "body_length_px", "x_min", "y_min", "x_max", "y_max")
  )]
  
  as.data.frame(df)
}

choose_bodylength_summary <- function(x, prefer = "median", trim = 0.10, min_n_for_trim = 200) {
  x <- x[is.finite(x) & x > 0]
  n <- length(x)
  if (n == 0L) return(list(method = NA_character_, value = NA_real_,
                            heavy_tail = NA, n = 0L))
  med <- median(x)
  iqr_val <- IQR(x)
  heavy_tail <- (iqr_val / med) > 0.35

  if (prefer == "median" || n < min_n_for_trim || heavy_tail) {
    list(method = "median", value = med, heavy_tail = heavy_tail, n = n)
  } else {
    list(method = "trimmed_mean", value = mean(x, trim = trim),
         heavy_tail = heavy_tail, n = n)
  }
}

flag_unstable_bodylength <- function(x, cv_thresh = 0.15, iqr_thresh = 0.20, min_n = 300) {
  x <- x[is.finite(x) & x > 0]
  n <- length(x)
  if (n == 0L) return(list(flag = TRUE, reason = "no_data", cv = NA_real_, iqr_med = NA_real_))
  mn <- mean(x)
  sd_x <- sd(x)
  med <- median(x)
  iqr_val <- IQR(x)
  cv <- if (mn > 0) sd_x / mn else NA_real_
  iqr_med <- if (med > 0) iqr_val / med else NA_real_

  if (n < min_n) {
    return(list(flag = TRUE, reason = sprintf("too_few_frames(n=%d)", n),
                cv = cv, iqr_med = iqr_med))
  }
  flagged <- (!is.na(cv) && cv > cv_thresh) || (!is.na(iqr_med) && iqr_med > iqr_thresh)
  reason <- if (flagged) {
    parts <- c()
    if (!is.na(cv) && cv > cv_thresh) parts <- c(parts, sprintf("cv=%.3f", cv))
    if (!is.na(iqr_med) && iqr_med > iqr_thresh) parts <- c(parts, sprintf("iqr_med=%.3f", iqr_med))
    paste(parts, collapse = "; ")
  } else "ok"

  list(flag = flagged, reason = reason, cv = cv, iqr_med = iqr_med)
}

## [audit C4] Replaced O(N^2) rowwise scan with pre-split list approach
compute_bodylength_by_fish <- function(dat_i, bboxes_df, length_unit) {
  dat_i <- as.data.frame(dat_i, stringsAsFactors = FALSE)
  bboxes_df <- as.data.frame(bboxes_df, stringsAsFactors = FALSE)
  
  joined <- dplyr::inner_join(
    dat_i[, c("session_id", "trial_id", "segment_id", "fish_id", "fish_in_session",
              "frame", "individual", "x_interp", "y_interp")],
    bboxes_df[, c("frame", "individual", "body_length_px")],
    by = c("frame", "individual")
  )
  
  if (nrow(joined) == 0L) return(NULL)
  
  joined$session_id      <- as.character(joined$session_id)
  joined$trial_id        <- as.integer(joined$trial_id)
  joined$segment_id      <- as.integer(joined$segment_id)
  joined$fish_id         <- as.character(joined$fish_id)
  joined$fish_in_session <- as.integer(joined$fish_in_session)
  joined$body_length_px  <- as.numeric(joined$body_length_px)
  
  bl_by_fish <- split(joined$body_length_px, joined$fish_id)
  
  summ <- joined %>%
    dplyr::group_by(fish_id) %>%
    dplyr::summarise(
      session_id        = session_id[1],
      trial_id          = trial_id[1],
      segment_id        = segment_id[1],
      fish_in_session   = fish_in_session[1],
      n_bbox_frames     = dplyr::n(),
      px_median         = stats::median(body_length_px, na.rm = TRUE),
      px_mean_trim      = mean(body_length_px, trim = 0.10, na.rm = TRUE),
      body_length_sd_px = stats::sd(body_length_px, na.rm = TRUE),
      .groups = "drop"
    )
  
  summ <- summ %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      bl_summary = list(choose_bodylength_summary(bl_by_fish[[fish_id]])),
      body_length_method     = bl_summary$method,
      body_length_px         = bl_summary$value,
      body_length_heavy_tail = bl_summary$heavy_tail,
      bl_unstable = list(flag_unstable_bodylength(bl_by_fish[[fish_id]])),
      body_length_unstable        = bl_unstable$flag,
      body_length_unstable_reason = bl_unstable$reason,
      body_length_cm    = body_length_px / length_unit,
      body_length_sd_cm = body_length_sd_px / length_unit
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-bl_summary, -bl_unstable)
  
  summ
}
make_body_length_summary_wide <- function(body_by_fish, n_fish_max = 4) {
  if (is.null(body_by_fish) || nrow(body_by_fish) == 0L) return(NULL)

  wide <- body_by_fish %>%
    filter(fish_in_session <= n_fish_max) %>%
    select(session_id, trial_id, segment_id, fish_in_session, body_length_cm) %>%
    pivot_wider(
      names_from  = fish_in_session,
      values_from = body_length_cm,
      names_prefix = "fish_"
    )

  for (i in seq_len(n_fish_max)) {
    col <- paste0("fish_", i)
    if (!col %in% names(wide)) wide[[col]] <- NA_real_
  }

  wide %>% select(session_id, trial_id, segment_id,
                  all_of(paste0("fish_", seq_len(n_fish_max))))
}

###############################################################################
##  7) PARSE TRIAL / SEGMENT                                                 ##
###############################################################################

parse_trial_segment <- function(session_id) {
  m <- str_match(session_id, "(?i)trial[_[:space:]]*(\\d+).*segment[_[:space:]]*(\\d+)")
  if (is.na(m[1, 1])) {
    tibble(trial_id = -1L, segment_id = -1L)
  } else {
    tibble(trial_id = as.integer(m[1, 2]), segment_id = as.integer(m[1, 3]))
  }
}

###############################################################################
##  8) HEAVY CORE: data.table INTERPOLATION (IN-PLACE)                       ##
###############################################################################

## [B1] Replaced O(N^2) which() loop with pre-computed .I split
dt_interpolate_short_gaps_inplace <- function(dt, fps, max_gap_seconds = 1) {
  stopifnot(data.table::is.data.table(dt))

  if (!"x_interp" %in% names(dt)) dt[, x_interp := x]
  if (!"y_interp" %in% names(dt)) dt[, y_interp := y]
  if (!"notes"    %in% names(dt)) dt[, notes := NA_character_]

  data.table::setorder(dt, individual, time)

  max_gap_frames <- max(0L, floor(fps * max_gap_seconds) - 1L)
  if (max_gap_frames == 0L) return(invisible(dt))

  ## Pre-compute row indices by individual (single pass, no repeated which())
  idx_by_ind <- split(seq_len(nrow(dt)), dt$individual)

  for (ind_name in names(idx_by_ind)) {
    idx <- idx_by_ind[[ind_name]]
    n   <- length(idx)
    if (n < 3L) next

    xv <- dt$x_interp[idx]
    yv <- dt$y_interp[idx]
    tv <- dt$time[idx]

    is_na <- is.na(xv) | is.na(yv)
    imask <- internal_na_mask(is_na)
    if (!any(imask)) next

    r <- rle(imask & is_na)
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    for (j in seq_along(r$lengths)) {
      if (!r$values[j]) next
      run_len <- r$lengths[j]
      if (run_len > max_gap_frames) next
      s <- starts[j]; e <- ends[j]
      before <- s - 1L; after <- e + 1L
      if (before < 1L || after > n) next
      if (is.na(xv[before]) || is.na(xv[after])) next
      ## [audit R4] Guard against identical timestamps
      if (tv[before] == tv[after]) next

      interp_x <- approx(c(tv[before], tv[after]),
                          c(xv[before], xv[after]),
                          xout = tv[s:e])$y
      interp_y <- approx(c(tv[before], tv[after]),
                          c(yv[before], yv[after]),
                          xout = tv[s:e])$y

      global_idx <- idx[s:e]
      data.table::set(dt, i = global_idx, j = "x_interp", value = interp_x)
      data.table::set(dt, i = global_idx, j = "y_interp", value = interp_y)
      data.table::set(dt, i = global_idx, j = "notes",    value = "interpolated")
    }
  }
  invisible(dt)
}

###############################################################################
##  8b) STEP 3.5: OUTLIER REJECTION (replaces §17 jump detection)            ##
###############################################################################
# Nulls raw x/y at frames that fail either test below. MUST run BEFORE
# dt_interpolate_short_gaps_inplace() and BEFORE speed is computed, so the
# resulting gaps are filled by the existing short-gap interpolator (which
# already correctly leaves gaps longer than max_gap_seconds -- e.g. an
# unrecoverable end-of-session tracking loss -- as genuine NA rather than
# extrapolating), and so speed_cm_s is never computed from a rejected point.
#
# Returns a one-row-per-fish summary data.frame (fish_id, n_frames,
# n_rejected_oob, n_rejected_speed, threshold_px_per_s, n_passes). `dt` is
# mutated in place via data.table reference semantics, same convention as
# dt_interpolate_short_gaps_inplace() -- the caller does not need to
# reassign it, only capture the returned summary.
dt_reject_outliers_inplace <- function(dt,
                                       frame_width = NA_real_,
                                       frame_height = NA_real_,
                                       median_body_length_px = NA_real_,
                                       pad_px = 5,
                                       mad_mult = 30,
                                       floor_bl_per_s = 8,
                                       max_passes = 5L) {
  stopifnot(data.table::is.data.table(dt))
  if (!all(c("x", "y", "individual", "frame", "time") %in% names(dt))) return(data.frame())

  data.table::setorder(dt, individual, frame)
  idx_by_ind <- split(seq_len(nrow(dt)), dt$individual)
  fish_id_of <- if ("fish_id" %in% names(dt)) {
    vapply(idx_by_ind, function(idx) as.character(dt$fish_id[idx[1L]]), character(1L))
  } else {
    names(idx_by_ind)
  }

  has_frame_box <- is.finite(frame_width) && is.finite(frame_height) &&
    frame_width > 0 && frame_height > 0
  has_bl_floor  <- is.finite(median_body_length_px) && median_body_length_px > 0

  # Per-frame flag, carried through speed_by_frame -> master_fish_by_frame as
  # `jump_flag` (see STEP 12 speed_cols and the retired-§17 placeholder), so
  # a downstream consumer can still tell which frames this gate corrected --
  # same schema meaning as before, populated by the new mechanism.
  if (!"outlier_flag" %in% names(dt)) dt[, outlier_flag := FALSE]

  n_oob   <- setNames(integer(length(idx_by_ind)), names(idx_by_ind))
  n_speed <- setNames(integer(length(idx_by_ind)), names(idx_by_ind))
  thr_used <- setNames(rep(NA_real_, length(idx_by_ind)), names(idx_by_ind))
  n_passes_taken <- 0L

  # ---- (a) out-of-frame: a hard physical bound, apply once up front --------
  if (has_frame_box) {
    xv <- dt$x; yv <- dt$y
    oob <- is.finite(xv) & is.finite(yv) &
      (xv < -pad_px | xv > frame_width + pad_px | yv < -pad_px | yv > frame_height + pad_px)
    if (any(oob)) {
      data.table::set(dt, i = which(oob), j = "x", value = NA_real_)
      data.table::set(dt, i = which(oob), j = "y", value = NA_real_)
      data.table::set(dt, i = which(oob), j = "outlier_flag", value = TRUE)
      for (ind_name in names(idx_by_ind)) {
        n_oob[[ind_name]] <- sum(oob[idx_by_ind[[ind_name]]])
      }
    }
  }

  # ---- (b) implausible displacement: pooled per session, iterated ----------
  for (pass in seq_len(max_passes)) {
    per_ind_spd <- vector("list", length(idx_by_ind))
    names(per_ind_spd) <- names(idx_by_ind)
    all_spd <- vector("list", length(idx_by_ind))

    for (ind_name in names(idx_by_ind)) {
      idx <- idx_by_ind[[ind_name]]
      xv <- dt$x[idx]; yv <- dt$y[idx]; tv <- dt$time[idx]
      vld <- is.finite(xv) & is.finite(yv)
      vi  <- which(vld)
      if (length(vi) < 4L) next
      dxv <- diff(xv[vi]); dyv <- diff(yv[vi]); dtm <- diff(tv[vi])
      dtm[!is.finite(dtm) | dtm <= 0] <- NA_real_
      spd <- sqrt(dxv^2 + dyv^2) / dtm
      per_ind_spd[[ind_name]] <- list(vi = vi, spd = spd)
      all_spd[[ind_name]] <- spd[is.finite(spd)]
    }

    pooled <- unlist(all_spd, use.names = FALSE)
    if (length(pooled) < 20L) break

    med <- stats::median(pooled, na.rm = TRUE)
    mad <- stats::median(abs(pooled - med), na.rm = TRUE) * 1.4826
    thr_mad <- med + mad_mult * mad

    any_flagged <- FALSE
    for (ind_name in names(per_ind_spd)) {
      pis <- per_ind_spd[[ind_name]]
      if (is.null(pis)) next
      thr <- thr_mad
      if (has_bl_floor) thr <- max(thr, floor_bl_per_s * median_body_length_px)
      thr_used[[ind_name]] <- thr

      bad_step <- is.finite(pis$spd) & (pis$spd > thr)
      if (!any(bad_step)) next
      any_flagged <- TRUE

      # Flag the ARRIVAL point of each jump (matches the retired §17
      # convention: the point being arrived AT is the suspect one).
      arrival_local <- pis$vi[which(bad_step) + 1L]
      global_idx <- idx_by_ind[[ind_name]][arrival_local]
      data.table::set(dt, i = global_idx, j = "x", value = NA_real_)
      data.table::set(dt, i = global_idx, j = "y", value = NA_real_)
      data.table::set(dt, i = global_idx, j = "outlier_flag", value = TRUE)
      n_speed[[ind_name]] <- n_speed[[ind_name]] + length(global_idx)
    }

    n_passes_taken <- pass
    if (!any_flagged) break
  }

  data.frame(
    fish_id                = unname(fish_id_of),
    n_frames                = vapply(idx_by_ind, length, integer(1L)),
    n_rejected_oob          = unname(n_oob),
    n_rejected_speed        = unname(n_speed),
    n_rejected_total        = unname(n_oob) + unname(n_speed),
    threshold_px_per_s      = round(unname(thr_used), 2),
    n_passes                = n_passes_taken,
    stringsAsFactors        = FALSE
  )
}

###############################################################################
##  9) HEAVY CORE: data.table SPEED (IN-PLACE)                               ##
###############################################################################

## [B2] Replaced per-individual for+which() loop with single dt[, := {...}, by=individual]
dt_compute_speed_cm_s_inplace <- function(dt, fps, length_unit,
                                           prob_min = NULL, use_time_dt = TRUE) {
  stopifnot(data.table::is.data.table(dt))

  new_cols <- c("vx_cm_s", "vy_cm_s", "speed_cm_s", "heading_rad",
                "turn_angle_rad", "accel_cm_s2")
  for (cc in new_cols) dt[, (cc) := NA_real_]

  if (!is_valid_positive(fps) || !is_valid_positive(length_unit)) return(invisible(dt))

  ## Defensive coercion (warn once per column)
  coerce_warned <- character()
  for (cc in c("time", "x", "y", "x_interp", "y_interp", "frame")) {
    if (cc %in% names(dt) && !is.numeric(dt[[cc]])) {
      if (!cc %in% coerce_warned) {
        warning("Coercing column '", cc, "' to numeric for speed computation.")
        coerce_warned <- c(coerce_warned, cc)
      }
      dt[, (cc) := as.numeric(get(cc))]
    }
  }

  data.table::setorder(dt, individual, frame)

  has_prob <- !is.null(prob_min) && "probability" %in% names(dt)
  lu <- length_unit  # capture in closure

  ## [B2] Single grouped := block — no explicit for loop or which()
  dt[, c("vx_cm_s", "vy_cm_s", "speed_cm_s", "heading_rad",
         "turn_angle_rad", "accel_cm_s2") := {

    nn <- .N
    if (nn < 2L) {
      list(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
    } else {
      dt_step <- if (use_time_dt) {
        c(NA_real_, diff(time))
      } else {
        dt_s <- rep(1 / fps, nn)
        dt_s[1] <- NA_real_
        dt_s
      }

      dx_cm <- c(NA_real_, diff(x_interp)) / lu
      dy_cm <- c(NA_real_, diff(y_interp)) / lu
      vx <- dx_cm / dt_step
      vy <- dy_cm / dt_step
      spd <- sqrt(vx^2 + vy^2)
      hdg <- atan2(vy, vx)

      ta <- c(NA_real_, diff(hdg))
      ta <- ((ta + pi) %% (2 * pi)) - pi
      ta[1] <- NA_real_

      acc <- c(NA_real_, diff(spd)) / dt_step
      acc[1] <- NA_real_

      bad <- is.na(x_interp) | is.na(y_interp) | is.na(dt_step) | dt_step <= 0
      bad[1] <- TRUE

      if (has_prob) {
        pr <- probability
        bad <- bad | (!is.na(pr) & pr < prob_min)
      }

      vx[bad]  <- NA_real_; vy[bad]  <- NA_real_
      spd[bad] <- NA_real_; hdg[bad] <- NA_real_
      ta[bad]  <- NA_real_; acc[bad] <- NA_real_

      list(vx, vy, spd, hdg, ta, acc)
    }
  }, by = individual]

  invisible(dt)
}

###############################################################################
##  10) GROUP DYNAMICS (swaRmverse-driven)                                   ##
###############################################################################

dt_group_dynamics <- function(dt, fps, length_unit,
                              pol_thresh = 0.7, min_dur_s = 1.0,
                              hierarchy_rule = "A") {
  stopifnot(data.table::is.data.table(dt))

  empty_result <- list(
    group_dynamics_ind_by_frame      = data.frame(),
    group_dynamics_pairwise_by_frame = data.frame(),
    group_dynamics_by_frame          = data.frame(),
    col_mot_events_per_session       = data.frame(),
    session_event_summary            = data.frame(),
    h_stab_per_session               = data.frame()
  )

  ## [audit R6] Skip group dynamics for single-fish sessions
  n_unique_fish <- length(unique(dt$fish_in_session))
  if (n_unique_fish < 2L) {
    ts_msg("Skipping group dynamics: only ", n_unique_fish, " fish")
    return(empty_result)
  }

  dt_cm <- copy(dt)
  dt_cm[, x_cm := x_interp / length_unit]
  dt_cm[, y_cm := y_interp / length_unit]

  ## [D1] Correct expected_n: uniqueN per session
  dt_cm[, expected_n := uniqueN(fish_in_session), by = session_id]

  ## [D2] Complete-frame filter using uniqueN per (session_id, frame)
  dt_cm[, n_present := uniqueN(fish_in_session), by = .(session_id, frame)]

  ## [C1] Key-based filtering (no paste-based string matching)
  dt_cm <- dt_cm[n_present == expected_n]
  dt_cm[, c("n_present") := NULL]

  if (nrow(dt_cm) == 0L) return(empty_result)

  period_ms <- max(1L, round(1000 / fps))

  ## Convert frame numbers to POSIXct timestamps so swaRmverse bins them correctly.
  ## raw_t must be actual time values (seconds from origin) that, when discretized
  ## by 'period', group all individuals at the same frame into one time step.
  ## We use frame/fps to produce seconds, then convert to POSIXct.
  frame_seconds <- dt_cm$frame / fps
  frame_posixct <- as.POSIXct(frame_seconds, origin = "1970-01-01", tz = "UTC")

  ## Use swaRmverse workflow — wrapped in suppressMessages/suppressWarnings
  ## to prevent "Id N has only one data point" flood that freezes the console
  sdata <- tryCatch({
    suppressMessages(suppressWarnings(
      set_data_format(
        raw_x   = dt_cm$x_cm,
        raw_y   = dt_cm$y_cm,
        raw_t   = frame_posixct,
        raw_id  = as.character(dt_cm$fish_in_session),
        raw_context = dt_cm$session_id,
        period  = paste0(period_ms, "ms"),
        origin  = "1970-01-01"
      )
    ))
  }, error = function(e) {
    ts_msg("swaRmverse::set_data_format failed: ", e$message)
    NULL
  })

  if (is.null(sdata)) return(empty_result)

  ## step2time = seconds per time step (1/fps converts frames to seconds)
  step2time_val <- 1.0 / fps

  sdata <- tryCatch(
    suppressMessages(suppressWarnings(add_velocities(sdata))),
    error = function(e) { ts_msg("add_velocities: ", e$message); sdata })
  pw    <- tryCatch(
    suppressMessages(suppressWarnings(pairwise_metrics(sdata))),
    error = function(e) { ts_msg("pairwise_metrics: ", e$message); NULL })
  gm    <- tryCatch(
    suppressMessages(suppressWarnings(
      group_metrics_per_set(sdata,
                            step2time = step2time_val,
                            geo = FALSE,
                            mov_av_time_window = max(1L, round(fps)))
    )),
    error = function(e) { ts_msg("group_metrics_per_set: ", e$message); NULL })

  ## [D3] Robust list->df flattener
  safe_as_df <- function(x) {
    if (is.null(x)) return(data.frame())
    if (is.data.frame(x)) {
      df <- x
    } else if (is.list(x) && length(x) > 0 && all(vapply(x, is.data.frame, logical(1)))) {
      df <- data.table::rbindlist(x, fill = TRUE)
      df <- as.data.frame(df)
    } else if (is.list(x)) {
      df <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE),
                     error = function(e) data.frame())
    } else {
      df <- data.frame()
    }
    for (cc in names(df)) {
      ## NOTE: POSIXt columns preserved here — needed for frame recovery from t column
      if (is.factor(df[[cc]])) df[[cc]] <- as.character(df[[cc]])
    }
    df
  }

  ind_df <- safe_as_df(sdata)
  pw_df  <- safe_as_df(pw)
  gm_df  <- safe_as_df(gm)

  ## Build metadata lookup keyed by frame
  meta <- unique(dt_cm[, .(session_id, trial_id, segment_id, time, frame = as.integer(frame))])
  ## Keep one row per frame (first time value)
  meta <- meta[, .SD[1L], by = frame]

  if (nrow(gm_df) > 0 && "set" %in% names(gm_df)) {
    ## t column from swaRmverse may be POSIXct; convert back to frame number
    if ("t" %in% names(gm_df)) {
      t_vals <- gm_df$t
      if (inherits(t_vals, "POSIXt")) {
        gm_df$frame <- as.integer(round(as.numeric(t_vals) * fps))
      } else {
        gm_df$frame <- as.integer(t_vals)
      }
    }
  }

  ## [D5] Key-based join for metadata
  if (nrow(gm_df) > 0 && "frame" %in% names(gm_df)) {
    gm_dt <- as.data.table(gm_df)
    gm_dt[, frame := as.integer(frame)]
    meta_dt <- as.data.table(meta)
    setkeyv(meta_dt, "frame")
    setkeyv(gm_dt, "frame")
    gm_df <- as.data.frame(meta_dt[gm_dt, on = "frame"])
  } else if (nrow(gm_df) > 0) {
    gm_df$session_id  <- dt_cm$session_id[1]
    gm_df$trial_id    <- dt_cm$trial_id[1]
    gm_df$segment_id  <- dt_cm$segment_id[1]
  }

  ## Ensure polarization and ang_momentum columns
  if (nrow(gm_df) > 0) {
    if (!"polarization" %in% names(gm_df)) {
      pol_candidates <- grep("pol", names(gm_df), value = TRUE, ignore.case = TRUE)
      if (length(pol_candidates) > 0) names(gm_df)[names(gm_df) == pol_candidates[1]] <- "polarization"
      else gm_df$polarization <- NA_real_
    }
    if (!"ang_momentum" %in% names(gm_df)) {
      am_candidates <- grep("rotat|ang.*mom|momentum", names(gm_df), value = TRUE, ignore.case = TRUE)
      if (length(am_candidates) > 0) names(gm_df)[names(gm_df) == am_candidates[1]] <- "ang_momentum"
      else gm_df$ang_momentum <- NA_real_
    }

    gm_df$col_mot_state <- !is.na(gm_df$polarization) & gm_df$polarization > pol_thresh
    gm_df$ang_momentum <- gm_df$ang_momentum * fps
  }

  ## --- [D4] Collective motion events — DISABLED ---
  col_events <- data.frame()
  session_evt_summ <- data.frame()

  ## DISABLED: collective-motion event detection removed from pipeline outputs.
  ## Re-enable by un-commenting the block below.
  # if (nrow(gm_df) > 0 && "frame" %in% names(gm_df) && "col_mot_state" %in% names(gm_df)) {
  #   gm_dt <- as.data.table(gm_df)
  #   setorder(gm_dt, frame)
  #   gm_dt[, frame_gap := c(0L, diff(frame))]
  #   gm_dt[, block_id := cumsum(frame_gap > 1L)]
  #   gm_dt[, run_id := rleid(col_mot_state), by = block_id]
  #   evt_dt <- gm_dt[col_mot_state == TRUE,
  #                    .(start_frame = min(frame),
  #                      end_frame   = max(frame),
  #                      duration_s  = max(time, na.rm = TRUE) - min(time, na.rm = TRUE),
  #                      session_id  = first(session_id),
  #                      trial_id    = first(trial_id),
  #                      segment_id  = first(segment_id)),
  #                    by = .(block_id, run_id)]
  #   evt_dt <- evt_dt[duration_s >= min_dur_s]
  #   evt_dt[, c("block_id", "run_id") := NULL]
  #   if (nrow(evt_dt) > 0) {
  #     col_events <- as.data.frame(evt_dt)
  #     session_evt_summ <- evt_dt[, .(n_events = .N, avg_duration = mean(duration_s)),
  #                                 by = .(session_id, trial_id, segment_id)]
  #     session_evt_summ <- as.data.frame(session_evt_summ)
  #   } else {
  #     session_evt_summ <- data.frame(
  #       session_id = gm_df$session_id[1], trial_id = gm_df$trial_id[1],
  #       segment_id = gm_df$segment_id[1], n_events = 0L, avg_duration = NA_real_,
  #       stringsAsFactors = FALSE)
  #   }
  # }

  ## --- [D5] Hierarchy stability — DISABLED ---
  h_stab_df <- data.frame()

  ## DISABLED: h-stability (dominance hierarchy) removed from pipeline outputs.
  ## Re-enable by un-commenting the block below.
  # if (hierarchy_rule == "A" && nrow(ind_df) > 0 && nrow(gm_df) > 0) {
  #   tryCatch({
  #     ind_dt <- as.data.table(ind_df)
  #     if ("head" %in% names(ind_dt)) {
  #       setnames(ind_dt, "head", "ind_heading", skip_absent = TRUE)
  #     } else if ("heading" %in% names(ind_dt)) {
  #       setnames(ind_dt, "heading", "ind_heading", skip_absent = TRUE)
  #     }
  #     if ("ind_heading" %in% names(ind_dt) && "t" %in% names(ind_dt)) {
  #       if (inherits(ind_dt$t, "POSIXt")) {
  #         ind_dt[, frame := as.integer(round(as.numeric(t) * fps))]
  #       } else {
  #         ind_dt[, frame := as.integer(t)]
  #       }
  #       grp_head <- ind_dt[, .(
  #         mean_cos = mean(cos(ind_heading), na.rm = TRUE),
  #         mean_sin = mean(sin(ind_heading), na.rm = TRUE)
  #       ), by = frame]
  #       grp_head[, group_heading := atan2(mean_sin, mean_cos)]
  #       pos_dt <- dt_cm[, .(frame, fish_in_session, x_cm, y_cm,
  #                           session_id, trial_id, segment_id)]
  #       setkeyv(grp_head, "frame")
  #       setkeyv(pos_dt, "frame")
  #       pos_dt <- grp_head[pos_dt, on = "frame", nomatch = 0L]
  #       pos_dt[, fish_in_session := as.integer(fish_in_session)]
  #       pos_dt[, proj := x_cm * cos(group_heading) + y_cm * sin(group_heading)]
  #       leaders <- pos_dt[, .(leader = fish_in_session[which.max(proj)]), by = frame]
  #       h_stab_df <- data.frame(
  #         session_id = dt_cm$session_id[1],
  #         trial_id   = dt_cm$trial_id[1],
  #         segment_id = dt_cm$segment_id[1],
  #         h_stab     = max(table(leaders$leader)) / nrow(leaders),
  #         rule       = "A",
  #         stringsAsFactors = FALSE
  #       )
  #     }
  #   }, error = function(e) {
  #     ts_msg("Hierarchy stability (rule A) failed: ", e$message)
  #   })
  # } else if (hierarchy_rule == "B" && nrow(ind_df) > 0) {
  #   tryCatch({
  #     ind_dt <- as.data.table(ind_df)
  #     if ("t" %in% names(ind_dt)) {
  #       if (inherits(ind_dt$t, "POSIXt")) {
  #         ind_dt[, frame := as.integer(round(as.numeric(t) * fps))]
  #       } else {
  #         ind_dt[, frame := as.integer(t)]
  #       }
  #     }
  #     if ("speed" %in% names(ind_dt) && "id" %in% names(ind_dt) && "frame" %in% names(ind_dt)) {
  #       mid <- median(ind_dt$frame, na.rm = TRUE)
  #       first_half  <- ind_dt[frame <= mid, .(v_first = mean(speed, na.rm = TRUE)), by = id]
  #       second_half <- ind_dt[frame > mid,  .(v_second = mean(speed, na.rm = TRUE)), by = id]
  #       merged <- merge(first_half, second_half, by = "id")
  #       if (nrow(merged) >= 3L) {
  #         rho <- cor(rank(merged$v_first), rank(merged$v_second), method = "spearman")
  #       } else {
  #         rho <- NA_real_
  #       }
  #       h_stab_df <- data.frame(
  #         session_id = dt_cm$session_id[1],
  #         trial_id   = dt_cm$trial_id[1],
  #         segment_id = dt_cm$segment_id[1],
  #         h_stab     = rho,
  #         rule       = "B",
  #         stringsAsFactors = FALSE
  #       )
  #     }
  #   }, error = function(e) {
  #     ts_msg("Hierarchy stability (rule B) failed: ", e$message)
  #   })
  # }

  list(
    group_dynamics_ind_by_frame      = ind_df,
    group_dynamics_pairwise_by_frame = pw_df,
    group_dynamics_by_frame          = gm_df,
    col_mot_events_per_session       = col_events,
    session_event_summary            = session_evt_summ,
    h_stab_per_session               = h_stab_df
  )
}

###############################################################################
##  11) ZONE CLASSIFICATION + OCCUPANCY                                      ##
###############################################################################

## [C2] Accepts data.table; uses locf_fill instead of zoo::na.locf
zone_occupancy_y_boundary <- function(df, length_unit, y0_px, flow_is_higher_y = TRUE) {
  if (!is.data.table(df)) df <- as.data.table(df)
  
  # valid only if interpolated position exists
  df[, valid := !is.na(x_interp) & !is.na(y_interp)]
  
  # assign directly to flow or calm; missing positions stay NA
  df[, zone := fcase(
    !valid, NA_character_,
    flow_is_higher_y & y_interp >= y0_px, "flow",
    flow_is_higher_y & y_interp <  y0_px, "calm",
    !flow_is_higher_y & y_interp >= y0_px, "calm",
    !flow_is_higher_y & y_interp <  y0_px, "flow",
    default = NA_character_
  )]
  
  # compute occupancy seconds using valid consecutive intervals
  df[, dt_sec := time - shift(time, 1L), by = individual]
  
  occ <- df[
    !is.na(dt_sec) & dt_sec > 0 & !is.na(zone),
    .(seconds = sum(dt_sec, na.rm = TRUE)),
    by = .(individual, zone)
  ]
  
  totals <- occ[, .(total_valid_seconds = sum(seconds, na.rm = TRUE)), by = individual]
  occ <- totals[occ, on = "individual"]
  occ[, prop_time := seconds / total_valid_seconds]
  
  list(
    df_z = as.data.frame(df),
    occ_seconds = as.data.frame(occ),
    y0_px = y0_px,
    flow_is_higher_y = flow_is_higher_y
  )
}
######################################################################################################################################

###############################################################################
## PRE-FLIGHT (ROBUST): check required files per session
## Strategy:
##   1) FAST PATH (“normal search route”): check canonical locations.
##   2) FALLBACK (only if fast path fails for a file):
##        - recursively search within the session folder
##        - tolerate minor folder-structure differences
##        - tolerate minor naming differences via patterns
###############################################################################

## ---- Canonical relative paths (fast path) ----
CANON_TRAJ_REL <- "trajectories/trajectories_csv/trajectories.csv"
CANON_ATTR_REL <- "trajectories/trajectories_csv/attributes.json"
CANON_AREAS_REL <- "trajectories/trajectories_csv/areas.csv"  # optional

## ---- Fallback patterns (used ONLY if canonical is missing) ----
## You can tighten/loosen these if you know your variants.
PAT_TRAJ <- "(?i)^trajectories(.*)?\\.csv$"
PAT_ATTR <- "(?i)^attributes(.*)?\\.json$"
PAT_AREAS <- "(?i)^areas(.*)?\\.csv$"
PAT_BBOX <- "(?i)_bboxes\\.csv$"   # for info only; not required for preflight

## Choose a “best” hit among multiple recursive matches
.choose_best_hit <- function(paths, prefer_dir_regex = "(?i)trajectories_csv") {
  if (length(paths) == 0L) return(NA_character_)
  paths <- norm_path(paths, mustWork = FALSE)
  
  ## Prefer hits under trajectories_csv, then shorter paths (closer to session root)
  score <- rep(0L, length(paths))
  score <- score + ifelse(grepl(prefer_dir_regex, paths, perl = TRUE), 100L, 0L)
  score <- score - nchar(paths)  # shorter is better
  paths[which.max(score)]
}

## Recursive search helper (ONLY invoked when canonical missing)
.find_in_session_recursive <- function(session_dir, filename_pattern, max_hits = 200L) {
  if (!dir.exists(session_dir)) return(character())
  hits <- list.files(
    session_dir,
    pattern = filename_pattern,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  if (length(hits) > max_hits) hits <- hits[seq_len(max_hits)]
  hits[file.exists(hits)]
}

## Core resolver for one “required file”
.resolve_required_file <- function(session_dir, canonical_rel, fallback_pattern) {
  canonical_full <- file.path(session_dir, canonical_rel)
  if (file.exists(canonical_full)) {
    return(list(path = norm_path(canonical_full), method = "canonical"))
  }
  
  ## Only now do we do the expensive recursive scan
  hits <- .find_in_session_recursive(session_dir, fallback_pattern)
  if (length(hits) == 0L) {
    return(list(path = NA_character_, method = "missing"))
  }
  
  best <- .choose_best_hit(hits)
  list(path = best, method = "fallback_recursive")
}

## Check one session folder
check_one_session_files_robust <- function(session_dir,
                                           require_areas = FALSE,
                                           include_bbox_info = TRUE) {
  sid <- basename(session_dir)
  
  traj  <- .resolve_required_file(session_dir, CANON_TRAJ_REL, PAT_TRAJ)
  attr  <- .resolve_required_file(session_dir, CANON_ATTR_REL, PAT_ATTR)
  sjson <- .resolve_required_file(session_dir, "session.json", "(?i)^session\\.json$")
  
  areas <- list(path = NA_character_, method = NA_character_)
  if (isTRUE(require_areas)) {
    areas <- .resolve_required_file(session_dir, CANON_AREAS_REL, PAT_AREAS)
  }
  
  ## Bboxes are not required for preflight, only reported if present
  bbox_path <- NA_character_
  bbox_method <- NA_character_
  if (isTRUE(include_bbox_info)) {
    if (exists("find_bboxes_in_session", mode = "function")) {
      bb <- find_bboxes_in_session(session_dir)
      if (!is.na(bb) && file.exists(bb)) {
        bbox_path <- norm_path(bb)
        bbox_method <- "find_bboxes_in_session"
      } else {
        bb_hits <- .find_in_session_recursive(session_dir, PAT_BBOX)
        bbox_path <- .choose_best_hit(bb_hits)
        bbox_method <- ifelse(!is.na(bbox_path), "fallback_recursive", "missing")
      }
    } else {
      bb_hits <- .find_in_session_recursive(session_dir, PAT_BBOX)
      bbox_path <- .choose_best_hit(bb_hits)
      bbox_method <- ifelse(!is.na(bbox_path), "fallback_recursive", "missing")
    }
  }
  
  ## Required logic:
  ## - trajectories.csv is always required
  ## - metadata is valid if EITHER attributes.json OR session.json exists
  missing <- c()
  if (is.na(traj$path)) missing <- c(missing, "trajectories.csv")
  if (is.na(attr$path) && is.na(sjson$path)) {
    missing <- c(missing, "attributes.json OR session.json")
  }
  if (isTRUE(require_areas) && is.na(areas$path)) {
    missing <- c(missing, "areas.csv")
  }
  
  tibble::tibble(
    session_id  = sid,
    session_dir = norm_path(session_dir),
    
    trajectories_path = traj$path,
    trajectories_found_by = traj$method,
    
    attributes_path = attr$path,
    attributes_found_by = attr$method,
    
    session_json_path = sjson$path,
    session_json_found_by = sjson$method,
    
    areas_path = if (isTRUE(require_areas)) areas$path else NA_character_,
    areas_found_by = if (isTRUE(require_areas)) areas$method else NA_character_,
    
    bboxes_path = if (isTRUE(include_bbox_info)) bbox_path else NA_character_,
    bboxes_found_by = if (isTRUE(include_bbox_info)) bbox_method else NA_character_,
    
    ok_required = length(missing) == 0L,
    missing_required = paste(missing, collapse = " | ")
  )
}

## Scan all sessions under base_dir
preflight_check_sessions_robust <- function(base_dir,
                                            session_dir_regex = "(?i)^session_",
                                            require_areas = FALSE,
                                            include_bbox_info = TRUE,
                                            print_ok_sessions = FALSE) {
  base_dir <- norm_path(base_dir, mustWork = TRUE)
  
  sess_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  sess_dirs <- sess_dirs[grepl(session_dir_regex, basename(sess_dirs), perl = TRUE)]
  
  if (length(sess_dirs) == 0L) {
    ts_msg("No session folders matched regex: ", session_dir_regex)
    return(invisible(list(report_all = tibble::tibble(), report_missing = tibble::tibble())))
  }
  
  rep <- dplyr::bind_rows(lapply(
    sess_dirs,
    check_one_session_files_robust,
    require_areas = require_areas,
    include_bbox_info = include_bbox_info
  )) %>%
    dplyr::arrange(session_id)
  
  missing <- rep %>% dplyr::filter(!ok_required)
  
  ts_msg("Pre-flight (robust): scanned ", nrow(rep), " sessions. ",
         nrow(missing), " sessions are missing REQUIRED files.")
  
  if (nrow(missing) > 0) {
    print(
      missing %>% dplyr::select(session_id, missing_required, trajectories_found_by, attributes_found_by),
      n = Inf
    )
  } else if (isTRUE(print_ok_sessions)) {
    print(rep %>% dplyr::select(session_id, ok_required), n = Inf)
  }
  
  invisible(list(report_all = rep, report_missing = missing))
}

###############################################################################
## USE IT (right before run_batch_folder)
###############################################################################

base_dir <- clipboard_dir()

pre <- preflight_check_sessions_robust(
  base_dir = base_dir,
  require_areas = FALSE,        # set TRUE if you want areas.csv required
  include_bbox_info = TRUE      # informative only; not required
)

## Stop if required files are missing
if (nrow(pre$report_missing) > 0) {
  stop("Pre-flight failed: some sessions are missing required files. See report above.")
}

######################################################################################################################################

###############################################################################
##  13) PER-SESSION PROCESSING                                               ##
###############################################################################

process_one_session <- function(session_id, files_tbl,
                                length_unit_fallback,
                                conda_exe, conda_env, script_path,
                                max_gap_seconds, speed_prob_min, use_time_dt,
                                prob_threshold, y0_px, flow_is_higher_y,
                                ensure_bboxes = TRUE, log_dir = NULL) {
  
  ## Helper: report step completion with row count
  nrow_safe <- function(x) {
    if (is.null(x)) return(0L)
    if (is.data.frame(x) || data.table::is.data.table(x)) return(nrow(x))
    0L
  }
  step_done <- function(step_num, label, result = NULL) {
    row_info <- if (!is.null(result)) paste0(" [", nrow_safe(result), " rows]") else ""
    ts_msg("  STEP ", sprintf("%02d", step_num), " completed: ", label, row_info)
  }
  
  ts_msg("========== Processing: ", session_id, " ==========")
  
  ## ---- STEP 01: Parse trial & segment IDs from session name ----
  ts <- parse_trial_segment(session_id)
  trial_id   <- ts$trial_id
  segment_id <- ts$segment_id
  step_done(1, "Parse trial/segment")
  
  ## ---- STEP 02: Read wide trajectories + attributes, pivot to long ----
  sess <- read_session_wide(session_id, files_tbl,
                            length_unit_fallback = length_unit_fallback)
  fps         <- sess$fps
  length_unit <- sess$length_unit
  n_fish      <- sess$n_fish

  # Frame geometry + median body length for STEP 3.5 (NA if this session's
  # attributes.json/session.json didn't carry them -- the gate degrades to
  # the displacement rule alone in that case).
  frame_width            <- sess$frame_width
  frame_height           <- sess$frame_height
  median_body_length_px  <- sess$median_body_length_px

  dat <- sess$trajectories_long
  dat[, session_id      := session_id]
  dat[, trial_id        := trial_id]
  dat[, segment_id      := segment_id]
  dat[, fish_in_session := individual]
  dat[, fish_id := paste0("T", trial_id, "_S", segment_id, "_F", fish_in_session)]
  step_done(2, sprintf("Read trajectories (fps=%.1f, lu=%.2f, n_fish=%d)", fps, length_unit, n_fish), dat)
  
  ## ---- STEP 03: Compute gap statistics BEFORE interpolation ----
  compute_gaps_dt <- function(dt_sub, stage_label) {
    dt_sub <- copy(dt_sub)
    if (stage_label == "after") {
      dt_sub[, gx := x_interp]
      dt_sub[, gy := y_interp]
    } else {
      dt_sub[, gx := x]
      dt_sub[, gy := y]
    }
    dt_sub[, is_na := is.na(gx) | is.na(gy)]
    
    if (!any(dt_sub$is_na)) return(data.table())
    
    setorder(dt_sub, fish_id, frame)
    dt_sub[, run_id := rleid(is_na), by = fish_id]
    
    na_runs <- dt_sub[is_na == TRUE, .(
      start_frame  = min(frame),
      end_frame    = max(frame),
      start_time   = min(time),
      end_time     = max(time),
      gap_frames   = .N,
      first_row_in_grp = min(.I),
      last_row_in_grp  = max(.I)
    ), by = .(fish_id, run_id)]
    
    if (nrow(na_runs) == 0L) return(data.table())
    
    fish_bounds <- dt_sub[, .(min_I = min(.I), max_I = max(.I)), by = fish_id]
    na_runs <- fish_bounds[na_runs, on = "fish_id"]
    na_runs[, internal_gap := (first_row_in_grp > min_I) & (last_row_in_grp < max_I)]
    na_runs[, gap_seconds := end_time - start_time]
    
    na_runs[, `:=`(session_id = session_id, trial_id = trial_id,
                   segment_id = segment_id, stage = stage_label)]
    
    na_runs[, .(session_id, trial_id, segment_id, fish_id,
                start_frame, end_frame, start_time, end_time,
                gap_frames, gap_seconds, internal_gap, stage)]
  }
  
  gaps_before <- compute_gaps_dt(dat, "before")
  step_done(3, "Gaps BEFORE interpolation", gaps_before)

  ## ---- STEP 03.5: Reject outlier trajectory points (BEFORE interpolation,
  ##      BEFORE speed is computed -- see JUMP_DETECTION_ENABLED comment) ----
  outlier_rejection_by_fish <- data.frame()
  if (isTRUE(if (exists("OUTLIER_REJECT_ENABLED", inherits = FALSE)) OUTLIER_REJECT_ENABLED else TRUE)) {
    .or_pad      <- if (exists("OUTLIER_FRAME_PAD_PX", inherits = FALSE)) OUTLIER_FRAME_PAD_PX else 5
    .or_madmult  <- if (exists("OUTLIER_MAD_MULT",     inherits = FALSE)) OUTLIER_MAD_MULT     else 30
    .or_floorbls <- if (exists("OUTLIER_FLOOR_BL_S",   inherits = FALSE)) OUTLIER_FLOOR_BL_S   else 8
    .or_passes   <- if (exists("OUTLIER_MAX_PASSES",   inherits = FALSE)) OUTLIER_MAX_PASSES   else 5L

    outlier_rejection_by_fish <- dt_reject_outliers_inplace(
      dat,
      frame_width            = frame_width,
      frame_height           = frame_height,
      median_body_length_px  = median_body_length_px,
      pad_px                 = .or_pad,
      mad_mult               = .or_madmult,
      floor_bl_per_s         = .or_floorbls,
      max_passes             = .or_passes
    )
    if (nrow(outlier_rejection_by_fish) > 0) {
      outlier_rejection_by_fish$session_id <- session_id
      outlier_rejection_by_fish$trial_id   <- trial_id
      outlier_rejection_by_fish$segment_id <- segment_id
    }
  }
  .n_rejected <- if (nrow(outlier_rejection_by_fish) > 0) sum(outlier_rejection_by_fish$n_rejected_total) else 0L
  step_done(35, sprintf("Outlier rejection (%d frames nulled)", .n_rejected), outlier_rejection_by_fish)

  ## ---- STEP 04: Interpolate short internal gaps in-place ----
  ## Fills the gaps STEP 03.5 just created (plus any that were already NA
  ## from tracking), same as before -- runs are only filled when they are
  ## <= max_gap_seconds AND bracketed by valid data on both sides, so the
  ## unrecoverable end-of-session tail an outlier run can leave behind stays
  ## genuinely missing rather than being extrapolated.
  dt_interpolate_short_gaps_inplace(dat, fps, max_gap_seconds)
  n_interp <- sum(dat$notes == "interpolated", na.rm = TRUE)
  step_done(4, sprintf("Interpolation (%d frames filled)", n_interp))
  
  ## ---- STEP 05: Compute gap statistics AFTER interpolation ----
  gaps_after <- compute_gaps_dt(dat, "after")
  gaps_all <- rbindlist(list(gaps_before, gaps_after), fill = TRUE)
  step_done(5, "Gaps AFTER interpolation", gaps_after)
  
  ## ---- STEP 06: Locate or generate bounding-box CSV ----
  session_dir <- files_tbl$session_dir[files_tbl$session_id == session_id]
  bb_path <- find_bboxes_in_session(session_dir)
  
  if (is.na(bb_path) && isTRUE(ensure_bboxes)) {
    bb_path <- ensure_bboxes_for_session(session_dir, conda_exe, conda_env, script_path,
                                         log_dir = log_dir)
  }
  step_done(6, paste0("Bboxes lookup (", ifelse(is.na(bb_path), "NOT FOUND", "found"), ")"))
  
  ## ---- STEP 07: Sanity-check bbox IDs vs trajectory fish count ----
  traj_path <- sess$trajectories_path
  sanity <- sanity_check_bboxes_vs_trajectories(
    bboxes_path = bb_path,
    trajectories_path = traj_path,
    allow_drop_zero = ALLOW_DROP_SPURIOUS_BBOX_ID_ZERO
  )
  
  bboxes_mismatch        <- sanity$mismatch
  bboxes_mismatch_reason <- sanity$reason
  bboxes_ids_raw         <- if (length(sanity$ids) > 0) paste(sanity$ids, collapse = ",") else NA_character_
  bboxes_ids_used        <- if (length(sanity$ids_used) > 0) paste(sanity$ids_used, collapse = ",") else NA_character_
  bboxes_drop_ids        <- if (length(sanity$drop_ids) > 0) paste(sanity$drop_ids, collapse = ",") else NA_character_
  bboxes_autocorrected   <- isTRUE(sanity$autocorrected)
  
  step_done(7, paste0("Sanity check (", bboxes_mismatch_reason, ")"))
  
  ## ---- STEP 08: Body length per fish (from bboxes) ----
  body_by_fish <- NULL
  if (isFALSE(bboxes_mismatch)) {
    if (bboxes_autocorrected) {
      ts_msg("Autocorrecting bbox identities for ", session_id,
             " by dropping ID(s): ", bboxes_drop_ids,
             " | raw IDs: ", bboxes_ids_raw,
             " | used IDs: ", bboxes_ids_used)
    }
    
    bboxes_df <- read_bboxes_safe(
      bb_path,
      drop_ids = sanity$drop_ids
    )
    
    if (!is.null(bboxes_df)) {
      body_by_fish <- compute_bodylength_by_fish(
        as.data.frame(dat[, .(
          session_id, trial_id, segment_id, fish_id,
          fish_in_session, frame, individual, x_interp, y_interp
        )]),
        bboxes_df,
        length_unit
      )
    }
  }
  step_done(8, "Body length", body_by_fish)
  
  ## ---- STEP 09: Zone occupancy ----
  zone_res <- zone_occupancy_y_boundary(
    dat, length_unit, y0_px, flow_is_higher_y
  )
  
  zone_by_frame <- zone_res$df_z[, intersect(
    names(zone_res$df_z),
    c("session_id", "trial_id", "segment_id", "fish_id", "fish_in_session",
      "frame", "time", "valid", "zone", "y0_px", "flow_is_higher_y")
  )]
  
  occ_seconds <- zone_res$occ_seconds
  occ_seconds$session_id  <- session_id
  occ_seconds$trial_id    <- trial_id
  occ_seconds$segment_id  <- segment_id
  occ_seconds$fish_id     <- paste0("T", trial_id, "_S", segment_id, "_F", occ_seconds$individual)
  step_done(9, "Zone occupancy", occ_seconds)
  
  ## ---- STEP 10: QC summary per fish ----
  qc <- dat[, .(
    fps                    = fps,
    length_unit_px_per_cm  = length_unit,
    length_unit_fallback_used = sess$length_unit_fallback_used,
    fraction_identified    = sum(!is.na(x) & !is.na(y)) / .N,
    interpolated_frames    = sum(notes == "interpolated", na.rm = TRUE),
    bboxes_identity_mismatch        = ifelse(is.na(bboxes_mismatch), NA, bboxes_mismatch),
    bboxes_identity_mismatch_reason = bboxes_mismatch_reason,
    bboxes_ids_raw                  = bboxes_ids_raw,
    bboxes_ids_used                 = bboxes_ids_used,
    bboxes_drop_ids                 = bboxes_drop_ids,
    bboxes_autocorrected            = bboxes_autocorrected
  ), by = .(fish_id, session_id, trial_id, segment_id, fish_in_session)]
  
  qc <- as.data.frame(qc)
  
  if (!is.null(body_by_fish) && nrow(body_by_fish) > 0) {
    morph_cols <- body_by_fish[, c("fish_id", "n_bbox_frames", "body_length_cm",
                                   "body_length_sd_cm", "body_length_unstable",
                                   "body_length_unstable_reason")]
    qc <- left_join(qc, morph_cols, by = "fish_id")
  }
  step_done(10, "QC by fish", qc)
  
  ## ---- STEP 11: Position by frame ----
  position_by_frame <- dat[, .(session_id, trial_id, segment_id, fish_id,
                               fish_in_session, frame, time, x, y,
                               x_interp, y_interp, notes)]
  step_done(11, "Position by frame", position_by_frame)
  
  ## ---- STEP 12: Speed, heading, acceleration ----
  dt_compute_speed_cm_s_inplace(dat, fps, length_unit,
                                prob_min = speed_prob_min, use_time_dt = use_time_dt)
  
  ## keep all computed kinematics
  speed_cols <- intersect(
    names(dat),
    c(
      "individual",
      "frame",
      "time",
      "x",
      "y",
      "probability",
      "session_id",
      "trial_id",
      "segment_id",
      "fish_in_session",
      "fish_id",
      "x_interp",
      "y_interp",
      "notes",
      "vx_cm_s",
      "vy_cm_s",
      "speed_cm_s",
      "heading_rad",
      "turn_angle_rad",
      "accel_cm_s2",
      "outlier_flag"
    )
  )
  
  speed_by_frame <- dat[, .SD, .SDcols = speed_cols]
  speed_by_frame[, fps := fps]
  speed_by_frame[, length_unit_px_per_cm := length_unit]
  
  speed_by_fish <- dat[, .(
    mean_speed_cm_s   = mean(speed_cm_s, na.rm = TRUE),
    median_speed_cm_s = median(speed_cm_s, na.rm = TRUE),
    p95_speed_cm_s    = as.numeric(quantile(speed_cm_s, 0.95, na.rm = TRUE)),
    n_speed           = sum(!is.na(speed_cm_s))
  ), by = .(fish_id, session_id, trial_id, segment_id, fish_in_session)]
  step_done(12, "Speed computation", speed_by_fish)
  
  ## ---- STEP 13: Group dynamics ----
  gd <- dt_group_dynamics(dat, fps, length_unit,
                          pol_thresh = 0.7, min_dur_s = 1.0,
                          hierarchy_rule = "A")
  
  group_by_session <- data.frame()
  if (nrow(gd$group_dynamics_by_frame) > 0) {
    gdf <- as.data.table(gd$group_dynamics_by_frame)
    group_by_session <- gdf[, .(
      n_frames                = .N,
      mean_polarization       = mean(polarization, na.rm = TRUE),
      mean_ang_momentum_cm2_s = mean(ang_momentum, na.rm = TRUE)
    ), by = .(session_id, trial_id, segment_id)]
    group_by_session <- as.data.frame(group_by_session)
  }
  step_done(13, "Group dynamics", gd$group_dynamics_by_frame)
  
  ts_msg("  All 13 steps completed for: ", session_id)
  
  list(
    speed_by_fish                    = as.data.frame(speed_by_fish),
    group_by_session                 = group_by_session,
    occupancy_seconds                = as.data.frame(occ_seconds),
    qc_by_fish                       = qc,
    gaps                             = as.data.frame(gaps_all),
    outlier_rejection_by_fish        = outlier_rejection_by_fish,
    body_by_fish                     = body_by_fish,
    position_by_frame                = as.data.frame(position_by_frame),
    speed_by_frame                   = as.data.frame(speed_by_frame),
    group_dynamics_ind_by_frame      = gd$group_dynamics_ind_by_frame,
    group_dynamics_pairwise_by_frame = gd$group_dynamics_pairwise_by_frame,
    group_dynamics_by_frame          = gd$group_dynamics_by_frame,
    col_mot_events_per_session       = gd$col_mot_events_per_session,
    session_event_summary            = gd$session_event_summary,
    h_stab_per_session               = gd$h_stab_per_session,
    zone_by_frame                    = zone_by_frame
  )
}


###############################################################################
##  14) BATCH RUNNER                                                         ##
###############################################################################

## [A3] Added ensure_bboxes toggle (default TRUE)
run_batch_folder <- function(base_dir,
                             out_dir = file.path(base_dir, "master_outputs"),
                             n_fish_max = DEFAULT_N_FISH_MAX,
                             max_gap_seconds = DEFAULT_MAX_GAP_SECONDS,
                             speed_prob_min = DEFAULT_SPEED_PROB_MIN,
                             use_time_dt = DEFAULT_USE_TIME_DT,
                             prob_threshold = DEFAULT_PROB_THRESHOLD,
                             y0_px = DEFAULT_Y0_PX,
                             flow_is_higher_y = DEFAULT_FLOW_HIGHER_Y,
                             conda_exe,
                             conda_env,
                             script_path,
                             ensure_bboxes = TRUE) {
  
  base_dir  <- norm_path(base_dir, mustWork = TRUE)
  
  ## ---- Wipe out_dir completely if it already exists -------------------------
  if (dir.exists(out_dir)) {
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
  ensure_dir(out_dir)
  out_dir     <- norm_path(out_dir, mustWork = TRUE)
  
  conda_exe   <- norm_path(conda_exe, mustWork = TRUE)
  script_path <- norm_path(script_path, mustWork = TRUE)
  
  ## ---- Create only requested directory structure ----------------------------
  ## R7: removed "events" and "group_plus" (collective-motion events and
  ## h-stability are DISABLED upstream — see ~lines 1374, 1409 — so those
  ## folders were always created but never written to). Re-add them here if
  ## the upstream blocks are ever re-enabled.
  subdirs <- c(
    "logs",
    "morph",
    "qc",
    "by_frame/position",
    "by_frame/speed",
    "by_frame/occupancy",
    "by_frame/group_ind",
    "by_frame/group_pairwise",
    "by_frame/group_group"
  )
  for (sd in subdirs) ensure_dir(file.path(out_dir, sd))
  
  log_dir <- file.path(out_dir, "logs")
  
  ## Index sessions
  files_index <- index_session_files_smart(base_dir)
  ts_msg("Found ", nrow(files_index), " sessions in ", base_dir)

  ## OPT-0: Apply scope filters (KEEP_TP1, DROP_PHARMA) from GlobalEnv --------
  ## Extract trial_id and segment_id from session_id strings so we can filter
  ## before entering the per-session loop.  This is the EARLIEST point at which
  ## skipped sessions are identified; the log lines below are what the user
  ## should check to confirm the filters are active.
  .m_opt0          <- stringr::str_match(files_index$session_id,
                       "(?i)trial[_[:space:]]*(\\d+).*segment[_[:space:]]*(\\d+)")
  .trial_id_opt0   <- suppressWarnings(as.integer(.m_opt0[, 2]))
  .segment_id_opt0 <- suppressWarnings(as.integer(.m_opt0[, 3]))
  .n_before_opt0   <- nrow(files_index)
  .opt0_changed    <- FALSE

  ## Drop Segment 1 (Timepoint 1) if KEEP_TP1 = FALSE
  if (exists("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE) &&
      !isTRUE(get("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE))) {
    .n_tp1  <- sum(.segment_id_opt0 == 1L, na.rm = TRUE)
    .keep   <- is.na(.segment_id_opt0) | .segment_id_opt0 != 1L
    files_index      <- files_index[.keep, , drop = FALSE]
    .trial_id_opt0   <- .trial_id_opt0[.keep]
    .segment_id_opt0 <- .segment_id_opt0[.keep]
    ts_msg("OPT-0: dropped ", .n_tp1, " Segment1 sessions (KEEP_TP1=FALSE); ",
           nrow(files_index), " sessions remain")
    .opt0_changed <- TRUE
  }

  ## Drop pharma-group trials (Morphine / Naloxone) if DROP_PHARMA = TRUE
  if (exists("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE) &&
      isTRUE(get("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE)) &&
      exists("PHARMA_GROUPS", envir = .GlobalEnv, inherits = FALSE)) {
    .pg   <- get("PHARMA_GROUPS", envir = .GlobalEnv, inherits = FALSE)
    .xlsx <- if (exists("TRIAL_SUMMARY_XLSX", envir = .GlobalEnv, inherits = FALSE))
               get("TRIAL_SUMMARY_XLSX", envir = .GlobalEnv, inherits = FALSE)
             else NA_character_
    if (!is.na(.xlsx) && file.exists(.xlsx)) {
      .tmeta <- tryCatch(
        readxl::read_excel(.xlsx) %>%
          dplyr::filter(grepl(paste0(.pg, collapse = "|"),
                              .data[["Group"]], ignore.case = TRUE)),
        error = function(e) {
          ts_msg("OPT-0 WARNING: could not read TRIAL_SUMMARY_XLSX for pharma filter: ",
                 conditionMessage(e))
          NULL
        }
      )
      if (!is.null(.tmeta) && nrow(.tmeta) > 0 && "Trial" %in% names(.tmeta)) {
        .pharma_ids <- unique(suppressWarnings(as.integer(.tmeta[["Trial"]])))
        .pharma_ids <- .pharma_ids[!is.na(.pharma_ids)]
        .n_pharma   <- sum(.trial_id_opt0 %in% .pharma_ids, na.rm = TRUE)
        .keep       <- is.na(.trial_id_opt0) | !.trial_id_opt0 %in% .pharma_ids
        files_index      <- files_index[.keep, , drop = FALSE]
        .trial_id_opt0   <- .trial_id_opt0[.keep]
        .segment_id_opt0 <- .segment_id_opt0[.keep]
        ts_msg("OPT-0: dropped ", .n_pharma, " pharma-group sessions",
               " (DROP_PHARMA=TRUE; groups: ", paste(.pg, collapse = ", "), "); ",
               nrow(files_index), " sessions remain")
        .opt0_changed <- TRUE
      } else {
        ts_msg("OPT-0: DROP_PHARMA=TRUE but no matching pharma rows found in metadata",
               " — check PHARMA_GROUPS and TRIAL_SUMMARY_XLSX")
      }
    } else {
      ts_msg("OPT-0 WARNING: DROP_PHARMA=TRUE but TRIAL_SUMMARY_XLSX not found",
             " — pharma filter skipped (path: ", .xlsx, ")")
    }
  }

  if (.opt0_changed) {
    ts_msg("OPT-0 summary: ", .n_before_opt0, " sessions found -> ",
           nrow(files_index), " sessions will be processed")
  }
  rm(.m_opt0, .trial_id_opt0, .segment_id_opt0, .n_before_opt0, .opt0_changed)
  ## --------------------------------------------------------------------------

  ## Compute fallback length_unit across sessions
  lu_vals <- numeric(0)
  for (i in seq_len(nrow(files_index))) {
    ap <- files_index$attributes_path[i]
    if (!is.na(ap) && file.exists(ap)) {
      att <- tryCatch(read_session_attributes(ap), error = function(e) NULL)
      if (!is.null(att) && is_valid_positive(att$length_unit)) {
        lu_vals <- c(lu_vals, att$length_unit)
      }
    }
  }
  length_unit_fallback <- if (length(lu_vals) > 0) mean(lu_vals) else NA_real_
  ts_msg("Length unit fallback = ", round(length_unit_fallback, 4))
  
  ## Iterate
  results <- list()
  log_rows <- list()
  
  n_sessions <- nrow(files_index)
  ts_msg("Starting batch: ", n_sessions, " sessions")
  pb <- txtProgressBar(min = 0, max = n_sessions, style = 3, char = "=")
  
  for (i in seq_len(n_sessions)) {
    sid <- files_index$session_id[i]
    t0  <- Sys.time()
    
    has_traj  <- !is.na(files_index$trajectories_path[i])
    has_attr  <- !is.na(files_index$attributes_path[i])
    has_sjson <- "session_json_path" %in% names(files_index) && !is.na(files_index$session_json_path[i])
    
    if (!has_traj || !(has_attr || has_sjson)) {
      ts_msg("SKIP ", sid, ": missing trajectories or both attributes.json and session.json")
      log_rows[[i]] <- data.frame(
        session_id = sid,
        status = "skipped",
        error = "missing_files",
        runtime_s = 0,
        stringsAsFactors = FALSE
      )
      setTxtProgressBar(pb, i)
      next
    }
    
    res <- tryCatch({
      out <- process_one_session(
        session_id           = sid,
        files_tbl            = files_index,
        length_unit_fallback = length_unit_fallback,
        conda_exe            = conda_exe,
        conda_env            = conda_env,
        script_path          = script_path,
        max_gap_seconds      = max_gap_seconds,
        speed_prob_min       = speed_prob_min,
        use_time_dt          = use_time_dt,
        prob_threshold       = prob_threshold,
        y0_px                = y0_px,
        flow_is_higher_y     = flow_is_higher_y,
        ensure_bboxes        = ensure_bboxes,
        log_dir              = log_dir
      )
      list(status = "ok", data = out, error = NA_character_)
    }, error = function(e) {
      err_msg <- e$message
      err_log <- file.path(log_dir, paste0("error_", sid, ".log"))
      tryCatch(writeLines(c(as.character(Sys.time()), err_msg), err_log), error = function(e2) NULL)
      ts_msg("ERROR in ", sid, " (see ", err_log, ")")
      list(status = "error", data = NULL, error = err_msg)
    })
    
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ts_msg(sid, " -> ", res$status, " (", round(elapsed, 1), "s)")
    
    log_rows[[i]] <- data.frame(
      session_id = sid,
      status = res$status,
      error = res$error,
      runtime_s = elapsed,
      stringsAsFactors = FALSE
    )
    
    if (res$status == "ok") results[[sid]] <- res$data
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  batch_log <- bind_rows(log_rows)
  
  ## Bind master tables
  safe_bind <- function(field) {
    dfs <- purrr::compact(purrr::map(results, field))
    dfs <- purrr::keep(dfs, ~ is.data.frame(.) && nrow(.) > 0)
    if (length(dfs) == 0) return(data.frame())
    bind_rows(dfs)
  }
  
  master_qc          <- safe_bind("qc_by_fish")
  master_outlier_rej <- safe_bind("outlier_rejection_by_fish")
  master_body_by_fish <- safe_bind("body_by_fish")
  master_position    <- safe_bind("position_by_frame")
  master_speed_frame <- safe_bind("speed_by_frame")
  master_zone_frame  <- safe_bind("zone_by_frame")
  master_gd_ind      <- safe_bind("group_dynamics_ind_by_frame")
  master_gd_pw       <- safe_bind("group_dynamics_pairwise_by_frame")
  master_gd_grp      <- safe_bind("group_dynamics_by_frame")
  master_col_events  <- safe_bind("col_mot_events_per_session")
  master_evt_summary <- safe_bind("session_event_summary")
  master_h_stab      <- safe_bind("h_stab_per_session")
  
  ## Write CSVs
  write_if <- function(df, path) {
    if (is.data.frame(df) && nrow(df) > 0) {
      readr::write_csv(df, path)
      ts_msg("Wrote: ", basename(path), " [", nrow(df), " rows x ", ncol(df), " cols]")
    }
  }
  
  ## ---- Only requested top-level outputs ------------------------------------
  write_if(master_qc,         file.path(out_dir, "qc", "master_qc_by_fish.csv"))
  write_if(master_outlier_rej, file.path(out_dir, "qc", "outlier_rejection_summary.csv"))
  readr::write_csv(batch_log, file.path(out_dir, "logs", "batch_log.csv"))
  ts_msg("Wrote: batch_log.csv [", nrow(batch_log), " rows]")
  
  ## ---- Morph ----------------------------------------------------------------
  if (is.data.frame(master_body_by_fish) && nrow(master_body_by_fish) > 0) {
    write_if(master_body_by_fish, file.path(out_dir, "morph", "master_bodylength_by_fish.csv"))
    bl_wide <- make_body_length_summary_wide(master_body_by_fish, n_fish_max)
    write_if(bl_wide, file.path(out_dir, "morph", "master_bodylength_summary_by_session.csv"))
  }
  
  ## ---- By-frame -------------------------------------------------------------
  write_if(master_position,    file.path(out_dir, "by_frame", "position", "master_position_by_frame.csv"))
  write_if(master_speed_frame, file.path(out_dir, "by_frame", "speed", "master_speed_by_frame.csv"))
  write_if(master_zone_frame,  file.path(out_dir, "by_frame", "occupancy", "master_zone_by_frame.csv"))
  
  ## ---- Group dynamics by-frame ---------------------------------------------
  write_if(master_gd_ind, file.path(out_dir, "by_frame", "group_ind", "master_group_dynamics_ind_by_frame.csv"))
  write_if(master_gd_pw,  file.path(out_dir, "by_frame", "group_pairwise", "master_group_dynamics_pairwise_by_frame.csv"))
  write_if(master_gd_grp, file.path(out_dir, "by_frame", "group_group", "master_group_dynamics_by_frame.csv"))
  
  ## ---- Group plus / Events --------------------------------------------------
  ## R7: write_if() calls removed — the producers (master_h_stab,
  ## master_evt_summary, master_col_events) are DISABLED upstream and always
  ## return empty data frames, so write_if was a no-op every run.
  ## Re-add the write_if calls + the matching subdir creation above if those
  ## upstream blocks are ever re-enabled.
  
  ts_msg("Batch complete. ", sum(batch_log$status == "ok"), "/",
         nrow(batch_log), " sessions processed successfully.")
  
  ts_msg("Batch complete. ", sum(batch_log$status == "ok"), "/",
         nrow(batch_log), " sessions processed successfully.")
  
  ## --------------------------------------------------------------------------
  ## Expose key master tables to Global Environment for post-run inspection
  ## --------------------------------------------------------------------------
  assign("master_qc",           master_qc,           envir = .GlobalEnv)
  assign("master_body_by_fish", master_body_by_fish, envir = .GlobalEnv)
  assign("master_position",     master_position,     envir = .GlobalEnv)
  assign("master_speed_frame",  master_speed_frame,  envir = .GlobalEnv)
  assign("master_zone_frame",   master_zone_frame,   envir = .GlobalEnv)
  assign("master_gd_ind",       master_gd_ind,       envir = .GlobalEnv)
  assign("master_gd_pw",        master_gd_pw,        envir = .GlobalEnv)
  assign("master_gd_grp",       master_gd_grp,       envir = .GlobalEnv)
  assign("master_col_events",   master_col_events,   envir = .GlobalEnv)
  assign("master_evt_summary",  master_evt_summary,  envir = .GlobalEnv)
  assign("master_h_stab",       master_h_stab,       envir = .GlobalEnv)
  
  invisible(list(log = batch_log, files_index = files_index, out_dir = out_dir))
}


###############################################################################
##  15) RUN                                                                  ##
###############################################################################

## [E] Optional: to run as background job, uncomment and use:
##   rstudioapi::jobRunScript("idtracker_batch_pipeline.R")
## Or from terminal:
##   Rscript --vanilla idtracker_batch_pipeline.R

base_dir <- clipboard_dir()

batch_out <- run_batch_folder(
  base_dir         = base_dir,
  n_fish_max       = 4L,
  max_gap_seconds  = 1,
  speed_prob_min   = NULL,
  use_time_dt      = TRUE,
  prob_threshold   = 0.80,
  y0_px            = mean(c(483, 477)),
  flow_is_higher_y = TRUE,
  conda_exe        = CONDA_EXE,
  conda_env        = CONDA_ENV,
  script_path      = BBOX_SCRIPT,
  ensure_bboxes    = TRUE
)

batch_out$log

###############################################################################
##  16) POST-RUN PATCH: BUILD master_fish_by_frame                           ##
##  Runs AFTER run_batch_folder() using in-memory objects if available,       ##
##  otherwise reads already-written master CSVs from disk.                    ##
###############################################################################

build_master_fish_by_frame_postrun <- function(batch_out) {
  if (is.null(batch_out) || is.null(batch_out$out_dir) || !dir.exists(batch_out$out_dir)) {
    stop("batch_out$out_dir is missing or does not exist.")
  }
  
  out_dir <- normalizePath(batch_out$out_dir, winslash = "/", mustWork = TRUE)
  
  read_if_exists <- function(path) {
    if (file.exists(path)) {
      suppressMessages(readr::read_csv(path, show_col_types = FALSE))
    } else {
      NULL
    }
  }
  
  first_nonnull_df <- function(...) {
    xs <- list(...)
    for (x in xs) {
      if (is.data.frame(x)) return(x)
    }
    NULL
  }
  
  as_int_if_present <- function(df, cols) {
    for (cc in intersect(cols, names(df))) df[[cc]] <- as.integer(df[[cc]])
    df
  }
  
  as_num_if_present <- function(df, cols) {
    for (cc in intersect(cols, names(df))) df[[cc]] <- suppressWarnings(as.numeric(df[[cc]]))
    df
  }
  
  as_chr_if_present <- function(df, cols) {
    for (cc in intersect(cols, names(df))) df[[cc]] <- as.character(df[[cc]])
    df
  }
  
  ## --------------------------------------------------------------------------
  ## Load source tables
  ## --------------------------------------------------------------------------
  speed_by_frame <- first_nonnull_df(
    if (exists("master_speed_frame", envir = .GlobalEnv)) get("master_speed_frame", envir = .GlobalEnv) else NULL,
    read_if_exists(file.path(out_dir, "by_frame", "speed", "master_speed_by_frame.csv"))
  )
  
  zone_by_frame <- first_nonnull_df(
    if (exists("master_zone_frame", envir = .GlobalEnv)) get("master_zone_frame", envir = .GlobalEnv) else NULL,
    read_if_exists(file.path(out_dir, "by_frame", "occupancy", "master_zone_by_frame.csv"))
  )
  
  body_by_fish <- first_nonnull_df(
    if (exists("master_body_by_fish", envir = .GlobalEnv)) get("master_body_by_fish", envir = .GlobalEnv) else NULL,
    read_if_exists(file.path(out_dir, "morph", "master_bodylength_by_fish.csv"))
  )
  
  qc_by_fish <- first_nonnull_df(
    if (exists("master_qc", envir = .GlobalEnv)) get("master_qc", envir = .GlobalEnv) else NULL,
    read_if_exists(file.path(out_dir, "qc", "master_qc_by_fish.csv"))
  )
  
  if (is.null(speed_by_frame) || nrow(speed_by_frame) == 0L) {
    stop("Could not find master_speed_by_frame in memory or on disk.")
  }
  
  ## --------------------------------------------------------------------------
  ## Normalize types
  ## --------------------------------------------------------------------------
  speed_by_frame <- speed_by_frame %>%
    as_chr_if_present(c("session_id", "fish_id", "notes")) %>%
    as_int_if_present(c("trial_id", "segment_id", "fish_in_session", "individual", "frame")) %>%
    as_num_if_present(c(
      "time", "x_interp", "y_interp", "probability",
      "vx_cm_s", "vy_cm_s", "speed_cm_s",
      "heading_rad", "turn_angle_rad", "accel_cm_s2"
    ))
  
  if (!is.null(zone_by_frame) && nrow(zone_by_frame) > 0) {
    zone_by_frame <- zone_by_frame %>%
      as_chr_if_present(c("session_id", "fish_id", "zone")) %>%
      as_int_if_present(c("trial_id", "segment_id", "fish_in_session", "frame")) %>%
      as_num_if_present(c("time"))
    
    if ("valid" %in% names(zone_by_frame)) {
      zone_by_frame$valid <- as.logical(zone_by_frame$valid)
    }
  }
  
  if (!is.null(body_by_fish) && nrow(body_by_fish) > 0) {
    body_by_fish <- body_by_fish %>%
      as_chr_if_present(c("session_id", "fish_id")) %>%
      as_int_if_present(c("trial_id", "segment_id", "fish_in_session")) %>%
      as_num_if_present(c("body_length_cm", "body_length_sd_cm"))
  }
  
  if (!is.null(qc_by_fish) && nrow(qc_by_fish) > 0) {
    qc_by_fish <- qc_by_fish %>%
      as_chr_if_present(c("session_id", "fish_id")) %>%
      as_int_if_present(c("trial_id", "segment_id", "fish_in_session")) %>%
      as_num_if_present(c("length_unit_px_per_cm"))
  }
  
  fish_frame_keys <- c("session_id", "trial_id", "segment_id", "fish_id", "fish_in_session", "frame")
  fish_keys <- c("session_id", "trial_id", "segment_id", "fish_id", "fish_in_session")
  
  ## --------------------------------------------------------------------------
  ## Start from speed_by_frame
  ## --------------------------------------------------------------------------
  speed_keep <- intersect(c(
    "session_id",
    "trial_id",
    "segment_id",
    "fish_id",
    "fish_in_session",
    "individual",
    "frame",
    "time",
    "x_interp",
    "y_interp",
    "notes",
    "probability",
    "vx_cm_s",
    "vy_cm_s",
    "speed_cm_s",
    "heading_rad",
    "turn_angle_rad",
    "accel_cm_s2",
    "outlier_flag"
  ), names(speed_by_frame))

  master_fish_by_frame <- speed_by_frame[, speed_keep, drop = FALSE]
  
  ## --------------------------------------------------------------------------
  ## Join zone_by_frame
  ## New schema: valid + zone only (no prob_threshold, no boundary metadata)
  ## --------------------------------------------------------------------------
  if (!is.null(zone_by_frame) && nrow(zone_by_frame) > 0) {
    zone_keep <- intersect(c(
      fish_frame_keys,
      "valid",
      "zone"
    ), names(zone_by_frame))
    
    if (length(intersect(fish_frame_keys, zone_keep)) == length(fish_frame_keys)) {
      zone_small <- zone_by_frame[, zone_keep, drop = FALSE] %>%
        dplyr::distinct()
      
      master_fish_by_frame <- dplyr::left_join(
        master_fish_by_frame,
        zone_small,
        by = fish_frame_keys
      )
    }
  }
  
  ## --------------------------------------------------------------------------
  ## Join body length by fish
  ## --------------------------------------------------------------------------
  if (!is.null(body_by_fish) && nrow(body_by_fish) > 0) {
    body_keep <- intersect(c(
      fish_keys,
      "body_length_cm",
      "body_length_sd_cm"
    ), names(body_by_fish))
    
    if (length(intersect(fish_keys, body_keep)) == length(fish_keys)) {
      body_small <- body_by_fish[, body_keep, drop = FALSE] %>%
        dplyr::distinct()
      
      master_fish_by_frame <- dplyr::left_join(
        master_fish_by_frame,
        body_small,
        by = fish_keys
      )
    }
  }
  
  ## --------------------------------------------------------------------------
  ## Compute locomotor_status_TS{0.5,1,1.5,2} from speed and body length
  ## A fish is "moving" when speed_cm_s > body_length_cm × threshold_multiplier.
  ## --------------------------------------------------------------------------
  if ("body_length_cm" %in% names(master_fish_by_frame) &&
      "speed_cm_s"     %in% names(master_fish_by_frame)) {
    for (.thr in c(0.5, 1.0, 1.5, 2.0)) {
      .label  <- paste0("TS", as.character(.thr))   # "TS0.5", "TS1", "TS1.5", "TS2"
      .st_col <- paste0("speed_threshold_", .label)
      .ls_col <- paste0("locomotor_status_", .label)
      master_fish_by_frame[[.st_col]] <- master_fish_by_frame$body_length_cm * .thr
      master_fish_by_frame[[.ls_col]] <- dplyr::case_when(
        is.na(master_fish_by_frame$speed_cm_s)           ~ NA_character_,
        is.na(master_fish_by_frame[[.st_col]])             ~ NA_character_,
        master_fish_by_frame$speed_cm_s >
          master_fish_by_frame[[.st_col]]                 ~ "moving",
        TRUE                                               ~ "immobile"
      )
    }
    rm(.thr, .label, .st_col, .ls_col)
  }

  ## --------------------------------------------------------------------------
  ## Join length_unit_px_per_cm from qc_by_fish
  ## --------------------------------------------------------------------------
  if (!is.null(qc_by_fish) && nrow(qc_by_fish) > 0) {
    qc_keep <- intersect(c(
      fish_keys,
      "length_unit_px_per_cm"
    ), names(qc_by_fish))
    
    if (length(intersect(fish_keys, qc_keep)) == length(fish_keys)) {
      qc_small <- qc_by_fish[, qc_keep, drop = FALSE] %>%
        dplyr::distinct()
      
      master_fish_by_frame <- dplyr::left_join(
        master_fish_by_frame,
        qc_small,
        by = fish_keys
      )
    }
  }
  
  ## --------------------------------------------------------------------------
  ## Enforce final requested schema
  ## --------------------------------------------------------------------------
  target_cols <- c(
    "session_id",
    "trial_id",
    "segment_id",
    "fish_id",
    "fish_in_session",
    "individual",
    "frame",
    "time",
    "x_interp",
    "y_interp",
    "notes",
    "probability",
    "vx_cm_s",
    "vy_cm_s",
    "speed_cm_s",
    "heading_rad",
    "turn_angle_rad",
    "accel_cm_s2",
    "valid",
    "zone",
    "length_unit_px_per_cm",
    "body_length_cm",
    "body_length_sd_cm",
    "speed_threshold_TS0.5", "locomotor_status_TS0.5",
    "speed_threshold_TS1",   "locomotor_status_TS1",
    "speed_threshold_TS1.5", "locomotor_status_TS1.5",
    "speed_threshold_TS2",   "locomotor_status_TS2",
    "outlier_flag"
  )
  
  missing_cols <- setdiff(target_cols, names(master_fish_by_frame))
  if (length(missing_cols) > 0) {
    for (cc in missing_cols) master_fish_by_frame[[cc]] <- NA
  }
  
  master_fish_by_frame <- master_fish_by_frame %>%
    dplyr::select(dplyr::all_of(target_cols)) %>%
    dplyr::arrange(session_id, fish_in_session, frame)
  
  out_path <- file.path(out_dir, "by_frame", "master_fish_by_frame.csv")
  readr::write_csv(master_fish_by_frame, out_path)
  
  assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)
  
  message(
    format(Sys.time(), "%H:%M:%S"),
    " | Wrote: ", out_path,
    " [", nrow(master_fish_by_frame), " rows x ", ncol(master_fish_by_frame), " cols]"
  )
  message(
    format(Sys.time(), "%H:%M:%S"),
    " | Loaded into Global Environment: master_fish_by_frame"
  )
  
  invisible(master_fish_by_frame)
}

## ---- RUN THE PATCH ONLY ----
master_fish_by_frame <- build_master_fish_by_frame_postrun(batch_out)

###############################################################################
##  17) JUMP DETECTION (RETIRED 2026-08-20 — see STEP 3.5)                   ##
###############################################################################
# This block used to run its OWN jump detection here, on master_fish_by_frame,
# ~700 lines after speed_cm_s had already been computed from the same
# uncorrected coordinates — so any point it nulled left speed_cm_s stale
# (verified: 47/47 corrected frames in run STEP2_output_20260524_160853 kept
# pre-correction speeds up to 1,008,298 cm/s), and its mean+SD threshold was
# itself inflated by the extreme frames it needed to catch (pooled SD 89
# cm/s -> a 906 cm/s threshold that a 26,677 cm/s frame sailed under). See
# LOGS.md 2026-08-20 and the JUMP_DETECTION_ENABLED comment in section 1b.
#
# Rejection now happens in STEP 3.5 (dt_reject_outliers_inplace), per
# session, on raw x/y, BEFORE interpolation and BEFORE speed is computed —
# so the correction actually reaches every downstream kinematic. This block
# is kept only to preserve the master_fish_by_frame schema that downstream
# scripts (02_individual_indicators.R and others) already validate
# against: x_interp_raw, y_interp_raw, jump_flag.
###############################################################################

.jd_enabled <- isTRUE(if (exists("JUMP_DETECTION_ENABLED", inherits = FALSE))
                         JUMP_DETECTION_ENABLED else FALSE)
if (.jd_enabled) {
  stop("JUMP_DETECTION_ENABLED <- TRUE is no longer supported: the legacy ",
       "mean+SD algorithm was retired 2026-08-20 (see the comment above and ",
       "LOGS.md) because it could not fix the artifacts it was built to ",
       "catch. Use STEP 3.5 / OUTLIER_REJECT_* instead.", call. = FALSE)
}

# x_interp at this point already reflects STEP 3.5's rejections plus STEP 04's
# interpolation over them, so "raw" here is trivially identical to it — kept
# only for schema compatibility with anything that reads x_interp_raw.
master_fish_by_frame$x_interp_raw <- master_fish_by_frame$x_interp
master_fish_by_frame$y_interp_raw <- master_fish_by_frame$y_interp
master_fish_by_frame$jump_flag    <- if ("outlier_flag" %in% names(master_fish_by_frame)) {
  ifelse(is.na(master_fish_by_frame$outlier_flag), FALSE, master_fish_by_frame$outlier_flag)
} else {
  FALSE
}
assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)
assign("jump_detection_summary", data.frame(),        envir = .GlobalEnv)

ts_msg("Jump detection (retired): jump_flag populated from STEP 3.5 outlier_flag (",
       sum(master_fish_by_frame$jump_flag), " of ", nrow(master_fish_by_frame), " frames).")

## optional quick check
dplyr::glimpse(master_fish_by_frame)