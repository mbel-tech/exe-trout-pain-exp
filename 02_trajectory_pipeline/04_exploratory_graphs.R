# =============================================================================
# Load group + event master outputs from master_outputs/ (path from clipboard)
# - Windows-safe clipboard path reading + normalization
# - Imports 4 CSVs into the Global Environment
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(officer)
  library(rvg)
  library(readr)
  library(tibble)
  library(purrr)
  library(future)
  library(future.apply)
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

# ---- pretty labels (R3: centralized — sourced by master) ----
# pretty_labels.R is sourced once by 00_run_trajectory_pipeline.R
# before any step runs, so `pretty_indicator_name()` and `GOT_INDICATOR_LABELS`
# are normally already in .GlobalEnv when this script starts.
#
# Fallback below is only for the rare case of running this script standalone
# (i.e. not via master). It looks once in this script's own directory.
if (!exists("pretty_indicator_name", mode = "function")) {
  .this_script_dir <- tryCatch({
    frames <- sys.frames()
    ofile_envs <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
    if (length(ofile_envs) > 0L) {
      normalizePath(dirname(get("ofile", envir = ofile_envs[[length(ofile_envs)]])),
                    winslash = "/", mustWork = FALSE)
    } else getwd()
  }, error = function(e) getwd())

  .pl_path <- file.path(.this_script_dir, "helpers/pretty_labels.R")
  if (file.exists(.pl_path)) {
    source(.pl_path, local = FALSE)
  } else {
    message("Note: pretty_labels.R not found at ", .pl_path,
            " — falling back to pretty_var_label() if available.")
    if (exists("pretty_var_label", mode = "function")) {
      pretty_indicator_name <- pretty_var_label
    } else {
      pretty_indicator_name <- function(x) gsub("_", " ", as.character(x)[1])
    }
  }
  rm(.this_script_dir, .pl_path)
}

# --- helpers --------------------------a-------------------------------------'

read_csv_optional <- function(path, object_name = NULL) {
  if (!file.exists(path)) {
    warning("Optional file not found", 
            if (!is.null(object_name)) paste0(" for `", object_name, "`"), 
            ": ", path,
            call. = FALSE)
    return(NULL)
  }
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

norm_path_win <- function(p, mustWork = TRUE) {
  p <- trimws(p)
  p <- gsub("\r", "", p)
  p <- gsub('^"|"$', "", p)
  p <- gsub("\\\\", "/", p)
  normalizePath(p, winslash = "/", mustWork = mustWork)
}

clipboard_dir <- function() {
  if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
    p <- get("PIPELINE_DATA_DIR", envir = .GlobalEnv)
    return(norm_path_win(p, mustWork = FALSE))
  }
  default_dir <- "data/tracking_sessions"
  if (dir.exists(default_dir)) {
    message("clipboard_dir: using default data directory: ", default_dir)
    return(norm_path_win(default_dir, mustWork = FALSE))
  }
  p <- readLines("clipboard", warn = FALSE)
  p <- trimws(p)
  p <- p[nzchar(p)]
  if (length(p) == 0) stop("Clipboard is empty. Copy the master_outputs folder path and retry.", call. = FALSE)
  p <- p[1]
  p <- norm_path_win(p, mustWork = TRUE)
  if (!dir.exists(p)) stop("Clipboard path is not a directory: ", p, call. = FALSE)
  p
}

read_csv_safe <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

# --- 1) base folder from clipboard -----------------------------------------

master_outputs_dir <- clipboard_dir()

# If user copied the *parent* directory, try to descend into master_outputs/
if (!grepl("(?i)/master_outputs$", master_outputs_dir)) {
  candidate <- file.path(master_outputs_dir, "master_outputs")
  candidate <- norm_path_win(candidate, mustWork = FALSE)
  if (dir.exists(candidate)) master_outputs_dir <- candidate
}

if (!dir.exists(master_outputs_dir)) {
  stop("master_outputs folder not found. Clipboard path resolved to: ", master_outputs_dir, call. = FALSE)
}

message("Using master_outputs dir: ", master_outputs_dir)

# --- 2) optionally read group/events outputs --------------------------------

group_dir  <- file.path(master_outputs_dir, "group")
events_dir <- file.path(master_outputs_dir, "events")

path_gdi <- file.path(group_dir,  "master_group_dynamics_ind_by_frame.csv")
path_gd  <- file.path(group_dir,  "master_group_dynamics_by_frame.csv")

path_ev  <- file.path(events_dir, "master_col_mot_events_per_session.csv")
path_sum <- file.path(events_dir, "master_session_event_summary.csv")

master_group_dynamics_ind_by_frame <- if (dir.exists(group_dir)) {
  read_csv_optional(path_gdi, "master_group_dynamics_ind_by_frame")
} else {
  warning("Optional folder missing: ", group_dir, call. = FALSE)
  NULL
}

master_group_dynamics_by_frame <- if (dir.exists(group_dir)) {
  read_csv_optional(path_gd, "master_group_dynamics_by_frame")
} else {
  NULL
}

master_col_mot_events_per_session <- if (dir.exists(events_dir)) {
  read_csv_optional(path_ev, "master_col_mot_events_per_session")
} else {
  warning("Optional folder missing: ", events_dir, call. = FALSE)
  NULL
}

master_session_event_summary <- if (dir.exists(events_dir)) {
  read_csv_optional(path_sum, "master_session_event_summary")
} else {
  NULL
}

# --- 3) put available objects into Global Environment -----------------------

if (!is.null(master_group_dynamics_ind_by_frame)) {
  assign("master_group_dynamics_ind_by_frame", master_group_dynamics_ind_by_frame, envir = .GlobalEnv)
  message("Imported: master_group_dynamics_ind_by_frame (", nrow(master_group_dynamics_ind_by_frame), " rows)")
} else {
  message("Skipped: master_group_dynamics_ind_by_frame")
}

if (!is.null(master_group_dynamics_by_frame)) {
  assign("master_group_dynamics_by_frame", master_group_dynamics_by_frame, envir = .GlobalEnv)
  message("Imported: master_group_dynamics_by_frame (", nrow(master_group_dynamics_by_frame), " rows)")
} else {
  message("Skipped: master_group_dynamics_by_frame")
}

if (!is.null(master_col_mot_events_per_session)) {
  assign("master_col_mot_events_per_session", master_col_mot_events_per_session, envir = .GlobalEnv)
  message("Imported: master_col_mot_events_per_session (", nrow(master_col_mot_events_per_session), " rows)")
} else {
  message("Skipped: master_col_mot_events_per_session")
}

if (!is.null(master_session_event_summary)) {
  assign("master_session_event_summary", master_session_event_summary, envir = .GlobalEnv)
  message("Imported: master_session_event_summary (", nrow(master_session_event_summary), " rows)")
} else {
  message("Skipped: master_session_event_summary")
}

################################################################################

# =============================================================================
# Assumption: `set` looks like "session_trialX_SegmentY"
# where X = trial_id and Y = segment_id (your timepoint key)
# This script:
# 1) loads trials summary GOT.xlsx as `trial_summary`
# 2) builds meta table (Trial + timepoint -> Group/Treatment/Date)
# 3) parses trial_id + segment_id from `set`
# 4) annotates the 4 group/events datasets in GlobalEnv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readxl)
})

# -----------------------------------------------------------------------------
# 0) GLOBAL SETTINGS
# -----------------------------------------------------------------------------
PT_ALPHA     <- 1.00
PT_SIZE      <- 3.6
JITTER_WIDTH <- 0.14
DODGE_WIDTH  <- 0.75

MEAN_CROSSBAR_WIDTH     <- 0.55
MEAN_CROSSBAR_LINE_SIZE <- 0.90
ERRORBAR_WIDTH          <- 0.18
ERRORBAR_LINE_SIZE      <- 0.70

FACET_NROW <- 1
PANEL_SPACING_LINES <- 1.2

OUT_P_LO <- 0.01
OUT_P_HI <- 0.99

AXIS_LINE_WIDTH <- 0.85
CAPTION_TEXT_SIZE <- 10.5

TREND_LINE_SIZE <- 1.2
TREND_POINT_SIZE <- 3.8
TREND_ERRORBAR_WIDTH <- 0.10

GROUP_COLORS <- c(
  "No choice"       = "#56B4E9",
  "Exercise choice" = "#D55E00",
  "Morphine"        = "#8DBB1F",
  "Naloxone"        = "#C084F0"
)

GROUP_ORDER <- c("No choice", "Exercise choice", "Naloxone", "Morphine")
TREATMENT_ORDER <- c("Control", "Exercise")

timestamp_tag <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
TS_NOW <- timestamp_tag()
# MAIN_EXPORT_DIR and SUMMARY_CSV_PATH are set in the main execution section below.

POS_JD <- position_jitterdodge(
  jitter.width  = JITTER_WIDTH,
  jitter.height = 0,
  dodge.width   = DODGE_WIDTH
)

POS_DODGE <- position_dodge(width = DODGE_WIDTH)

BASE_THEME <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    axis.line.x = element_line(color = "black", linewidth = 1.05),
    axis.line.y = element_line(color = "black", linewidth = 1.05),
    axis.ticks  = element_line(color = "black", linewidth = 0.8),
    
    axis.text.x   = element_text(size = 16, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 12, color = "black"),
    axis.title.x  = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y  = element_text(size = 16, face = "bold", color = "black"),
    strip.text    = element_text(size = 14, face = "bold", color = "black"),
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 12),
    
    panel.spacing = grid::unit(PANEL_SPACING_LINES, "lines"),
    plot.margin   = ggplot2::margin(t = 8, r = 10, b = 18, l = 10),
    
    plot.caption = element_text(
      size = CAPTION_TEXT_SIZE,
      hjust = 0,
      margin = margin(t = 18),
      color = "black",
      lineheight = 1.05
    ),
    plot.caption.position = "plot"
  )

# branch_specs / subbranch_specs are defined later (after data load) so they
# can respond to the upstream KEEP_TP1 / DROP_PHARMA switches.  These early
# placeholder values ensure downstream code that references them before the
# real definitions doesn't crash; they are overwritten unconditionally below.
branch_specs    <- list()
subbranch_specs <- list()

# -----------------------------------------------------------------------------
# 1) Load trial summary to GlobalEnv
# -----------------------------------------------------------------------------
TRIAL_SUMMARY_PATH <- if (exists("TRIAL_SUMMARY_XLSX", envir = .GlobalEnv, inherits = FALSE)) {
  get("TRIAL_SUMMARY_XLSX", envir = .GlobalEnv)
} else {
  resolve_dataset_path(
    "data/trials summary GOT.xlsx"
  )
}

trial_summary <- readxl::read_excel(TRIAL_SUMMARY_PATH)
assign("trial_summary", trial_summary, envir = .GlobalEnv)

# -----------------------------------------------------------------------------
# 2) Helpers
# -----------------------------------------------------------------------------
norm_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("\\s+", "", x)
  x
}

# Parse "session_trialX_SegmentY" robustly (case-insensitive, flexible separators/spaces)
extract_trial_segment_from_set <- function(set_vec) {
  set_vec <- as.character(set_vec)
  m <- stringr::str_match(set_vec, "(?i)trial\\s*(\\d+).*segment\\s*(\\d+)")
  tibble(
    trial_id   = suppressWarnings(as.integer(m[, 2])),
    segment_id = suppressWarnings(as.integer(m[, 3]))
  )
}

# -----------------------------------------------------------------------------
# 3) Meta builder + annotator (your original logic, unchanged)
# -----------------------------------------------------------------------------
build_trial_meta_robust <- function(trial_summary,
                                    meta_trial_col = "Trial",
                                    meta_time_col  = "timepoint",
                                    meta_cols      = c("Group", "Treatment", "Date")) {
  required <- c(meta_trial_col, meta_time_col, meta_cols)
  missing <- setdiff(required, names(trial_summary))
  if (length(missing) > 0) stop("trial_summary missing required cols: ", paste(missing, collapse = ", "), call. = FALSE)
  
  trial_summary %>%
    transmute(
      .trial_key = norm_key(.data[[meta_trial_col]]),
      .time_key  = norm_key(.data[[meta_time_col]]),
      Group      = as.character(.data[["Group"]]),
      Treatment  = as.character(.data[["Treatment"]]),
      trial_date = .data[["Date"]]
    ) %>%
    filter(nzchar(.trial_key), nzchar(.time_key)) %>%
    distinct(.trial_key, .time_key, .keep_all = TRUE)
}

annotate_with_meta_robust <- function(df,
                                      meta_tbl,
                                      df_trial_col = "trial_id",
                                      df_time_col_candidates = c("segment_id", "timepoint"),
                                      rename_time_to = "timepoint") {
  if (is.null(df)) return(NULL)
  
  time_col <- df_time_col_candidates[df_time_col_candidates %in% names(df)][1]
  if (!(df_trial_col %in% names(df)) || is.na(time_col)) {
    warning("Cannot annotate: df missing trial/time keys (trial_id + segment_id/timepoint).", call. = FALSE)
    return(df)
  }
  
  out <- df %>%
    mutate(
      .trial_key = norm_key(.data[[df_trial_col]]),
      .time_key  = norm_key(.data[[time_col]])
    )
  
  if (time_col != rename_time_to && !(rename_time_to %in% names(out))) {
    out <- rename(out, !!rename_time_to := !!time_col)
  }
  
  out %>%
    left_join(meta_tbl, by = c(".trial_key" = ".trial_key", ".time_key" = ".time_key")) %>%
    select(-.trial_key, -.time_key)
}

# -----------------------------------------------------------------------------
# 4) Build meta
# -----------------------------------------------------------------------------
meta_tbl <- build_trial_meta_robust(
  trial_summary = trial_summary,
  meta_trial_col = "Trial",
  meta_time_col  = "timepoint",
  meta_cols      = c("Group", "Treatment", "Date")
)

# -----------------------------------------------------------------------------
# 5) Ensure trial_id + segment_id exist by parsing from `set`,
#    then annotate each dataset in GlobalEnv (overwrites objects)
# -----------------------------------------------------------------------------
ensure_keys_from_set <- function(df) {
  if (is.null(df) || !is.data.frame(df)) return(df)
  if (!("set" %in% names(df))) return(df)
  
  # If already present, keep; else parse
  if (!("trial_id" %in% names(df)) || !("segment_id" %in% names(df))) {
    ts <- extract_trial_segment_from_set(df$set)
    if (!("trial_id" %in% names(df)))   df$trial_id   <- ts$trial_id
    if (!("segment_id" %in% names(df))) df$segment_id <- ts$segment_id
  }
  df
}

annotate_object_in_global <- function(obj_name) {
  if (!exists(obj_name, envir = .GlobalEnv)) {
    warning("Global object not found: ", obj_name, call. = FALSE)
    return(invisible(NULL))
  }
  df <- get(obj_name, envir = .GlobalEnv)
  df <- ensure_keys_from_set(df)
  df2 <- annotate_with_meta_robust(
    df,
    meta_tbl,
    df_trial_col = "trial_id",
    df_time_col_candidates = c("segment_id", "timepoint"),
    rename_time_to = "timepoint"
  )
  assign(obj_name, df2, envir = .GlobalEnv)
  invisible(df2)
}

optional_objects_to_annotate <- c(
  "master_group_dynamics_ind_by_frame",
  "master_group_dynamics_by_frame",
  "master_col_mot_events_per_session",
  "master_session_event_summary"
)

invisible(lapply(optional_objects_to_annotate, annotate_object_in_global))

message("Done: parsed trial_id/segment_id from `set` and annotated Group/Treatment/trial_date.")

######################################################################################

################################################################################
# =============================================================================
# ACTIVITY INDICATORS GRAPH PATCH — CONSOLIDATED / CORRECTED VERSION
# INPUT REQUIRED IN GLOBAL ENV:
#   master_fish_by_frame
#
# EXPECTED COLUMNS IN master_fish_by_frame:
#   session_id, trial_id, timepoint, fish_id, fish_in_session, individual,
#   frame, time, x_interp, y_interp, notes, probability, vx_cm_s, vy_cm_s,
#   speed_cm_s, heading_rad, turn_angle_rad, accel_cm_s2, valid, zone,
#   length_unit_px_per_cm, Group, Treatment, trial_date,
#   speed_threshold_TS0.5, locomotor_status_TS0.5,
#   speed_threshold_TS1, locomotor_status_TS1,
#   speed_threshold_TS1.5, locomotor_status_TS1.5,
#   speed_threshold_TS2, locomotor_status_TS2,
#   body_length_cm, body_length_sd_cm, speed_BL_s
#
# OUTPUTS:
#   activity_indicators_graphs_<timestamp>/
#     - fish_activity_summary.csv
#     - All_groups/
#     - Exe_choice/
#
# FEATURES:
# - Standard scatter graphs with:
#   * x = Treatment
#   * points colored by Group
#   * mean crossbar + SEM
#   * no grid
#   * thick black x/y axis lines
#   * captions including sample-size / outlier note
# - Additional timepoint-trend graphs for NON-centered indicators:
#   * x = timepoint
#   * one mean point per Treatment × Group × timepoint
#   * SEM error bars
#   * color = Group
#   * shape = Treatment
#   * linetype = Treatment (Control dashed, Exercise solid)
#   * caption = significant Pearson correlations only
# - Robust label normalization for Group / Treatment
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(officer)
  library(rvg)
  library(readr)
  library(tibble)
  library(purrr)
  library(cowplot)
})

# -----------------------------------------------------------------------------
# 0) GLOBAL SETTINGS
# -----------------------------------------------------------------------------
PT_ALPHA     <- 0.75
PT_SIZE      <- 3.6
JITTER_WIDTH <- 0.14
DODGE_WIDTH  <- 0.75

MEAN_CROSSBAR_WIDTH     <- 0.55
MEAN_CROSSBAR_LINE_SIZE <- 0.90
ERRORBAR_WIDTH          <- 0.18
ERRORBAR_LINE_SIZE      <- 0.70

TREND_LINE_SIZE       <- 1.2
TREND_POINT_SIZE      <- 3.8
TREND_ERRORBAR_WIDTH  <- 0.10

FACET_NROW <- 1
PANEL_SPACING_LINES <- 1.2

OUT_P_LO <- 0.01
OUT_P_HI <- 0.99

N_LABEL_SIZE   <- 3.8
N_LABEL_PREFIX <- "N="

CAPTION_TEXT_SIZE <- 10.5

ANOVA_CAPTION_TEXT_SIZE <- 13.5
ANOVA_CAPTION_LINEHEIGHT <- 1.08

GROUP_COLORS <- c(
  "No choice"       = "#56B4E9",
  "Exercise choice" = "#D55E00",
  "Morphine"        = "#8DBB1F",
  "Naloxone"        = "#C084F0"
)

GROUP_ORDER <- c("No choice", "Exercise choice", "Naloxone", "Morphine")
TREATMENT_ORDER <- c("Control", "Exercise")

# OPT-0: build branch list dynamically based on upstream filter switches.
# When DROP_PHARMA=TRUE pharma groups are already absent from the data, so the
# "All_groups" branch is redundant — only keep "Exe_choice" (Exercise choice +
# No choice).  When DROP_PHARMA is FALSE (or absent) both branches run.
.drop_pharma_active <- exists("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE) &&
  isTRUE(get("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE))

branch_specs <- if (.drop_pharma_active) {
  message("OPT-0 (STEP4): DROP_PHARMA=TRUE -> running Exe_choice branch only")
  list(
    list(
      branch_name  = "Exe_choice",
      group_filter = c("Exercise choice", "No choice")
    )
  )
} else {
  list(
    list(
      branch_name  = "All_groups",
      group_filter = NULL
    ),
    list(
      branch_name  = "Exe_choice",
      group_filter = c("Exercise choice", "No choice")
    )
  )
}

# OPT-0: exclude Timepoint1 subbranch when KEEP_TP1=FALSE because those rows
# have been removed from the data upstream.
.keep_tp1_active <- !exists("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE) ||
  isTRUE(get("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE))

subbranch_specs <- if (.keep_tp1_active) {
  list(
    list(subbranch_name = "All_timepoints", timepoint_filter = NULL),
    list(subbranch_name = "Timepoints_23",  timepoint_filter = c(2, 3)),
    list(subbranch_name = "Timepoint1",     timepoint_filter = c(1))
  )
} else {
  message("OPT-0 (STEP4): KEEP_TP1=FALSE -> dropping Timepoint1 subbranch")
  list(
    list(subbranch_name = "All_timepoints", timepoint_filter = NULL),
    list(subbranch_name = "Timepoints_23",  timepoint_filter = c(2, 3))
  )
}

timestamp_tag <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
TS_NOW <- timestamp_tag()
MAIN_EXPORT_DIR  <- .got_make_output_dir("STEP4_output")
SUMMARY_CSV_PATH <- file.path(MAIN_EXPORT_DIR, "fish_activity_summary.csv")

POS_JD <- position_jitterdodge(
  jitter.width  = JITTER_WIDTH,
  jitter.height = 0,
  dodge.width   = DODGE_WIDTH
)

POS_DODGE <- position_dodge(width = DODGE_WIDTH)

BASE_THEME <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 1.05),
    axis.line.y = element_line(color = "black", linewidth = 1.05),
    axis.ticks  = element_line(color = "black", linewidth = 0.8),
    axis.text.x   = element_text(size = 16, face = "bold", color = "black"),
    axis.text.y   = element_text(size = 12, color = "black"),
    axis.title.x  = element_text(size = 16, face = "bold", color = "black"),
    axis.title.y  = element_text(size = 16, face = "bold", color = "black"),
    strip.text    = element_text(size = 14, face = "bold", color = "black"),
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 12),
    panel.spacing = grid::unit(PANEL_SPACING_LINES, "lines"),
    plot.margin   = ggplot2::margin(t = 8, r = 10, b = 18, l = 10),
    plot.caption  = element_text(
      size = CAPTION_TEXT_SIZE,
      hjust = 0,
      margin = margin(t = 18),
      color = "black",
      lineheight = 1.05
    ),
    plot.caption.position = "plot"
  )

# -----------------------------------------------------------------------------
# 1) CHECK INPUT
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# 1) CHECK INPUT
#    - use existing GlobalEnv object if present
#    - otherwise search recursively under clipboard-resolved folder
# -----------------------------------------------------------------------------
find_master_fish_file <- function(base_dir) {
  if (!dir.exists(base_dir)) {
    stop("Base directory does not exist: ", base_dir, call. = FALSE)
  }
  
  # exact filename first
  hits_exact <- list.files(
    path = base_dir,
    pattern = "^master_fish_by_frame\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(hits_exact) > 0) {
    return(normalizePath(hits_exact[1], winslash = "/", mustWork = TRUE))
  }
  
  # fallback: any csv containing the target string
  hits_loose <- list.files(
    path = base_dir,
    pattern = "master_fish_by_frame.*\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(hits_loose) > 0) {
    return(normalizePath(hits_loose[1], winslash = "/", mustWork = TRUE))
  }
  
  NA_character_
}

load_master_fish_from_disk <- function(base_dir) {
  fish_file <- find_master_fish_file(base_dir)
  
  if (!is.character(fish_file) || is.na(fish_file) || !nzchar(fish_file)) {
    stop(
      "Could not find `master_fish_by_frame.csv` anywhere under: ",
      base_dir,
      call. = FALSE
    )
  }
  
  message("Loading master_fish_by_frame from: ", fish_file)
  
  df <- readr::read_csv(
    fish_file,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  assign("master_fish_by_frame", df, envir = .GlobalEnv)
  df
}

if (exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
  master_fish_by_frame <- get("master_fish_by_frame", envir = .GlobalEnv)
  message("Using existing GlobalEnv object: master_fish_by_frame (", nrow(master_fish_by_frame), " rows)")
} else {
  # Priority 1: canonical master_outputs/by_frame (persistent store, not invalidated by partial STEP2 re-runs)
  master_outputs_by_frame_file <- file.path(master_outputs_dir, "by_frame", "master_fish_by_frame.csv")
  if (!file.exists(master_outputs_by_frame_file)) {
    master_outputs_by_frame_file <- file.path(master_outputs_dir, "master_fish_by_frame.csv")
  }
  if (file.exists(master_outputs_by_frame_file)) {
    message("Loading master_fish_by_frame from master_outputs: ", master_outputs_by_frame_file)
    master_fish_by_frame <- readr::read_csv(master_outputs_by_frame_file, show_col_types = FALSE, progress = FALSE)
    assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)
  }

  # Priority 2: most recent STEP2_output (fallback when master_outputs is unavailable)
  if (!exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
    step2_dir <- file.path(dirname(dirname(master_outputs_dir)), "STEP2_output")
    if (dir.exists(step2_dir)) {
      step2_candidates <- list.files(step2_dir, pattern = "^STEP2_output_", full.names = TRUE)
      if (length(step2_candidates) > 0) {
        step2_latest <- step2_candidates[which.max(file.info(step2_candidates)$mtime)]
        step2_fish_file <- list.files(step2_latest, pattern = "^master_fish_by_frame\\.csv$", full.names = TRUE, recursive = TRUE)[1]
        if (!is.na(step2_fish_file)) {
          message("Loading master_fish_by_frame from STEP2_output: ", step2_fish_file)
          master_fish_by_frame <- readr::read_csv(step2_fish_file, show_col_types = FALSE, progress = FALSE)
          assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)
        }
      }
    }
  }

  # Priority 3: search under clipboard path parent
  if (!exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
    search_root <- master_outputs_dir

    master_fish_by_frame <- tryCatch(
      load_master_fish_from_disk(search_root),
      error = function(e1) {
        parent_root <- normalizePath(file.path(search_root, ".."), winslash = "/", mustWork = FALSE)

        if (!identical(parent_root, search_root) && dir.exists(parent_root)) {
          message("Not found under master_outputs; trying parent folder: ", parent_root)
          return(load_master_fish_from_disk(parent_root))
        }

        stop(e1$message, call. = FALSE)
      }
    )
  }
}

required_cols <- c(
  "trial_id", "timepoint", "fish_id", "frame", "time",
  "x_interp", "y_interp",
  "speed_cm_s", "heading_rad", "turn_angle_rad", "accel_cm_s2", "speed_BL_s",
  "length_unit_px_per_cm", "body_length_cm", "zone",
  "Group", "Treatment", "trial_date",
  "locomotor_status_TS0.5", "locomotor_status_TS1",
  "locomotor_status_TS1.5", "locomotor_status_TS2"
)

missing_cols <- setdiff(required_cols, names(master_fish_by_frame))
if (length(missing_cols) > 0) {
  stop(
    "master_fish_by_frame is missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# OPT-0 downstream row filter ------------------------------------------------
{
  .n0 <- nrow(master_fish_by_frame)
  if (exists("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE) &&
      !isTRUE(get("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE))) {
    master_fish_by_frame <- master_fish_by_frame %>%
      dplyr::filter(as.integer(timepoint) != 1L)
    message("OPT-0 (STEP4): removed ", .n0 - nrow(master_fish_by_frame),
            " TP1 rows; ", nrow(master_fish_by_frame), " remain")
    .n0 <- nrow(master_fish_by_frame)
  }
  if (exists("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE) &&
      isTRUE(get("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE)) &&
      exists("PHARMA_GROUPS", envir = .GlobalEnv, inherits = FALSE)) {
    .pg <- get("PHARMA_GROUPS", envir = .GlobalEnv, inherits = FALSE)
    master_fish_by_frame <- master_fish_by_frame %>%
      dplyr::filter(!Group %in% .pg)
    message("OPT-0 (STEP4): removed ", .n0 - nrow(master_fish_by_frame),
            " pharma rows (Group %in% ", paste(.pg, collapse = ","),
            "); ", nrow(master_fish_by_frame), " remain")
  }
  assign("master_fish_by_frame", master_fish_by_frame, envir = .GlobalEnv)
  rm(.n0)
}

# -----------------------------------------------------------------------------
# 2) HELPERS
# -----------------------------------------------------------------------------

has_at_least_two_timepoints <- function(df) {
  tp <- unique(as.integer(as.character(df$timepoint)))
  tp <- tp[is.finite(tp)]
  length(tp) >= 2
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_quantile <- function(x, probs) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(
    x,
    probs  = probs,
    na.rm  = TRUE,
    names  = FALSE,
    type   = 7
  ))
}

calc_dt <- function(time_vec) {
  d <- c(NA_real_, diff(time_vec))
  d[!is.finite(d) | d < 0] <- NA_real_
  d
}

is_moving_flag <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y %in% c("moving", "move", "1", "true", "yes")
}

sanitize_zone_name <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(trimws(x))] <- "missing"
  x <- trimws(x)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[nchar(x) == 0] <- "missing"
  x
}

normalize_group_label <- function(x) {
  y <- trimws(as.character(x))
  y_low <- tolower(y)
  
  dplyr::case_when(
    y_low %in% c("exercise choice", "exercise_choice") ~ "Exercise choice",
    y_low %in% c("no choice", "no_choice") ~ "No choice",
    y_low %in% c("morphine") ~ "Morphine",
    y_low %in% c("naloxone") ~ "Naloxone",
    TRUE ~ y
  )
}

normalize_treatment_label <- function(x) {
  y <- trimws(as.character(x))
  y_low <- tolower(y)
  
  dplyr::case_when(
    y_low %in% c("control", "ctrl") ~ "Control",
    y_low %in% c("exercise", "swim", "training") ~ "Exercise",
    TRUE ~ y
  )
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x)
}

safe_sem <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}

pretty_var_label <- function(v) {
  map <- c(
    distance_m                  = "Distance travelled (m)",
    distance_BL                 = "Distance travelled (body lengths)",
    
    speed_cm_s_mean            = "Mean speed (cm/s)",
    speed_cm_s_median          = "Median speed (cm/s)",
    speed_cm_s_max             = "Max speed (cm/s)",
    
    speed_BL_s_mean            = "Mean speed (BL/s)",
    speed_BL_s_median          = "Median speed (BL/s)",
    speed_BL_s_max             = "Max speed (BL/s)",
    
    accel_cm_s2_mean           = "Mean acceleration (cm/s²)",
    accel_cm_s2_median         = "Median acceleration (cm/s²)",
    accel_cm_s2_max            = "Max acceleration (cm/s²)",
    
    time_moving_TS0.5_s        = "Time moving TS0.5 (s)",
    time_moving_TS1_s          = "Time moving TS1 (s)",
    time_moving_TS1.5_s        = "Time moving TS1.5 (s)",
    time_moving_TS2_s          = "Time moving TS2 (s)",
    
    n_zone_transitions         = "Zone transitions (count)",
    
    flow_time_centered_300_s         = "Net time* spent in flow (s)",
    time_moving_centered_TS0_5_300_s = "Net time* spent moving (s)",
    time_moving_centered_TS1_300_s   = "Net time* spent moving (s)",
    time_moving_centered_TS1_5_300_s = "Net time* spent moving (s)",
    time_moving_centered_TS2_300_s   = "Net time* spent moving (s)"
  )
  
  if (v %in% names(map)) return(unname(map[[v]]))
  
  if (grepl("^time_spent_zone_", v)) {
    z <- gsub("^time_spent_zone_", "", v)
    z <- gsub("_s$", "", z)
    z <- gsub("_", " ", z)
    return(paste0("Time spent in zone: ", str_to_title(z), " (s)"))
  }
  
  str_to_title(gsub("_", " ", v))
}

is_time_spent_zone_indicator <- function(v) {
  grepl("^time_spent_zone_.*_s$", v)
}

is_time_moving_indicator <- function(v) {
  grepl("^time_moving_TS", v)
}

is_centered_flow_indicator <- function(v) {
  identical(v, "flow_time_centered_300_s")
}

is_centered_moving_indicator <- function(v) {
  grepl("^time_moving_centered_TS", v)
}

is_centered_indicator <- function(v) {
  is_centered_flow_indicator(v) || is_centered_moving_indicator(v)
}

wrap_caption <- function(x, width = 115) {
  paste(strwrap(x, width = width), collapse = "\n")
}

normalize_stats_branch_name <- function(x) {
  x <- as.character(x)[1]
  dplyr::case_when(
    x == "All_groups" ~ "all_groups",
    x == "Exe_choice" ~ "exe_choice",
    TRUE ~ tolower(trimws(x))
  )
}

normalize_stats_subbranch_name <- function(x) {
  x <- as.character(x)[1]
  dplyr::case_when(
    x == "All_timepoints" ~ "all_timepoints",
    x == "Timepoints_23"  ~ "timepoints23",
    x == "Timepoint1"     ~ "timepoint1",
    TRUE ~ tolower(gsub("[^a-z0-9]+", "", trimws(x)))
  )
}

get_anova_string_from_results <- function(branch_name, subbranch_name, indicator,
                                          results_obj = NULL) {
  if (is.null(results_obj)) {
    # Priority 1: in-memory all_branch_results (set by STEP 5 in same session)
    if (exists("all_branch_results", envir = .GlobalEnv, inherits = FALSE)) {
      results_obj <- get("all_branch_results", envir = .GlobalEnv)
    } else if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
      # R6 priority 2: canonical stats cache from a prior STEP 5 run
      .pdd <- get("PIPELINE_DATA_DIR", envir = .GlobalEnv)
      .cached <- file.path(.pdd, "dataset", "stats_cache", "all_branch_results.rds")
      if (file.exists(.cached)) {
        results_obj <- tryCatch(readRDS(.cached), error = function(e) NULL)
      }
    }
    if (is.null(results_obj)) return(NA_character_)
  }
  
  if (!is.list(results_obj)) return(NA_character_)
  
  stats_branch    <- normalize_stats_branch_name(branch_name)
  stats_subbranch <- normalize_stats_subbranch_name(subbranch_name)
  
  if (!stats_branch %in% names(results_obj)) return(NA_character_)
  if (!stats_subbranch %in% names(results_obj[[stats_branch]])) return(NA_character_)
  
  x <- results_obj[[stats_branch]][[stats_subbranch]][["anova_string"]]
  
  if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) return(NA_character_)
  if (!all(c("indicator", "anova_string") %in% names(x))) return(NA_character_)
  
  hit <- x[x$indicator == indicator, , drop = FALSE]
  if (nrow(hit) == 0) return(NA_character_)
  
  out <- as.character(hit$anova_string[1])
  out <- trimws(out)
  
  if (!nzchar(out) || is.na(out)) return(NA_character_)
  out
}

compose_scatter_caption <- function(base_caption, anova_string) {
  if (is.null(anova_string) || is.na(anova_string) || !nzchar(trimws(anova_string))) {
    return(base_caption)
  }
  paste0(base_caption, "\n\nANOVA: ", anova_string)
}

fmt_sig3 <- function(x) {
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  out[ok] <- format(signif(x[ok], 3), scientific = FALSE, trim = TRUE)
  out[!ok] <- "NA"
  out
}

make_zone_wide <- function(df_long) {
  id_cols <- c("fish_id", "trial_id", "timepoint",
               "Group", "Treatment", "trial_date")
  
  if (nrow(df_long) == 0) {
    return(tibble(
      fish_id    = character(),
      trial_id   = numeric(),
      timepoint  = numeric(),
      Group      = character(),
      Treatment  = character(),
      trial_date = as.Date(character())
    ))
  }
  
  value_tbl <- df_long %>%
    transmute(
      fish_id, trial_id, timepoint, Group, Treatment, trial_date,
      var_name = paste0("time_spent_zone_", zone_clean, "_s"),
      value    = time_spent_s
    ) %>%
    tidyr::pivot_wider(
      names_from  = var_name,
      values_from = value,
      values_fill = list(value = 0)
    )
  
  sd_tbl <- df_long %>%
    transmute(
      fish_id, trial_id, timepoint, Group, Treatment, trial_date,
      var_name = paste0("time_spent_zone_", zone_clean, "_s_SD"),
      value    = time_spent_s_SD
    ) %>%
    tidyr::pivot_wider(
      names_from  = var_name,
      values_from = value,
      values_fill = list(value = 0)
    )
  
  sem_tbl <- df_long %>%
    transmute(
      fish_id, trial_id, timepoint, Group, Treatment, trial_date,
      var_name = paste0("time_spent_zone_", zone_clean, "_s_SEM"),
      value    = time_spent_s_SEM
    ) %>%
    tidyr::pivot_wider(
      names_from  = var_name,
      values_from = value,
      values_fill = list(value = 0)
    )
  
  value_tbl %>%
    full_join(sd_tbl,  by = id_cols) %>%
    full_join(sem_tbl, by = id_cols)
}

# -----------------------------------------------------------------------------
# 3) PREP BASE DATA
# -----------------------------------------------------------------------------
df0 <- master_fish_by_frame %>%
  mutate(
    timepoint              = safe_num(timepoint),
    frame                  = safe_num(frame),
    time                   = safe_num(time),
    x_interp               = safe_num(x_interp),
    y_interp               = safe_num(y_interp),
    speed_cm_s             = safe_num(speed_cm_s),
    heading_rad            = safe_num(heading_rad),
    turn_angle_rad         = safe_num(turn_angle_rad),
    accel_cm_s2            = safe_num(accel_cm_s2),
    speed_BL_s             = safe_num(speed_BL_s),
    length_unit_px_per_cm  = safe_num(length_unit_px_per_cm),
    body_length_cm         = safe_num(body_length_cm),
    zone_clean             = sanitize_zone_name(zone),
    Group                  = normalize_group_label(Group),
    Treatment              = normalize_treatment_label(Treatment)
  ) %>%
  filter(timepoint %in% (if (.keep_tp1_active) c(1, 2, 3) else c(2, 3))) %>%
  arrange(trial_id, timepoint, fish_id, frame)

# -----------------------------------------------------------------------------
# 4) FISH-LEVEL DISTANCE + SPEED + ANGLE + ACCEL SUMMARY
# -----------------------------------------------------------------------------
fish_summary_core <- df0 %>%
  group_by(fish_id, trial_id, timepoint, Group, Treatment, trial_date) %>%
  arrange(frame, .by_group = TRUE) %>%
  mutate(
    dx_px   = x_interp - lag(x_interp),
    dy_px   = y_interp - lag(y_interp),
    step_px = sqrt(dx_px^2 + dy_px^2),
    
    step_cm = ifelse(
      is.finite(step_px) &
        is.finite(length_unit_px_per_cm) &
        length_unit_px_per_cm > 0,
      step_px / length_unit_px_per_cm,
      NA_real_
    ),
    
    step_m = ifelse(is.finite(step_cm), step_cm / 100, NA_real_),
    
    step_BL = ifelse(
      is.finite(step_cm) &
        is.finite(body_length_cm) &
        body_length_cm > 0,
      step_cm / body_length_cm,
      NA_real_
    ),
    
    prev_zone_real = !is.na(lag(zone_clean)) & lag(zone_clean) != "missing",
    curr_zone_real = !is.na(zone_clean) & zone_clean != "missing",
    
    zone_changed = dplyr::case_when(
      prev_zone_real & curr_zone_real & lag(zone_clean) != zone_clean ~ 1L,
      TRUE ~ 0L
    )
  ) %>%
  summarise(
    n_speed_BL_s_valid = sum(is.finite(speed_BL_s)),
    
    distance_m     = sum(step_m, na.rm = TRUE),
    distance_m_SD  = safe_sd(step_m),
    distance_m_SEM = safe_sem(step_m),
    
    distance_BL     = sum(step_BL, na.rm = TRUE),
    distance_BL_SD  = safe_sd(step_BL),
    distance_BL_SEM = safe_sem(step_BL),
    
    speed_cm_s_mean      = safe_mean(speed_cm_s),
    speed_cm_s_mean_SD   = safe_sd(speed_cm_s),
    speed_cm_s_mean_SEM  = safe_sem(speed_cm_s),
    
    speed_cm_s_median     = safe_median(speed_cm_s),
    speed_cm_s_median_SD  = safe_sd(speed_cm_s),
    speed_cm_s_median_SEM = safe_sem(speed_cm_s),
    
    speed_cm_s_max     = safe_max(speed_cm_s),
    speed_cm_s_max_SD  = safe_sd(speed_cm_s),
    speed_cm_s_max_SEM = safe_sem(speed_cm_s),
    
    speed_BL_s_mean      = safe_mean(speed_BL_s),
    speed_BL_s_mean_SD   = safe_sd(speed_BL_s),
    speed_BL_s_mean_SEM  = safe_sem(speed_BL_s),
    
    speed_BL_s_median     = safe_median(speed_BL_s),
    speed_BL_s_median_SD  = safe_sd(speed_BL_s),
    speed_BL_s_median_SEM = safe_sem(speed_BL_s),
    
    speed_BL_s_max     = safe_max(speed_BL_s),
    speed_BL_s_max_SD  = safe_sd(speed_BL_s),
    speed_BL_s_max_SEM = safe_sem(speed_BL_s),
    
    accel_cm_s2_mean      = safe_mean(accel_cm_s2),
    accel_cm_s2_mean_SD   = safe_sd(accel_cm_s2),
    accel_cm_s2_mean_SEM  = safe_sem(accel_cm_s2),
    
    accel_cm_s2_median     = safe_median(accel_cm_s2),
    accel_cm_s2_median_SD  = safe_sd(accel_cm_s2),
    accel_cm_s2_median_SEM = safe_sem(accel_cm_s2),
    
    accel_cm_s2_max     = safe_max(accel_cm_s2),
    accel_cm_s2_max_SD  = safe_sd(accel_cm_s2),
    accel_cm_s2_max_SEM = safe_sem(accel_cm_s2),
    
    n_zone_transitions     = sum(zone_changed, na.rm = TRUE),
    n_zone_transitions_SD  = safe_sd(zone_changed),
    n_zone_transitions_SEM = safe_sem(zone_changed),
    
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 5) FISH-LEVEL TIME MOVING SUMMARY
# -----------------------------------------------------------------------------
fish_summary_move <- df0 %>%
  group_by(fish_id, trial_id, timepoint, Group, Treatment, trial_date) %>%
  arrange(frame, .by_group = TRUE) %>%
  mutate(
    dt = calc_dt(time),
    
    move_TS0.5 = is_moving_flag(locomotor_status_TS0.5),
    move_TS1   = is_moving_flag(locomotor_status_TS1),
    move_TS1.5 = is_moving_flag(locomotor_status_TS1.5),
    move_TS2   = is_moving_flag(locomotor_status_TS2),
    
    contrib_TS0.5 = ifelse(move_TS0.5, dt, 0),
    contrib_TS1   = ifelse(move_TS1,   dt, 0),
    contrib_TS1.5 = ifelse(move_TS1.5, dt, 0),
    contrib_TS2   = ifelse(move_TS2,   dt, 0)
  ) %>%
  summarise(
    time_moving_TS0.5_s     = sum(contrib_TS0.5, na.rm = TRUE),
    time_moving_TS0.5_s_SD  = safe_sd(contrib_TS0.5),
    time_moving_TS0.5_s_SEM = safe_sem(contrib_TS0.5),
    
    time_moving_TS1_s     = sum(contrib_TS1, na.rm = TRUE),
    time_moving_TS1_s_SD  = safe_sd(contrib_TS1),
    time_moving_TS1_s_SEM = safe_sem(contrib_TS1),
    
    time_moving_TS1.5_s     = sum(contrib_TS1.5, na.rm = TRUE),
    time_moving_TS1.5_s_SD  = safe_sd(contrib_TS1.5),
    time_moving_TS1.5_s_SEM = safe_sem(contrib_TS1.5),
    
    time_moving_TS2_s     = sum(contrib_TS2, na.rm = TRUE),
    time_moving_TS2_s_SD  = safe_sd(contrib_TS2),
    time_moving_TS2_s_SEM = safe_sem(contrib_TS2),
    
    .groups = "drop"
  )

# -----------------------------------------------------------------------------
# 6) FISH-LEVEL TIME SPENT PER ZONE SUMMARY
# -----------------------------------------------------------------------------
fish_zone_long <- df0 %>%
  group_by(fish_id, trial_id, timepoint, Group, Treatment, trial_date) %>%
  arrange(frame, .by_group = TRUE) %>%
  mutate(
    dt = calc_dt(time)
  ) %>%
  group_by(fish_id, trial_id, timepoint, Group, Treatment, trial_date, zone_clean) %>%
  summarise(
    time_spent_s     = sum(dt, na.rm = TRUE),
    time_spent_s_SD  = safe_sd(dt),
    time_spent_s_SEM = safe_sem(dt),
    .groups = "drop"
  )

fish_zone_summary <- make_zone_wide(fish_zone_long)

# -----------------------------------------------------------------------------
# 7) FINAL FISH SUMMARY
# -----------------------------------------------------------------------------
fish_activity_summary <- fish_summary_core %>%
  left_join(
    fish_summary_move,
    by = c("fish_id", "trial_id", "timepoint", "Group", "Treatment", "trial_date")
  ) %>%
  left_join(
    fish_zone_summary,
    by = c("fish_id", "trial_id", "timepoint", "Group", "Treatment", "trial_date")
  ) %>%
  mutate(
    Group     = normalize_group_label(Group),
    Treatment = normalize_treatment_label(Treatment),
    timepoint = as.integer(timepoint)
  )

zone_time_cols <- names(fish_activity_summary)[grepl("^time_spent_zone_.*_s$", names(fish_activity_summary))]

flow_zone_cols <- zone_time_cols[grepl("(^time_spent_zone_flow_s$)|(^time_spent_zone_.*flow.*_s$)", zone_time_cols)]

if (length(flow_zone_cols) == 0) {
  stop(
    "No flow zone time column found in fish_activity_summary. Available zone-time columns: ",
    paste(zone_time_cols, collapse = ", "),
    call. = FALSE
  )
}

if (length(flow_zone_cols) > 1) {
  message(
    "Multiple possible flow zone columns found: ",
    paste(flow_zone_cols, collapse = ", "),
    ". Using the first one: ", flow_zone_cols[1]
  )
}

FLOW_ZONE_COL <- flow_zone_cols[1]

fish_activity_summary <- fish_activity_summary %>%
  mutate(
    flow_time_centered_300_s = .data[[FLOW_ZONE_COL]] - 300,
    time_moving_centered_TS0_5_300_s = time_moving_TS0.5_s - 300,
    time_moving_centered_TS1_300_s   = time_moving_TS1_s   - 300,
    time_moving_centered_TS1_5_300_s = time_moving_TS1.5_s - 300,
    time_moving_centered_TS2_300_s   = time_moving_TS2_s   - 300
  )

unknown_groups <- setdiff(unique(fish_activity_summary$Group), GROUP_ORDER)
unknown_treatments <- setdiff(unique(fish_activity_summary$Treatment), TREATMENT_ORDER)

if (length(unknown_groups) > 0) {
  warning("Unexpected Group labels found: ", paste(unknown_groups, collapse = ", "), call. = FALSE)
}
if (length(unknown_treatments) > 0) {
  warning("Unexpected Treatment labels found: ", paste(unknown_treatments, collapse = ", "), call. = FALSE)
}

fish_activity_summary <- fish_activity_summary %>%
  filter(Group %in% GROUP_ORDER, Treatment %in% TREATMENT_ORDER) %>%
  mutate(
    Group = factor(Group, levels = GROUP_ORDER),
    Treatment = factor(Treatment, levels = TREATMENT_ORDER),
    timepoint = factor(timepoint, levels = c(1, 2, 3))
  )

assign("fish_activity_summary", fish_activity_summary, envir = .GlobalEnv)
readr::write_csv(fish_activity_summary, SUMMARY_CSV_PATH)

# ---- R1: also write to canonical dataset/ folder ----------------------------
# Stable, fixed-name output that STEP 5 and any future re-analysis script can
# read without needing to discover the latest STEP4_output_<TS>/ folder. This
# eliminates the "most-recent-by-mtime" failure mode where a partial STEP 4
# leaves a tiny fish_activity_summary.csv that STEP 5 then loads.
.canonical_dataset_dir <- file.path(PIPELINE_DATA_DIR, "dataset")
dir.create(.canonical_dataset_dir, recursive = TRUE, showWarnings = FALSE)
.canonical_summary_path <- file.path(.canonical_dataset_dir,
                                      "fish_activity_summary.csv")
readr::write_csv(fish_activity_summary, .canonical_summary_path)
message("R1: also wrote canonical copy: ", .canonical_summary_path)
rm(.canonical_dataset_dir, .canonical_summary_path)

fish_activity_summary <- fish_activity_summary %>%
  filter(n_speed_BL_s_valid > 0) %>%
  select(-n_speed_BL_s_valid)

dropped_n <- fish_summary_core %>%
  filter(n_speed_BL_s_valid == 0) %>%
  nrow()

if (dropped_n > 0) {
  message("Dropped ", dropped_n,
          " fish summaries with no valid speed_BL_s values before plotting.")
}

message("Groups in fish_activity_summary: ", paste(unique(as.character(fish_activity_summary$Group)), collapse = ", "))
message("Treatments in fish_activity_summary: ", paste(unique(as.character(fish_activity_summary$Treatment)), collapse = ", "))
message("Timepoints in fish_activity_summary: ", paste(unique(as.character(fish_activity_summary$timepoint)), collapse = ", "))

# -----------------------------------------------------------------------------
# 8) PLOTTING HELPERS
# -----------------------------------------------------------------------------
filter_outliers_percentile <- function(df_points, y_col, dataset_name, series_tag) {
  y <- df_points[[y_col]]
  qs <- safe_quantile(y, c(OUT_P_LO, OUT_P_HI))
  lo <- qs[1]
  hi <- qs[2]
  
  if (!is.finite(lo) || !is.finite(hi) || lo >= hi) {
    return(list(
      df_kept        = df_points,
      df_removed_log = df_points[0, , drop = FALSE],
      bounds         = c(lo = lo, hi = hi)
    ))
  }
  
  keep    <- is.finite(y) & y >= lo & y <= hi
  removed <- df_points[!keep & is.finite(y), , drop = FALSE]
  
  removed_log <- if (nrow(removed) > 0) {
    removed %>%
      mutate(
        dataset      = dataset_name,
        indicator    = y_col,
        value        = .data[[y_col]],
        outlier_rule = paste0("percentile_", OUT_P_LO, "_", OUT_P_HI),
        p_low        = lo,
        p_high       = hi,
        series       = series_tag,
        removed_at   = Sys.time()
      )
  } else {
    removed
  }
  
  list(
    df_kept        = df_points[keep | !is.finite(y), , drop = FALSE],
    df_removed_log = removed_log,
    bounds         = c(lo = lo, hi = hi)
  )
}

make_outlier_caption <- function(raw_df, removed_log_df, y_col) {
  base_txt <- paste0(
    "N = 12 for each treatment × group × timepoint combination, ",
    "except for trial 10, N = 11"
  )
  
  if (is.null(removed_log_df) || nrow(removed_log_df) == 0) {
    return(wrap_caption(base_txt))
  }
  
  outlier_counts <- removed_log_df %>%
    count(timepoint, Treatment, Group, name = "n_outliers") %>%
    mutate(
      timepoint = as.character(timepoint),
      Treatment = as.character(Treatment),
      Group = as.character(Group)
    )
  
  raw_counts <- raw_df %>%
    filter(is.finite(.data[[y_col]])) %>%
    count(timepoint, Treatment, Group, name = "n_before") %>%
    mutate(
      timepoint = as.character(timepoint),
      Treatment = as.character(Treatment),
      Group = as.character(Group)
    )
  
  cap_tbl <- raw_counts %>%
    left_join(outlier_counts, by = c("timepoint", "Treatment", "Group")) %>%
    mutate(
      n_outliers = dplyr::coalesce(n_outliers, 0L),
      n_after = n_before - n_outliers
    ) %>%
    filter(n_outliers > 0) %>%
    mutate(
      txt = paste0(
        "TP", timepoint, ": ",
        Treatment, " × ", Group,
        " N = ", n_after
      )
    )
  
  if (nrow(cap_tbl) == 0) {
    return(wrap_caption(base_txt))
  }
  
  extra_txt <- paste(cap_tbl$txt, collapse = ", ")
  wrap_caption(paste0(base_txt, "; after outlier removal: ", extra_txt))
}

prepare_plot_data <- function(df_points, y_col, dataset_name, series_tag) {
  filt <- filter_outliers_percentile(
    df_points    = df_points,
    y_col        = y_col,
    dataset_name = dataset_name,
    series_tag   = series_tag
  )
  
  kept <- filt$df_kept %>%
    filter(is.finite(.data[[y_col]]))
  
  if (nrow(kept) == 0) {
    return(list(
      kept               = kept,
      removed_log        = filt$df_removed_log,
      mean_sem           = tibble(),
      y_limits           = c(NA_real_, NA_real_),
      caption_text       = NA_character_,
      net_caption_text   = NA_character_
    ))
  }
  
  centered_var    <- is_centered_indicator(y_col)
  zone_time_var   <- is_time_spent_zone_indicator(y_col)
  moving_time_var <- is_time_moving_indicator(y_col)
  
  yy   <- kept[[y_col]]
  ymin <- min(yy, na.rm = TRUE)
  ymax <- max(yy, na.rm = TRUE)
  
  if (centered_var) {
    lower <- -300
    upper <- 300
  } else if (zone_time_var || moving_time_var) {
    lower <- min(0, ymin)
    upper <- 600
  } else {
    yrng  <- ymax - ymin
    lower <- ymin - 0.05 * (yrng + 1e-9)
    upper <- ymax + 0.05 * (yrng + 1e-9)
  }
  
  mean_sem <- kept %>%
    group_by(timepoint, Treatment, Group) %>%
    summarise(
      n    = n(),
      mean = mean(.data[[y_col]], na.rm = TRUE),
      sd   = stats::sd(.data[[y_col]], na.rm = TRUE),
      sem  = ifelse(n > 0, sd / sqrt(n), NA_real_),
      .groups = "drop"
    )
  
  caption_text <- make_outlier_caption(
    raw_df = df_points,
    removed_log_df = filt$df_removed_log,
    y_col = y_col
  )
  
  net_caption_text <- if (centered_var) {
    "(*) net with respect to 300 s (5 min) i.e. half-time of each timepoint."
  } else {
    NA_character_
  }
  
  list(
    kept             = kept,
    removed_log      = filt$df_removed_log,
    mean_sem         = mean_sem,
    y_limits         = c(lower, upper),
    caption_text     = caption_text,
    net_caption_text = net_caption_text
  )
}

# -----------------------------------------------------------------------------
# 8B) TIMEPOINT-TREND GRAPH HELPERS
# -----------------------------------------------------------------------------
make_combo_offsets <- function(df) {
  combos <- df %>%
    distinct(Group, Treatment) %>%
    mutate(
      Group = as.character(Group),
      Treatment = as.character(Treatment)
    ) %>%
    arrange(
      factor(Group, levels = GROUP_ORDER),
      factor(Treatment, levels = TREATMENT_ORDER)
    )
  
  n <- nrow(combos)
  
  if (n == 0) {
    combos$offset <- numeric(0)
  } else if (n == 1) {
    combos$offset <- 0
  } else {
    combos$offset <- seq(-0.24, 0.24, length.out = n)
  }
  
  combos
}

make_timepoint_mean_sem <- function(df_points, y_col) {
  df_points %>%
    filter(is.finite(.data[[y_col]])) %>%
    mutate(
      timepoint_num = as.integer(as.character(timepoint)),
      Group = as.character(Group),
      Treatment = as.character(Treatment)
    ) %>%
    group_by(timepoint_num, Group, Treatment) %>%
    summarise(
      n = n(),
      mean = mean(.data[[y_col]], na.rm = TRUE),
      sd   = stats::sd(.data[[y_col]], na.rm = TRUE),
      sem  = ifelse(n > 0, sd / sqrt(n), NA_real_),
      .groups = "drop"
    ) %>%
    left_join(make_combo_offsets(df_points), by = c("Group", "Treatment")) %>%
    mutate(
      x_plot = timepoint_num + offset,
      timepoint = factor(timepoint_num, levels = c(1, 2, 3)),
      Group = factor(Group, levels = GROUP_ORDER),
      Treatment = factor(Treatment, levels = TREATMENT_ORDER)
    )
}

make_corr_caption <- function(df_points, y_col) {
  corr_tbl <- df_points %>%
    filter(is.finite(.data[[y_col]])) %>%
    mutate(
      timepoint_num = as.integer(as.character(timepoint)),
      Group = as.character(Group),
      Treatment = as.character(Treatment)
    ) %>%
    group_by(Treatment, Group) %>%
    group_modify(~{
      dd <- .x %>%
        filter(is.finite(timepoint_num), is.finite(.data[[y_col]]))
      
      if (nrow(dd) < 3 || dplyr::n_distinct(dd$timepoint_num) < 2) {
        return(tibble(
          estimate = NA_real_,
          p_value  = NA_real_,
          n        = nrow(dd)
        ))
      }
      
      ct <- suppressWarnings(stats::cor.test(
        x = dd$timepoint_num,
        y = dd[[y_col]],
        method = "pearson"
      ))
      
      tibble(
        estimate = unname(ct$estimate),
        p_value  = ct$p.value,
        n        = nrow(dd)
      )
    }) %>%
    ungroup() %>%
    filter(is.finite(estimate), is.finite(p_value), p_value < 0.05) %>%
    mutate(
      line = paste0(
        Treatment, " × ", Group,
        ": r(", pmax(n - 2, 1), ") = ", fmt_sig3(estimate),
        ", p = ", fmt_sig3(p_value)
      )
    )
  
  if (nrow(corr_tbl) == 0) {
    return("No significant Pearson correlations between timepoint and indicator within treatment × group.")
  }
  
  paste(corr_tbl$line, collapse = "\n")
}

plot_indicator_timepoint_trend <- function(prepped, y_col, dataset_title, series_tag) {
  dfp <- prepped$kept
  if (nrow(dfp) == 0) return(NULL)
  
  centered_var <- is_centered_indicator(y_col)
  if (centered_var) return(NULL)
  
  # skip trend plots if fewer than 2 distinct timepoints survive
  if (!has_at_least_two_timepoints(dfp)) {
    message("[", series_tag, "] skipping trend graph for ", y_col, " (< 2 distinct timepoints)")
    return(NULL)
  }
  
  ms <- make_timepoint_mean_sem(dfp, y_col)
  if (nrow(ms) == 0) return(NULL)
  
  # identify which Treatment × Group combos have >= 2 timepoints
  ms_line <- ms %>%
    mutate(
      Group = as.character(Group),
      Treatment = as.character(Treatment)
    ) %>%
    group_by(Group, Treatment) %>%
    filter(sum(is.finite(mean)) >= 2, dplyr::n_distinct(timepoint) >= 2) %>%
    ungroup() %>%
    mutate(
      Group = factor(Group, levels = GROUP_ORDER),
      Treatment = factor(Treatment, levels = TREATMENT_ORDER)
    )
  
  y_lab <- pretty_indicator_name(y_col)
  ttl   <- paste0(dataset_title, " | ", y_lab, " | Timepoint trends | ", series_tag)
  
  zone_time_var   <- is_time_spent_zone_indicator(y_col)
  moving_time_var <- is_time_moving_indicator(y_col)
  
  yy   <- c(ms$mean - ms$sem, ms$mean + ms$sem)
  yy   <- yy[is.finite(yy)]
  ymin <- min(yy, na.rm = TRUE)
  ymax <- max(yy, na.rm = TRUE)
  
  if (zone_time_var || moving_time_var) {
    lower <- min(0, ymin)
    upper <- 600
  } else {
    yrng  <- ymax - ymin
    lower <- ymin - 0.05 * (yrng + 1e-9)
    upper <- ymax + 0.05 * (yrng + 1e-9)
  }
  
  corr_caption <- make_corr_caption(dfp, y_col)
  
  p <- ggplot(
    ms,
    aes(
      x = x_plot,
      y = mean,
      color = Group,
      shape = Treatment,
      group = interaction(Group, Treatment),
      linetype = Treatment
    )
  ) +
    # draw lines only for combos with >= 2 timepoints
    geom_line(
      data = ms_line,
      linewidth = TREND_LINE_SIZE
    ) +
    geom_errorbar(
      aes(ymin = mean - sem, ymax = mean + sem),
      width = TREND_ERRORBAR_WIDTH,
      linewidth = ERRORBAR_LINE_SIZE
    ) +
    geom_point(size = TREND_POINT_SIZE) +
    labs(
      title = ttl,
      x = "Interval",
      y = y_lab,
      color = "Group",
      shape = "Treatment",
      linetype = "Treatment",
      caption = corr_caption
    ) +
    scale_color_manual(
      values = GROUP_COLORS,
      breaks = GROUP_ORDER,
      drop = TRUE
    ) +
    scale_shape_manual(
      values = setNames(c(17, 16), TREATMENT_ORDER),
      breaks = TREATMENT_ORDER,
      drop = TRUE
    ) +
    scale_linetype_manual(
      values = setNames(c("dashed", "solid"), TREATMENT_ORDER),
      breaks = TREATMENT_ORDER,
      drop = TRUE
    ) +
    scale_x_continuous(
      breaks = c(1, 2, 3),
      labels = c("Interval 1", "Interval 2", "Interval 3"),
      limits = c(0.5, 3.5)
    ) +
    BASE_THEME
  
  if (zone_time_var || moving_time_var) {
    p <- p +
      geom_hline(
        yintercept = 300,
        linetype   = "dashed",
        color      = "red",
        linewidth  = 1.1,
        alpha      = 0.85
      ) +
      scale_y_continuous(
        breaks = seq(0, 600, 100),
        limits = c(lower, upper),
        expand = expansion(mult = c(0.02, 0.04))
      )
  } else {
    p <- p +
      coord_cartesian(ylim = c(lower, upper)) +
      scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.04))
      )
  }
  
  p
}

plot_indicator_scatter <- function(prepped, y_col, dataset_title, series_tag, anova_string = NA_character_) {
  dfp <- prepped$kept
  if (nrow(dfp) == 0) return(NULL)
  
  ms    <- prepped$mean_sem
  lower <- prepped$y_limits[1]
  upper <- prepped$y_limits[2]
  
  y_lab <- pretty_indicator_name(y_col)
  ttl   <- paste0(dataset_title, " | ", y_lab, " | ", series_tag)
  
  centered_flow_var   <- is_centered_flow_indicator(y_col)
  centered_moving_var <- is_centered_moving_indicator(y_col)
  zone_time_var       <- is_time_spent_zone_indicator(y_col)
  moving_time_var     <- is_time_moving_indicator(y_col)
  
  full_caption <- if (!is.na(prepped$net_caption_text)) {
    paste0(prepped$net_caption_text, "\n\n", prepped$caption_text)
  } else {
    prepped$caption_text
  }
  
  p <- ggplot(dfp, aes(x = Treatment, y = .data[[y_col]]))
  
  if (centered_flow_var) {
    p <- p +
      annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0, ymax = 300,
        fill = "lightskyblue2",
        alpha = 0.22
      )
  }
  
  if (centered_moving_var) {
    p <- p +
      annotate(
        "rect",
        xmin = -Inf, xmax = Inf,
        ymin = 0, ymax = 300,
        fill = "gold1",
        alpha = 0.20
      )
  }
  
  p <- p +
    geom_crossbar(
      data = ms,
      aes(x = Treatment, y = mean, ymin = mean, ymax = mean, group = Group),
      position    = POS_DODGE,
      inherit.aes = FALSE,
      color       = "black",
      width       = MEAN_CROSSBAR_WIDTH,
      linewidth   = MEAN_CROSSBAR_LINE_SIZE
    ) +
    geom_errorbar(
      data = ms,
      aes(x = Treatment, ymin = mean - sem, ymax = mean + sem, group = Group),
      position    = POS_DODGE,
      width       = ERRORBAR_WIDTH,
      color       = "black",
      linewidth   = ERRORBAR_LINE_SIZE,
      inherit.aes = FALSE
    ) +
    geom_point(
      aes(color = Group, shape = Treatment, group = Group),
      position = POS_JD,
      alpha    = PT_ALPHA,
      size     = PT_SIZE
    ) +
    facet_wrap(
      ~ timepoint,
      nrow = FACET_NROW,
      labeller = as_labeller(function(x) paste0("Interval ", x))
    ) +
    labs(
      title   = ttl,
      x       = NULL,
      y       = y_lab,
      color   = "Group",
      shape   = "Treatment",
      caption = full_caption
    ) +
    scale_color_manual(
      values = GROUP_COLORS,
      breaks = GROUP_ORDER,
      drop = TRUE
    ) +
    scale_shape_discrete(drop = TRUE) +
    BASE_THEME
  
  if (centered_flow_var || centered_moving_var) {
    p <- p +
      geom_hline(
        yintercept = 0,
        linetype   = "solid",
        color      = "black",
        linewidth  = 0.8
      ) +
      scale_y_continuous(
        limits = c(-300, 300),
        breaks = seq(-300, 300, by = 100),
        expand = expansion(mult = c(0, 0.02))
      )
  } else if (zone_time_var || moving_time_var) {
    p <- p +
      geom_hline(
        yintercept = 300,
        linetype   = "dashed",
        color      = "red",
        linewidth  = 1.1,
        alpha      = 0.85
      ) +
      scale_y_continuous(
        limits = c(lower, upper),
        breaks = seq(0, 600, 100),
        expand = expansion(mult = c(0.02, 0.04))
      )
  } else {
    p <- p +
      coord_cartesian(ylim = c(lower, upper)) +
      scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.04))
      )
  }
  
  # If no ANOVA string is available, return the original scatter plot
  if (is.na(anova_string) || !nzchar(trimws(anova_string))) {
    return(p)
  }
  
  anova_caption_text <- wrap_caption(paste0("ANOVA: ", anova_string), width = 120)
  
  anova_grob <- cowplot::ggdraw() +
    cowplot::draw_label(
      anova_caption_text,
      x = 0, y = 1,
      hjust = 0, vjust = 1,
      size = ANOVA_CAPTION_TEXT_SIZE,
      lineheight = ANOVA_CAPTION_LINEHEIGHT
    ) +
    theme(plot.margin = margin(t = 0, r = 10, b = 2, l = 10))
  
  cowplot::plot_grid(
    p,
    anova_grob,
    ncol = 1,
    align = "v",
    axis = "lr",
    rel_heights = c(1, 0.10)
  )
}

ppt_add_plot <- function(ppt, plot_obj, slide_title = NULL) {
  ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
  if (!is.null(slide_title)) {
    ppt <- ph_with(ppt, slide_title, location = ph_location_type(type = "title"))
  }
  ppt <- ph_with(
    ppt,
    rvg::dml(ggobj = plot_obj),
    location = ph_location_type(type = "body")
  )
  ppt
}

# -----------------------------------------------------------------------------
# 9) INDICATOR SELECTION
# -----------------------------------------------------------------------------
indicator_vars <- names(fish_activity_summary)[vapply(fish_activity_summary, is.numeric, logical(1))]
indicator_vars <- setdiff(indicator_vars, "trial_id")
indicator_vars <- indicator_vars[!grepl("(_SD|_SEM)$", indicator_vars)]
indicator_vars <- setdiff(indicator_vars, c("n_speed_BL_s_valid"))
indicator_vars <- unique(indicator_vars)

# Plotting principle guard: every y-axis variable must come from a raw
# (untransformed) column in fish_activity_summary — never from model
# predictions, EMMs, or Box-Cox-transformed values.
# Inference lives in Word reports and CSV tables; plots show raw data only.
.bad_plotvars <- indicator_vars[!indicator_vars %in% names(fish_activity_summary)]
if (length(.bad_plotvars) > 0) {
  stop(
    "Raw-data assertion failed: the following plot variables are not present ",
    "in fish_activity_summary (they may be model-derived or transformed):\n  ",
    paste(.bad_plotvars, collapse = ", "),
    "\nAll plot y-axes must come from raw fish_activity_summary columns.",
    call. = FALSE
  )
}
rm(.bad_plotvars)

pretty_slide_title <- function(v, trend = FALSE) {
  slide_map <- c(
    flow_time_centered_300_s         = "Net time spent in flow | centered at 300 s",
    time_moving_centered_TS0_5_300_s = "Net time spent moving | centered at 300 s | TS0.5",
    time_moving_centered_TS1_300_s   = "Net time spent moving | centered at 300 s | TS1",
    time_moving_centered_TS1_5_300_s = "Net time spent moving | centered at 300 s | TS1.5",
    time_moving_centered_TS2_300_s   = "Net time spent moving | centered at 300 s | TS2"
  )
  
  out <- if (v %in% names(slide_map)) unname(slide_map[[v]]) else pretty_indicator_name(v)
  
  if (isTRUE(trend)) {
    out <- paste0(out, " | Interval trend")
  }
  
  out
}

# -----------------------------------------------------------------------------
# 10) RUN BRANCHES × SUBBRANCHES
# -----------------------------------------------------------------------------
run_activity_branch_subbranch <- function(df_summary,
                                          branch_name,
                                          group_filter,
                                          subbranch_name,
                                          timepoint_filter,
                                          main_export_dir) {
  df_use <- df_summary

  if (!is.null(group_filter)) {
    df_use <- df_use %>% filter(Group %in% group_filter)
  }

  if (!is.null(timepoint_filter)) {
    df_use <- df_use %>% filter(as.integer(as.character(timepoint)) %in% timepoint_filter)
  }

  df_use <- df_use %>%
    mutate(
      Group = normalize_group_label(Group),
      Treatment = normalize_treatment_label(Treatment),
      timepoint = as.integer(as.character(timepoint))
    ) %>%
    filter(
      Group %in% GROUP_ORDER,
      Treatment %in% TREATMENT_ORDER,
      timepoint %in% c(1, 2, 3)
    ) %>%
    mutate(
      Group = factor(Group, levels = GROUP_ORDER),
      Treatment = factor(Treatment, levels = TREATMENT_ORDER),
      timepoint = factor(timepoint, levels = c(1, 2, 3))
    ) %>%
    droplevels()

  if (nrow(df_use) == 0) {
    warning(
      "No rows left after filtering for branch/subbranch: ",
      branch_name, " / ", subbranch_name,
      call. = FALSE
    )
    return(invisible(list(
      ppt_path = NA_character_,
      outlier_csv_path = NA_character_
    )))
  }

  message(
    "[", branch_name, " / ", subbranch_name, "] rows = ", nrow(df_use),
    " | groups = ", paste(unique(as.character(df_use$Group)), collapse = ", "),
    " | treatments = ", paste(unique(as.character(df_use$Treatment)), collapse = ", "),
    " | timepoints = ", paste(unique(as.character(df_use$timepoint)), collapse = ", ")
  )

  branch_dir <- file.path(main_export_dir, branch_name)
  subbranch_dir <- file.path(branch_dir, subbranch_name)

  if (!dir.exists(branch_dir)) dir.create(branch_dir, recursive = TRUE)
  if (!dir.exists(subbranch_dir)) dir.create(subbranch_dir, recursive = TRUE)

  png_dir <- file.path(subbranch_dir, "interval_level_png")
  if (!dir.exists(png_dir)) dir.create(png_dir, recursive = TRUE)

  ppt_path <- file.path(
    subbranch_dir,
    paste0("activity_indicators_graphs_", branch_name, "_", subbranch_name, "_", TS_NOW, ".pptx")
  )

  outlier_csv_path <- file.path(
    subbranch_dir,
    paste0("activity_outlier_log_", branch_name, "_", subbranch_name, ".csv")
  )

  series_tag <- paste(branch_name, subbranch_name, sep = " | ")

  subbranch_has_trend <- (!is.null(timepoint_filter) && length(timepoint_filter) >= 2) ||
    is.null(timepoint_filter)

  # =====================================================================
  # PHASE A: PARALLEL PLOT GENERATION
  # =====================================================================
  # Build task combinations: for each indicator, generate both scatter plot
  # and (optionally) trend plot. Each task returns a list with plot info.
  #
  # Task structure: list(indicator, plot_type) where plot_type is "scatter" or "trend"
  # This allows us to parallelize the expensive ggplot rendering + RVG serialization.

  task_list <- list()
  for (v in indicator_vars) {
    task_list[[length(task_list) + 1]] <- list(indicator = v, plot_type = "scatter")
    if (!is_centered_indicator(v) && subbranch_has_trend) {
      task_list[[length(task_list) + 1]] <- list(indicator = v, plot_type = "trend")
    }
  }

  # OPT-4: Use the n_workers_step4 set at the top level of the script
  effective_workers <- if (n_workers_step4 > 0L) n_workers_step4 else 1L
  message("[OPT-4] Parallelizing ", length(task_list), " plots (scatter + trend) across ",
          effective_workers, " worker(s)...")

  # Run parallel plot generation
  plot_list <- future.apply::future_lapply(
    task_list,
    function(task_spec) {
      v <- task_spec$indicator
      plot_type <- task_spec$plot_type

      # Prepare data for this indicator
      prepped <- prepare_plot_data(
        df_points    = df_use,
        y_col        = v,
        dataset_name = "fish_activity_summary",
        series_tag   = series_tag
      )

      # Skip if no data
      if (nrow(prepped$kept) == 0) {
        return(NULL)
      }

      # Get ANOVA string for annotation
      anova_string_this <- get_anova_string_from_results(
        branch_name    = branch_name,
        subbranch_name = subbranch_name,
        indicator      = v
      )

      # Generate the appropriate plot
      if (plot_type == "scatter") {
        p <- plot_indicator_scatter(
          prepped       = prepped,
          y_col         = v,
          dataset_title = "Fish activity summary",
          series_tag    = series_tag,
          anova_string  = anova_string_this
        )
      } else if (plot_type == "trend") {
        p <- plot_indicator_timepoint_trend(
          prepped       = prepped,
          y_col         = v,
          dataset_title = "Fish activity summary",
          series_tag    = series_tag
        )
      } else {
        return(NULL)
      }

      if (is.null(p)) {
        return(NULL)
      }

      # Serialize to DML (officer-compatible format)
      dml_obj <- rvg::dml(ggobj = p)

      # Generate slide title
      slide_title <- paste0(
        branch_name, " | ", subbranch_name, " | ",
        pretty_slide_title(v, trend = (plot_type == "trend"))
      )

      # Return metadata + serialized plot
      list(
        indicator = v,
        plot_type = plot_type,
        dml = dml_obj,
        gg  = p,
        slide_title = slide_title,
        removed_log = prepped$removed_log
      )
    },
    future.packages = c("ggplot2", "rvg", "dplyr", "tidyr", "officer", "cowplot", "stringr", "stats"),
    future.globals = list(
      df_use = df_use,
      branch_name = branch_name,
      subbranch_name = subbranch_name,
      series_tag = series_tag,
      is_centered_indicator = is_centered_indicator,
      is_centered_flow_indicator = is_centered_flow_indicator,
      is_centered_moving_indicator = is_centered_moving_indicator,
      is_time_spent_zone_indicator = is_time_spent_zone_indicator,
      is_time_moving_indicator = is_time_moving_indicator,
      prepare_plot_data = prepare_plot_data,
      plot_indicator_scatter = plot_indicator_scatter,
      plot_indicator_timepoint_trend = plot_indicator_timepoint_trend,
      get_anova_string_from_results = get_anova_string_from_results,
      pretty_slide_title = pretty_slide_title,
      pretty_indicator_name = pretty_indicator_name,
      pretty_var_label = pretty_var_label,
      wrap_caption = wrap_caption,
      filter_outliers_percentile = filter_outliers_percentile,
      make_outlier_caption = make_outlier_caption,
      make_combo_offsets = make_combo_offsets,
      make_timepoint_mean_sem = make_timepoint_mean_sem,
      make_corr_caption = make_corr_caption,
      safe_quantile = safe_quantile,
      has_at_least_two_timepoints = has_at_least_two_timepoints,
      fmt_sig3 = fmt_sig3,
      PT_ALPHA = PT_ALPHA,
      PT_SIZE = PT_SIZE,
      JITTER_WIDTH = JITTER_WIDTH,
      DODGE_WIDTH = DODGE_WIDTH,
      GROUP_COLORS = GROUP_COLORS,
      GROUP_ORDER = GROUP_ORDER,
      TREATMENT_ORDER = TREATMENT_ORDER,
      BASE_THEME = BASE_THEME,
      POS_JD = POS_JD,
      POS_DODGE = POS_DODGE,
      MEAN_CROSSBAR_WIDTH = MEAN_CROSSBAR_WIDTH,
      MEAN_CROSSBAR_LINE_SIZE = MEAN_CROSSBAR_LINE_SIZE,
      ERRORBAR_WIDTH = ERRORBAR_WIDTH,
      ERRORBAR_LINE_SIZE = ERRORBAR_LINE_SIZE,
      FACET_NROW = FACET_NROW,
      PANEL_SPACING_LINES = PANEL_SPACING_LINES,
      OUT_P_LO = OUT_P_LO,
      OUT_P_HI = OUT_P_HI,
      N_LABEL_SIZE = N_LABEL_SIZE,
      N_LABEL_PREFIX = N_LABEL_PREFIX,
      CAPTION_TEXT_SIZE = CAPTION_TEXT_SIZE,
      ANOVA_CAPTION_TEXT_SIZE = ANOVA_CAPTION_TEXT_SIZE,
      ANOVA_CAPTION_LINEHEIGHT = ANOVA_CAPTION_LINEHEIGHT,
      TREND_LINE_SIZE = TREND_LINE_SIZE,
      TREND_POINT_SIZE = TREND_POINT_SIZE,
      TREND_ERRORBAR_WIDTH = TREND_ERRORBAR_WIDTH,
      GOT_INDICATOR_LABELS = if (exists("GOT_INDICATOR_LABELS")) GOT_INDICATOR_LABELS else character(0)
    ),
    future.seed = TRUE
  )

  # =====================================================================
  # PHASE B: SERIAL POWERPOINT ASSEMBLY
  # =====================================================================
  # Insert all DML objects into a fresh PowerPoint deck in deterministic order.
  # This avoids contention on the mutable ppt object.

  ppt <- read_pptx()
  outlier_log <- list()

  # Filter out NULL results and process in order
  valid_plots <- Filter(Negate(is.null), plot_list)

  message("[OPT-4] Inserting ", length(valid_plots), " plots into PowerPoint...")

  for (plot_item in valid_plots) {
    # Add slide with title and content
    ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
    ppt <- ph_with(ppt, plot_item$slide_title, location = ph_location_type(type = "title"))
    ppt <- ph_with(
      ppt,
      plot_item$dml,
      location = ph_location_type(type = "body")
    )

    # Save PNG (interval-level output)
    if (!is.null(plot_item$gg)) {
      safe_indic <- gsub("[^A-Za-z0-9_]", "_", plot_item$indicator)
      png_file   <- file.path(png_dir, paste0(safe_indic, "_", plot_item$plot_type, ".png"))
      tryCatch(
        ggplot2::ggsave(png_file, plot = plot_item$gg, width = 10, height = 7,
                        dpi = 150, bg = "white"),
        error = function(e) warning("PNG save failed for ", plot_item$indicator, ": ",
                                    conditionMessage(e), call. = FALSE)
      )
    }

    # Collect outlier logs
    if (nrow(plot_item$removed_log) > 0) {
      outlier_log[[length(outlier_log) + 1]] <- plot_item$removed_log
    }
  }

  # =====================================================================
  # POST-PROCESSING
  # =====================================================================

  outlier_log_df <- if (length(outlier_log) > 0) bind_rows(outlier_log) else tibble()

  if (nrow(outlier_log_df) > 0) {
    readr::write_csv(outlier_log_df, outlier_csv_path)
  }

  # Save PPT
  print(ppt, target = ppt_path)

  message("Done: ", branch_name, " / ", subbranch_name)
  message("  PPT: ", ppt_path)
  message(
    "  Outlier log: ",
    ifelse(file.exists(outlier_csv_path), outlier_csv_path, "No outliers removed")
  )

  invisible(list(
    ppt_path = ppt_path,
    outlier_csv_path = ifelse(file.exists(outlier_csv_path), outlier_csv_path, NA_character_),
    png_dir = png_dir
  ))
}

# -----------------------------------------------------------------------------
# PRE-LOOP CHECKS
# -----------------------------------------------------------------------------
activity_export_results <- list()

cat("\n===== PRE-LOOP CHECKS =====\n")
cat("nrow(fish_activity_summary): ", nrow(fish_activity_summary), "\n")
cat("Groups:\n")
print(table(as.character(fish_activity_summary$Group), useNA = "ifany"))
cat("Treatments:\n")
print(table(as.character(fish_activity_summary$Treatment), useNA = "ifany"))
cat("Timepoints:\n")
print(table(as.character(fish_activity_summary$timepoint), useNA = "ifany"))

cat("indicator_vars length: ", length(indicator_vars), "\n")
print(indicator_vars)

cat("branch_specs length: ", length(branch_specs), "\n")
print(lapply(branch_specs, `[[`, "branch_name"))

cat("subbranch_specs length: ", length(subbranch_specs), "\n")
print(lapply(subbranch_specs, `[[`, "subbranch_name"))

# OPT-4: Set up parallel workers for plot generation (moved to top level for proper scoping)
n_workers_step4 <- if (exists("STEP4_PARALLEL_WORKERS"))
  as.integer(STEP4_PARALLEL_WORKERS) else 0L
if (n_workers_step4 > 0L) {
  message("[OPT-4] Setting up parallel plan with ", n_workers_step4, " workers")
  future::plan(future::multisession(workers = n_workers_step4))
}

# -----------------------------------------------------------------------------
# RUN ALL EXPORTS
# -----------------------------------------------------------------------------
for (b in branch_specs) {
  for (s in subbranch_specs) {
    res_key <- paste(b$branch_name, s$subbranch_name, sep = "__")
    
    activity_export_results[[res_key]] <- run_activity_branch_subbranch(
      df_summary       = fish_activity_summary,
      branch_name      = b$branch_name,
      group_filter     = b$group_filter,
      subbranch_name   = s$subbranch_name,
      timepoint_filter = s$timepoint_filter,
      main_export_dir  = MAIN_EXPORT_DIR
    )
  }
}

assign("activity_export_results", activity_export_results, envir = .GlobalEnv)

# OPT-4: Reset parallel plan to sequential so STEP5 starts with a clean slate.
if (n_workers_step4 > 0L) {
  future::plan(future::sequential())
  message("[OPT-4] Parallel plan reset to sequential after STEP4 plot generation.")
}

.got_write_run_settings(
  MAIN_EXPORT_DIR,
  TRIAL_SUMMARY_PATH = TRIAL_SUMMARY_PATH
)
.got_write_session_info(MAIN_EXPORT_DIR)

message("Done.")
message("Main output folder: ", MAIN_EXPORT_DIR)
message("Summary CSV:        ", SUMMARY_CSV_PATH)

getwd()
MAIN_EXPORT_DIR
dir.exists(MAIN_EXPORT_DIR)
list.files(MAIN_EXPORT_DIR, recursive = TRUE, full.names = TRUE)
