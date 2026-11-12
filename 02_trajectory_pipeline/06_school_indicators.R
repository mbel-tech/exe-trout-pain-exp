# =============================================================================
# GROUP DYNAMICS PIPELINE — GOT ACTIVITY ANALYSIS
# =============================================================================
# This script computes group-level behavioural indicators from the tracking
# data produced by Steps 1 and 2.  It follows the same statistical and
# visualisation conventions as Steps 3–5.
#
# INDICATORS COMPUTED
# -------------------
# Fish-level (one value per fish per session timepoint):
#   1. Mean nearest-neighbour distance (NND, cm)
#      How close each fish was on average to its nearest companion.
#
# Session-level (one value per whole group per session timepoint):
#   2. Mean group polarization (0–1)
#      Alignment of swimming directions. 0 = random; 1 = all aligned.
#   3. Mean angular momentum (0–1)
#      Coordinated rotation of the group.
#   4. Mean inter-individual distance (IID, cm)
#      Average distance between all pairs of fish.
#   5. Zone synchrony (0–1)
#      Fraction of frames in which ALL fish were in the same zone.
#   6. Group centroid speed (cm/s)
#      Speed of the geometric centre of the group.
#   7. Number of collective motion events
#      Bouts where the whole group swam together (polarization > 0.7 for ≥ 1 s).
#   8. Mean collective motion event duration (s)
#   9. Leadership index
#      Fraction of frames led by the most consistent front fish.
#
# REQUIRED OBJECTS IN THE R ENVIRONMENT (supplied by pipeline Steps 1–2):
#   master_fish_by_frame
#   master_group_dynamics_by_frame  (or master_gd_grp)
#   master_group_dynamics_pairwise_by_frame  (or master_gd_pw)
#   master_group_dynamics_ind_by_frame  (or master_gd_ind)
#   master_col_mot_events_per_session  (or master_col_events)
#   master_session_event_summary  (or master_evt_summary)
#   master_h_stab  (optional — recomputed if absent)
#
# OUTPUT FOLDER:
#   GROUP_DYNAMICS_output_<timestamp>/
# =============================================================================


# =============================================================================
# 0) CONFIG — only section you need to edit
# =============================================================================
options(stringsAsFactors = FALSE)
options(na.action = "na.fail")
set.seed(123)

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
# Normally available in .GlobalEnv via 00_run_trajectory_pipeline.R.
# Fallback below only fires when this script is run standalone.
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
            " — using built-in labels only.")
  }
  rm(.this_script_dir, .pl_path)
}

# OPT-0: derive scope switches from master pipeline globals when available.
USE_TP1 <- if (exists("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE)) {
  isTRUE(get("KEEP_TP1", envir = .GlobalEnv, inherits = FALSE))
} else {
  FALSE   # standalone default
}

EXCLUDE_PHARMA_GROUPS <- if (exists("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE)) {
  isTRUE(get("DROP_PHARMA", envir = .GlobalEnv, inherits = FALSE))
} else {
  TRUE    # standalone default
}

# Include all experimental groups?  Derived from EXCLUDE_PHARMA_GROUPS.
USE_ALL_GROUPS <- !EXCLUDE_PHARMA_GROUPS

# Groups kept when USE_ALL_GROUPS = FALSE
GROUPS_RESTRICTED <- c("Exercise choice", "No choice")

# Statistical significance threshold
ALPHA <- 0.05

# ---------- output folders ----------
RUN_STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUT_DIR   <- .got_make_output_dir("STEP6_output")

DIR_QC      <- file.path(OUT_DIR, "00_qc")
DIR_DATA    <- file.path(OUT_DIR, "01_data")
DIR_AICC    <- file.path(OUT_DIR, "02_model_selection")
DIR_ANOVA   <- file.path(OUT_DIR, "03_anova")
DIR_POSTHOC <- file.path(OUT_DIR, "04_posthoc")
DIR_DIAG    <- file.path(OUT_DIR, "05_diagnostics")
DIR_GRAPHS  <- file.path(OUT_DIR, "06_graphs")
DIR_REPORTS <- file.path(OUT_DIR, "07_reports")

# OUT_DIR is already created by .got_make_output_dir(); create subfolders only.
for (d in c(DIR_QC, DIR_DATA, DIR_AICC, DIR_ANOVA,
            DIR_POSTHOC, DIR_DIAG, DIR_GRAPHS, DIR_REPORTS)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Derived constants
ANALYSIS_TIMEPOINTS <- if (isTRUE(USE_TP1)) c(1L, 2L, 3L) else c(2L, 3L)
TIMEPOINT_LABEL     <- if (isTRUE(USE_TP1)) "TP1_TP3" else "TP2_TP3"
GROUP_LABEL         <- if (isTRUE(USE_ALL_GROUPS)) "ALL_GROUPS" else "ExeChoice_NoChoice"


# =============================================================================
# 1) PACKAGES
# =============================================================================
required_pkgs <- c(
  "dplyr", "tidyr", "tibble", "stringr",
  "lme4", "lmerTest", "MuMIn", "car", "emmeans", "multcomp",
  "glmmTMB",
  "ggplot2", "officer", "rvg", "flextable",
  "data.table"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_pkgs) > 0L) {
  stop("Install required packages first:\n  ",
       paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr);   library(tidyr);   library(tibble);  library(stringr)
  library(lme4);    library(MuMIn);   library(car);
  library(emmeans); library(multcomp)
  library(ggplot2); library(officer); library(rvg);     library(flextable)
  library(data.table)
})

# Prevent name-masking surprises
select    <- dplyr::select
filter    <- dplyr::filter
mutate    <- dplyr::mutate
summarise <- dplyr::summarise

# Sum-to-zero contrasts — required for Type III ANOVA
options(contrasts = c("contr.sum", "contr.poly"))


# =============================================================================
# 2) LOAD DATA FROM GLOBAL ENVIRONMENT
# =============================================================================
# Each object below was created by Step 1 (idtracker parsing) or Step 2
# (activity annotation).  The master pipeline guarantees they exist before
# this step runs.

get_global <- function(primary, alternatives = character(0)) {
  # Try the primary name, then any alternative names provided.
  for (nm in c(primary, alternatives)) {
    if (exists(nm, envir = .GlobalEnv, inherits = FALSE)) {
      obj <- get(nm, envir = .GlobalEnv)
      if (is.data.frame(obj) && nrow(obj) > 0L) {
        message("Using GlobalEnv object: ", nm, " [", nrow(obj), " rows]")
        return(obj)
      }
    }
  }
  warning(
    "Object not found or empty in GlobalEnv: ",
    paste(c(primary, alternatives), collapse = " / "),
    "\nSome indicators that rely on this object will be skipped.",
    call. = FALSE
  )
  NULL
}

# Core tracking table — required for zone synchrony and centroid speed
mfbf <- get_global("master_fish_by_frame")
if (is.null(mfbf)) {
  # Disk fallback: when starting from STEP 5+, STEP 2 never ran so the object
  # is not in GlobalEnv. Search the most recent STEP2_output run first (using
  # the pipeline script directory, which is the canonical output location),
  # then fall back to getwd(), then to the master_outputs location.
  # Priority: .pipeline_dir > getwd() > PIPELINE_DATA_DIR/master_outputs
  .pick_latest_step2 <- function(root) {
    if (!dir.exists(root)) return(character(0L))
    dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
    if (length(dirs) == 0L) return(character(0L))
    # Sort by directory *name* (not full path) to get the newest timestamp
    dirs[order(basename(dirs), decreasing = TRUE)]
  }
  .step2_dirs <- if (exists(".pipeline_dir", inherits = TRUE)) {
    .pick_latest_step2(file.path(.pipeline_dir, "STEP2_output"))
  } else character(0L)
  # Fall back to getwd() if .pipeline_dir yielded nothing
  if (length(.step2_dirs) == 0L)
    .step2_dirs <- .pick_latest_step2(file.path(getwd(), "STEP2_output"))
  .mfbf_candidates <- c(
    if (length(.step2_dirs) > 0L)
      file.path(.step2_dirs[1L], "master_fish_by_frame.csv"),
    file.path(PIPELINE_DATA_DIR, "master_outputs", "by_frame", "master_fish_by_frame.csv")
  )
  .mfbf_candidates <- .mfbf_candidates[file.exists(.mfbf_candidates)]
  if (length(.mfbf_candidates) > 0L) {
    .csv_path <- .mfbf_candidates[1L]
    message("[STEP 6] master_fish_by_frame not in GlobalEnv — loading from disk:\n  ", .csv_path)
    mfbf <- data.table::fread(.csv_path, data.table = FALSE)
    assign("master_fish_by_frame", mfbf, envir = .GlobalEnv)
    message("[STEP 6] Loaded ", nrow(mfbf), " rows (", ncol(mfbf), " cols) from disk.")
  } else {
    stop(
      "master_fish_by_frame is required but not found in the R environment or on disk.\n",
      "Run STEP 2 before STEP 6, or ensure master_fish_by_frame.csv exists in STEP2_output/.",
      call. = FALSE
    )
  }
}

# Group dynamics tables produced by swaRmverse in Step 1
gd_grp  <- get_global("master_group_dynamics_by_frame",
                       "master_gd_grp")
gd_pw   <- get_global("master_group_dynamics_pairwise_by_frame",
                       "master_gd_pw")
gd_ind  <- get_global("master_group_dynamics_ind_by_frame",
                       "master_gd_ind")

# Disk fallbacks for group dynamics tables (needed when starting from STEP 5+)
.gd_by_frame_root <- file.path(PIPELINE_DATA_DIR, "master_outputs", "by_frame")
.load_gd_csv <- function(obj, subdir, filename, env_name) {
  if (!is.null(obj)) return(obj)
  path1 <- file.path(.gd_by_frame_root, subdir, filename)
  if (file.exists(path1)) {
    message("[STEP 6] Loading ", env_name, " from disk:\n  ", path1)
    tbl <- data.table::fread(path1, data.table = FALSE)
    assign(env_name, tbl, envir = .GlobalEnv)
    message("[STEP 6] Loaded ", nrow(tbl), " rows.")
    return(tbl)
  }
  obj
}
gd_grp <- .load_gd_csv(gd_grp, "group_group",    "master_group_dynamics_by_frame.csv",           "master_group_dynamics_by_frame")
gd_pw  <- .load_gd_csv(gd_pw,  "group_pairwise", "master_group_dynamics_pairwise_by_frame.csv",  "master_group_dynamics_pairwise_by_frame")
gd_ind <- .load_gd_csv(gd_ind, "group_ind",      "master_group_dynamics_ind_by_frame.csv",        "master_group_dynamics_ind_by_frame")
## DISABLED: collective motion and h-stability removed from pipeline outputs.
# evts    <- get_global("master_col_mot_events_per_session",
#                       "master_col_events")
# evtsumm <- get_global("master_session_event_summary",
#                       "master_evt_summary")
# h_stab  <- get_global("master_h_stab")
evts    <- NULL
evtsumm <- NULL
h_stab  <- NULL


# =============================================================================
# 3) METADATA LOOKUP
# =============================================================================
# The group dynamics tables are keyed by trial_id + segment_id (timepoint).
# We extract Group/Treatment/trial_date from master_fish_by_frame, which has
# those columns after Step 2.

required_meta <- c("trial_id", "timepoint", "Group", "Treatment", "trial_date")
missing_meta  <- setdiff(required_meta, names(mfbf))
if (length(missing_meta) > 0L) {
  stop("master_fish_by_frame is missing metadata columns: ",
       paste(missing_meta, collapse = ", "), call. = FALSE)
}

meta_lookup <- mfbf %>%
  mutate(timepoint = as.integer(as.character(timepoint))) %>%
  distinct(trial_id, timepoint, Group, Treatment, trial_date) %>%
  filter(!is.na(trial_id), !is.na(timepoint))


# =============================================================================
# 4) SMALL UTILITY FUNCTIONS
# =============================================================================

safe_mean <- function(x) { x <- x[is.finite(x)]; if (!length(x)) NA_real_ else mean(x) }
safe_sd   <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2L) NA_real_ else sd(x) }
safe_sem  <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2L) NA_real_ else sd(x)/sqrt(length(x)) }

# Find the first matching column name from a list of candidates
pick_col <- function(df, candidates) {
  intersect(candidates, names(df))[1L]
}


# =============================================================================
# 5: PRE-PROCESS PAIRWISE TABLE
# =============================================================================
# swaRmverse encodes session as "1970-01-01_session_trialN_SegmentM" in the
# 'set' column.  When trial_id / segment_id are absent, parse them out so that
# all downstream indicator blocks can use standard column names.

if (!is.null(gd_pw) && nrow(gd_pw) > 0L &&
    "set" %in% names(gd_pw) && !"trial_id" %in% names(gd_pw)) {
  gd_pw$trial_id   <- suppressWarnings(as.integer(
    sub(".*_trial(\\d+)_Segment\\d+$", "\\1", gd_pw$set)))
  gd_pw$segment_id <- suppressWarnings(as.integer(
    sub(".*_Segment(\\d+)$", "\\1", gd_pw$set)))
  n_ok <- sum(!is.na(gd_pw$trial_id))
  message("[STEP 6] Parsed 'set' -> trial_id/segment_id for pairwise table: ",
          n_ok, " / ", nrow(gd_pw), " rows resolved.")
}


# =============================================================================
# 5A) INDICATOR: NEAREST-NEIGHBOUR DISTANCE (NND)
# =============================================================================
# NND for each fish in each frame = shortest distance to any other fish.
# We average that across all frames within a fish × session × timepoint.
# Source: master_group_dynamics_pairwise_by_frame

nnd_fish_tp <- NULL

if (!is.null(gd_pw) && nrow(gd_pw) > 0L) {

  # Identify column names — swaRmverse versions use different names
  # 'nnd' is the per-fish nearest-neighbour distance already computed by swaRmverse
  id1_col  <- pick_col(gd_pw, c("id", "id1", "id_i", "ind"))
  dist_col <- pick_col(gd_pw, c("nnd", "dist", "distance", "dist_cm", "IID"))
  t_col    <- pick_col(gd_pw, c("frame", "t"))
  seg_col  <- pick_col(gd_pw, c("segment_id", "timepoint"))

  if (!is.na(id1_col) && !is.na(dist_col) &&
      !is.na(t_col)   && !is.na(seg_col)  &&
      "trial_id" %in% names(gd_pw)) {

    # For every fish × frame × session: minimum distance = NND
    nnd_by_frame <- gd_pw %>%
      rename(fish_id = !!id1_col,
             dist_pw = !!dist_col,
             tp_key  = !!seg_col) %>%
      filter(is.finite(dist_pw), dist_pw >= 0) %>%
      group_by(trial_id, tp_key, fish_id, .data[[t_col]]) %>%
      summarise(nnd_cm = min(dist_pw, na.rm = TRUE), .groups = "drop")

    # Average NND per fish × session × timepoint
    nnd_fish_tp <- nnd_by_frame %>%
      group_by(trial_id, tp_key, fish_id) %>%
      summarise(
        mean_nnd_cm  = safe_mean(nnd_cm),
        sd_nnd_cm    = safe_sd(nnd_cm),
        n_frames_nnd = sum(is.finite(nnd_cm)),
        .groups = "drop"
      ) %>%
      rename(timepoint = tp_key) %>%
      mutate(timepoint = as.integer(timepoint)) %>%
      left_join(meta_lookup, by = c("trial_id", "timepoint"))

    # swaRmverse reports distances in the same units as input positions (pixels).
    # Convert to cm using the trial-level px_per_cm from master_fish_by_frame.
    if ("length_unit_px_per_cm" %in% names(mfbf)) {
      .px_lut <- mfbf %>%
        mutate(timepoint = as.integer(as.character(timepoint))) %>%
        group_by(trial_id, timepoint) %>%
        summarise(px_per_cm = mean(as.numeric(length_unit_px_per_cm), na.rm = TRUE),
                  .groups = "drop")
      nnd_fish_tp <- nnd_fish_tp %>%
        left_join(.px_lut, by = c("trial_id", "timepoint")) %>%
        mutate(
          mean_nnd_cm = mean_nnd_cm / px_per_cm,
          sd_nnd_cm   = sd_nnd_cm   / px_per_cm
        ) %>%
        select(-px_per_cm)
    }

    message("NND computed: ", nrow(nnd_fish_tp), " fish x timepoint rows")
    write.csv(nnd_fish_tp,
              file.path(DIR_DATA, "NND_fish_by_timepoint.csv"),
              row.names = FALSE)

  } else {
    warning("Pairwise table missing required columns. NND skipped.", call. = FALSE)
  }

} else {
  message("Pairwise group-dynamics table not available. NND skipped.")
}


# =============================================================================
# 5B) INDICATOR: INTER-INDIVIDUAL DISTANCE (IID)
# =============================================================================
# IID computed from fish positions in master_fish_by_frame (x_interp, y_interp).
# Session-level: mean of all pairwise Euclidean distances per session x timepoint
#   (averaged over frames).
# Fish-level: per-fish mean distance to all other fish, averaged over frames
#   (used in the individual-level analysis in section 8B).
# data.table is used for efficient per-frame grouped distance computation.

iid_session_tp <- NULL
iid_fish_tp    <- NULL

iid_required <- c("trial_id","timepoint","frame","fish_id",
                   "x_interp","y_interp","length_unit_px_per_cm")

if (all(iid_required %in% names(mfbf))) {

  dt_iid <- data.table::as.data.table(mfbf)[
    is.finite(suppressWarnings(as.numeric(x_interp))) &
    is.finite(suppressWarnings(as.numeric(y_interp))) &
    suppressWarnings(as.numeric(length_unit_px_per_cm)) > 0,
    .(trial_id,
      timepoint = as.integer(as.character(timepoint)),
      frame, fish_id,
      x     = as.numeric(x_interp),
      y     = as.numeric(y_interp),
      px_cm = as.numeric(length_unit_px_per_cm))
  ]

  # Per frame: compute distance matrix; extract fish-level row means and
  # session-level mean of lower triangle.  A single dist() call per frame group.
  iid_by_frame <- dt_iid[, {
    n <- .N
    if (n < 2L) {
      .(fish_id          = fish_id,
        fish_iid_cm      = rep(NA_real_, n),
        session_iid_cm   = rep(NA_real_, n),
        n_fish           = rep(n, n))
    } else {
      dm       <- as.matrix(stats::dist(cbind(x, y)))
      diag(dm) <- NA_real_
      px_avg   <- mean(px_cm)
      .(fish_id        = fish_id,
        fish_iid_cm    = rowMeans(dm, na.rm = TRUE) / px_avg,
        session_iid_cm = rep(mean(dm[lower.tri(dm)], na.rm = TRUE) / px_avg, n),
        n_fish         = rep(n, n))
    }
  }, by = .(trial_id, timepoint, frame)]

  # Session-level IID: first value per frame (session_iid_cm is identical for
  # all fish in the same frame), then average across frames.
  iid_sess_raw <- unique(iid_by_frame[
    !is.na(session_iid_cm),
    .(session_iid_cm = session_iid_cm[1L], n_fish = n_fish[1L]),
    by = .(trial_id, timepoint, frame)])

  iid_session_tp <- as.data.frame(
    merge(
      iid_sess_raw[n_fish >= 2L, .(
        mean_iid_cm  = mean(session_iid_cm, na.rm = TRUE),
        sd_iid_cm    = sd(session_iid_cm,   na.rm = TRUE),
        n_frames_iid = .N
      ), by = .(trial_id, timepoint)],
      data.table::as.data.table(meta_lookup),
      by = c("trial_id", "timepoint"), all.x = TRUE)
  )

  # Fish-level IID: average fish_iid_cm across frames per fish x session x timepoint.
  iid_fish_tp <- as.data.frame(
    merge(
      iid_by_frame[!is.na(fish_iid_cm), .(
        mean_iid_cm  = mean(fish_iid_cm, na.rm = TRUE),
        sd_iid_cm    = sd(fish_iid_cm,   na.rm = TRUE),
        n_frames_iid = .N
      ), by = .(trial_id, timepoint, fish_id)],
      data.table::as.data.table(meta_lookup),
      by = c("trial_id", "timepoint"), all.x = TRUE)
  )

  message("IID computed: ", nrow(iid_session_tp), " session x timepoint rows; ",
          nrow(iid_fish_tp), " fish x timepoint rows")
  write.csv(iid_session_tp,
            file.path(DIR_DATA, "IID_session_by_timepoint.csv"), row.names = FALSE)
  write.csv(iid_fish_tp,
            file.path(DIR_DATA, "IID_fish_by_timepoint.csv"),    row.names = FALSE)

} else {
  warning("IID skipped — master_fish_by_frame missing columns: ",
          paste(setdiff(iid_required, names(mfbf)), collapse = ", "), call. = FALSE)
}


# =============================================================================
# 5C) INDICATORS: GROUP POLARIZATION AND ANGULAR MOMENTUM
# =============================================================================
# Polarization (0–1): how well fish align their swimming directions.
# Angular momentum (0–1): how much the group rotates together.
# Both are computed per frame in Step 1; here we average over each session.
# Source: master_group_dynamics_by_frame

pol_ang_session_tp <- NULL

if (!is.null(gd_grp) && nrow(gd_grp) > 0L) {

  seg_grp <- pick_col(gd_grp, c("segment_id","timepoint"))

  if (!is.na(seg_grp) && "trial_id" %in% names(gd_grp)) {

    # Resolve column names defensively
    pol_col <- if ("polarization" %in% names(gd_grp)) "polarization" else
      grep("pol", names(gd_grp), value = TRUE, ignore.case = TRUE)[1L]
    ang_col <- if ("ang_momentum" %in% names(gd_grp)) "ang_momentum" else
      grep("rotat|ang.*mom|momentum", names(gd_grp), value = TRUE, ignore.case = TRUE)[1L]

    df_grp <- gd_grp %>%
      rename(tp_key = !!seg_grp) %>%
      mutate(tp_key = as.integer(tp_key))

    # Add polarization column, or NA if not found
    if (!is.na(pol_col) && pol_col %in% names(df_grp)) {
      df_grp <- df_grp %>% mutate(polarization = as.numeric(.data[[pol_col]]))
    } else {
      df_grp <- df_grp %>% mutate(polarization = NA_real_)
    }

    # Add angular momentum column, or NA if not found
    if (!is.na(ang_col) && ang_col %in% names(df_grp)) {
      df_grp <- df_grp %>% mutate(ang_momentum = as.numeric(.data[[ang_col]]))
    } else {
      df_grp <- df_grp %>% mutate(ang_momentum = NA_real_)
    }

    pol_ang_session_tp <- df_grp %>%
      group_by(trial_id, tp_key) %>%
      summarise(
        mean_polarization = safe_mean(polarization),
        sd_polarization   = safe_sd(polarization),
        mean_ang_momentum = safe_mean(ang_momentum),
        sd_ang_momentum   = safe_sd(ang_momentum),
        n_frames_pol      = sum(is.finite(polarization)),
        .groups = "drop"
      ) %>%
      rename(timepoint = tp_key) %>%
      left_join(meta_lookup, by = c("trial_id", "timepoint"))

    message("Polarization / angular momentum: ", nrow(pol_ang_session_tp), " rows")
    write.csv(pol_ang_session_tp,
              file.path(DIR_DATA, "polarization_angmom_session_by_timepoint.csv"),
              row.names = FALSE)
  }

} else {
  message("Group-level dynamics table not available. Polarization skipped.")
}


# =============================================================================
# 5D) INDICATORS: COLLECTIVE MOTION EVENTS — DISABLED
# =============================================================================
## DISABLED: collective-motion event indicators removed from pipeline outputs.
## Re-enable by un-commenting the block below and restoring evts/evtsumm loading.
events_session_tp <- NULL
# build_events <- function(df, seg_candidates) {
#   seg <- pick_col(df, seg_candidates)
#   if (is.na(seg) || !("trial_id" %in% names(df))) return(NULL)
#   df %>%
#     rename(tp_key = !!seg) %>%
#     mutate(timepoint = as.integer(tp_key)) %>%
#     left_join(meta_lookup, by = c("trial_id", "timepoint"))
# }
# if (!is.null(evtsumm) && nrow(evtsumm) > 0L) {
#   seg_ev <- pick_col(evtsumm, c("segment_id","timepoint"))
#   if (!is.na(seg_ev) && "trial_id" %in% names(evtsumm)) {
#     events_raw <- evtsumm %>%
#       rename(tp_key = !!seg_ev) %>%
#       mutate(
#         timepoint        = as.integer(tp_key),
#         n_col_mot_events = as.integer(n_events),
#         mean_event_dur_s = as.numeric(avg_duration)
#       ) %>%
#       select(trial_id, timepoint, n_col_mot_events, mean_event_dur_s)
#     events_session_tp <- meta_lookup %>%
#       left_join(events_raw, by = c("trial_id", "timepoint")) %>%
#       mutate(n_col_mot_events = if_else(is.na(n_col_mot_events), 0L, n_col_mot_events))
#     message("Collective motion events: ", nrow(events_session_tp), " rows")
#     write.csv(events_session_tp,
#               file.path(DIR_DATA, "col_motion_events_session_by_timepoint.csv"),
#               row.names = FALSE)
#   }
# } else if (!is.null(evts) && nrow(evts) > 0L) {
#   seg_ev2 <- pick_col(evts, c("segment_id","timepoint"))
#   if (!is.na(seg_ev2) && "trial_id" %in% names(evts) && "duration_s" %in% names(evts)) {
#     events_session_tp <- evts %>%
#       rename(tp_key = !!seg_ev2) %>%
#       mutate(timepoint = as.integer(tp_key)) %>%
#       group_by(trial_id, timepoint) %>%
#       summarise(
#         n_col_mot_events = n(),
#         mean_event_dur_s = safe_mean(duration_s),
#         .groups = "drop"
#       ) %>%
#       right_join(meta_lookup, by = c("trial_id", "timepoint")) %>%
#       mutate(n_col_mot_events = if_else(is.na(n_col_mot_events), 0L, n_col_mot_events))
#     message("Collective motion events rebuilt from raw table: ",
#             nrow(events_session_tp), " rows")
#     write.csv(events_session_tp,
#               file.path(DIR_DATA, "col_motion_events_session_by_timepoint.csv"),
#               row.names = FALSE)
#   }
# } else {
#   message("Collective motion events table not available. Skipped.")
# }


# =============================================================================
# 5E) INDICATOR: ZONE SYNCHRONY
# =============================================================================
# Fraction of frames in which ALL fish were simultaneously in the same zone.
# A value of 1 means the group always moved between zones together.
# Source: master_fish_by_frame

synchrony_session_tp <- NULL

sync_required <- c("trial_id","timepoint","frame","fish_id","zone")

if (all(sync_required %in% names(mfbf))) {

  sync_by_frame <- mfbf %>%
    filter(!is.na(zone), nzchar(trimws(zone))) %>%
    mutate(
      timepoint   = as.integer(as.character(timepoint)),
      zone_lower  = tolower(trimws(zone))
    ) %>%
    group_by(trial_id, timepoint, frame) %>%
    summarise(
      n_distinct_zones = n_distinct(zone_lower),
      .groups = "drop"
    ) %>%
    mutate(synchronous = n_distinct_zones == 1L)

  synchrony_session_tp <- sync_by_frame %>%
    group_by(trial_id, timepoint) %>%
    summarise(
      zone_synchrony     = mean(synchronous, na.rm = TRUE),
      n_frames_sync      = n(),
      .groups = "drop"
    ) %>%
    left_join(meta_lookup, by = c("trial_id", "timepoint"))

  message("Zone synchrony: ", nrow(synchrony_session_tp), " rows")
  write.csv(synchrony_session_tp,
            file.path(DIR_DATA, "zone_synchrony_session_by_timepoint.csv"),
            row.names = FALSE)

} else {
  warning("Zone synchrony skipped — missing columns: ",
          paste(setdiff(sync_required, names(mfbf)), collapse = ", "), call. = FALSE)
}


# =============================================================================
# 5F) INDICATOR: GROUP CENTROID SPEED
# =============================================================================
# Mean (x, y) position of all fish per frame = group centroid.
# Centroid speed = frame-to-frame displacement converted to cm/s.
# Captures how fast the group as a whole moved through the arena.
# Source: master_fish_by_frame

centroid_session_tp <- NULL

cent_required <- c("trial_id","timepoint","frame","fish_id",
                   "x_interp","y_interp","length_unit_px_per_cm","time")

if (all(cent_required %in% names(mfbf))) {

  # Step 1 — compute mean centroid position per frame
  centroid_by_frame <- mfbf %>%
    mutate(
      timepoint             = as.integer(as.character(timepoint)),
      x_interp              = as.numeric(x_interp),
      y_interp              = as.numeric(y_interp),
      length_unit_px_per_cm = as.numeric(length_unit_px_per_cm),
      time                  = as.numeric(time)
    ) %>%
    filter(is.finite(x_interp), is.finite(y_interp), length_unit_px_per_cm > 0) %>%
    group_by(trial_id, timepoint, frame) %>%
    summarise(
      cx_px     = mean(x_interp, na.rm = TRUE),
      cy_px     = mean(y_interp, na.rm = TRUE),
      px_per_cm = mean(length_unit_px_per_cm, na.rm = TRUE),
      time_s    = first(time),
      .groups = "drop"
    )

  # Step 2 — frame-to-frame centroid displacement → cm/s
  centroid_speed_by_frame <- centroid_by_frame %>%
    arrange(trial_id, timepoint, frame) %>%
    group_by(trial_id, timepoint) %>%
    mutate(
      dcx   = cx_px  - lag(cx_px),
      dcy   = cy_px  - lag(cy_px),
      dt    = time_s - lag(time_s),
      # Convert pixel displacement to cm then divide by dt to get cm/s
      centroid_speed_cm_s = if_else(
        is.finite(dcx) & is.finite(dcy) & is.finite(dt) & dt > 0 & px_per_cm > 0,
        sqrt(dcx^2 + dcy^2) / px_per_cm / dt,
        NA_real_
      )
    ) %>%
    ungroup()

  # Step 3 — average over all frames per session × timepoint
  centroid_session_tp <- centroid_speed_by_frame %>%
    group_by(trial_id, timepoint) %>%
    summarise(
      mean_centroid_speed_cm_s = safe_mean(centroid_speed_cm_s),
      sd_centroid_speed_cm_s   = safe_sd(centroid_speed_cm_s),
      n_frames_centroid        = sum(is.finite(centroid_speed_cm_s)),
      .groups = "drop"
    ) %>%
    left_join(meta_lookup, by = c("trial_id", "timepoint"))

  message("Group centroid speed: ", nrow(centroid_session_tp), " rows")
  write.csv(centroid_session_tp,
            file.path(DIR_DATA, "centroid_speed_session_by_timepoint.csv"),
            row.names = FALSE)

} else {
  warning("Centroid speed skipped — missing columns: ",
          paste(setdiff(cent_required, names(mfbf)), collapse = ", "), call. = FALSE)
}


# =============================================================================
# 5G) INDICATOR: CONVEX HULL AREA
# =============================================================================
# Per frame: compute the convex hull of all fish positions and its area
# (Shoelace formula).  Average over frames per session x timepoint.
# Units: cm² (converted from pixels using length_unit_px_per_cm).
# Requires >= 3 fish per frame to form a non-degenerate hull.
# Source: master_fish_by_frame

hull_session_tp <- NULL

hull_required <- c("trial_id","timepoint","frame","fish_id",
                   "x_interp","y_interp","length_unit_px_per_cm")

if (all(hull_required %in% names(mfbf))) {

  dt_hull <- data.table::as.data.table(mfbf)[
    is.finite(suppressWarnings(as.numeric(x_interp))) &
    is.finite(suppressWarnings(as.numeric(y_interp))) &
    suppressWarnings(as.numeric(length_unit_px_per_cm)) > 0,
    .(trial_id,
      timepoint = as.integer(as.character(timepoint)),
      frame,
      x     = as.numeric(x_interp),
      y     = as.numeric(y_interp),
      px_cm = as.numeric(length_unit_px_per_cm))
  ]

  hull_by_frame <- dt_hull[, {
    n <- .N
    if (n < 3L) {
      .(hull_area_cm2 = NA_real_, n_fish = n)
    } else {
      hull_idx <- tryCatch(grDevices::chull(x, y), error = function(e) integer(0))
      if (length(hull_idx) < 3L) {
        .(hull_area_cm2 = 0, n_fish = n)
      } else {
        hx <- x[hull_idx]; hy <- y[hull_idx]; nh <- length(hx)
        area_px2 <- 0.5 * abs(
          sum(hx * c(hy[-1L], hy[1L]) - c(hx[-1L], hx[1L]) * hy))
        .(hull_area_cm2 = area_px2 / (mean(px_cm)^2), n_fish = n)
      }
    }
  }, by = .(trial_id, timepoint, frame)]

  hull_session_tp <- as.data.frame(
    merge(
      hull_by_frame[n_fish >= 3L & is.finite(hull_area_cm2), .(
        mean_hull_area_cm2 = mean(hull_area_cm2, na.rm = TRUE),
        sd_hull_area_cm2   = sd(hull_area_cm2,   na.rm = TRUE),
        n_frames_hull      = .N
      ), by = .(trial_id, timepoint)],
      data.table::as.data.table(meta_lookup),
      by = c("trial_id", "timepoint"), all.x = TRUE)
  )

  message("Convex hull area: ", nrow(hull_session_tp), " rows")
  write.csv(hull_session_tp,
            file.path(DIR_DATA, "hull_area_session_by_timepoint.csv"),
            row.names = FALSE)

} else {
  warning("Convex hull area skipped — master_fish_by_frame missing columns: ",
          paste(setdiff(hull_required, names(mfbf)), collapse = ", "), call. = FALSE)
}


# =============================================================================
# 5H) INDICATOR: LEADERSHIP INDEX — DISABLED
# =============================================================================
## DISABLED: leadership index (h-stability) removed from pipeline outputs.
## Re-enable by un-commenting the block below and restoring h_stab loading.
leadership_session_tp <- NULL
# if (!is.null(h_stab) && nrow(h_stab) > 0L && "h_stab" %in% names(h_stab)) {
#   seg_h <- pick_col(h_stab, c("segment_id","timepoint"))
#   if (!is.na(seg_h) && "trial_id" %in% names(h_stab)) {
#     leadership_session_tp <- h_stab %>%
#       rename(tp_key = !!seg_h) %>%
#       mutate(timepoint = as.integer(tp_key)) %>%
#       select(trial_id, timepoint, leadership_index = h_stab) %>%
#       left_join(meta_lookup, by = c("trial_id", "timepoint"))
#     message("Leadership index loaded from master_h_stab: ",
#             nrow(leadership_session_tp), " rows")
#     write.csv(leadership_session_tp,
#               file.path(DIR_DATA, "leadership_session_by_timepoint.csv"),
#               row.names = FALSE)
#   }
# } else if ("heading_rad" %in% names(mfbf)) {
#   message("master_h_stab not available — recomputing leadership from positions.")
#   tryCatch({
#     grp_heading <- mfbf %>%
#       mutate(timepoint = as.integer(as.character(timepoint)),
#              heading_rad = as.numeric(heading_rad)) %>%
#       filter(is.finite(heading_rad)) %>%
#       group_by(trial_id, timepoint, frame) %>%
#       summarise(
#         group_heading = atan2(mean(sin(heading_rad), na.rm = TRUE),
#                               mean(cos(heading_rad), na.rm = TRUE)),
#         .groups = "drop"
#       )
#     leaders <- mfbf %>%
#       mutate(timepoint = as.integer(as.character(timepoint)),
#              x_interp  = as.numeric(x_interp),
#              y_interp  = as.numeric(y_interp)) %>%
#       filter(is.finite(x_interp), is.finite(y_interp)) %>%
#       left_join(grp_heading, by = c("trial_id", "timepoint", "frame")) %>%
#       filter(is.finite(group_heading)) %>%
#       mutate(proj = x_interp * cos(group_heading) + y_interp * sin(group_heading)) %>%
#       group_by(trial_id, timepoint, frame) %>%
#       slice_max(proj, n = 1L, with_ties = FALSE) %>%
#       ungroup()
#     leadership_session_tp <- leaders %>%
#       group_by(trial_id, timepoint) %>%
#       summarise(leadership_index = max(table(fish_id)) / n(), .groups = "drop") %>%
#       left_join(meta_lookup, by = c("trial_id", "timepoint"))
#     message("Leadership index recomputed: ", nrow(leadership_session_tp), " rows")
#     write.csv(leadership_session_tp,
#               file.path(DIR_DATA, "leadership_session_by_timepoint.csv"),
#               row.names = FALSE)
#   }, error = function(e) {
#     warning("Leadership recomputation failed: ", e$message, call. = FALSE)
#     leadership_session_tp <<- NULL
#   })
# } else {
#   message("Leadership index skipped (heading_rad absent and master_h_stab unavailable).")
# }


# =============================================================================
# 6) MERGE SESSION-LEVEL SUMMARY TABLE
# =============================================================================
# Combine all session-level indicators into one table, then restrict to the
# chosen timepoints and groups.

# Lightweight join helper: adds only the new indicator columns from 'new_df'
add_indicator <- function(base, new_df, new_cols, by = c("trial_id","timepoint")) {
  if (is.null(new_df) || nrow(new_df) == 0L) return(base)
  avail <- intersect(new_cols, names(new_df))
  if (length(avail) == 0L) return(base)
  left_join(base, new_df %>% select(all_of(c(by, avail))), by = by)
}

session_gd_summary <- meta_lookup %>%
  filter(timepoint %in% ANALYSIS_TIMEPOINTS) %>%
  { if (!isTRUE(USE_ALL_GROUPS)) filter(., Group %in% GROUPS_RESTRICTED) else . } %>%
  add_indicator(pol_ang_session_tp,   c("mean_polarization","mean_ang_momentum")) %>%
  add_indicator(iid_session_tp,       c("mean_iid_cm")) %>%
  add_indicator(synchrony_session_tp, c("zone_synchrony")) %>%
  add_indicator(centroid_session_tp,  c("mean_centroid_speed_cm_s")) %>%
  add_indicator(hull_session_tp,      c("mean_hull_area_cm2")) %>%
  ## DISABLED: add_indicator(events_session_tp,    c("n_col_mot_events","mean_event_dur_s")) %>%
  ## DISABLED: add_indicator(leadership_session_tp,c("leadership_index")) %>%
  mutate(
    timepoint = factor(timepoint, levels = ANALYSIS_TIMEPOINTS),
    Group     = factor(Group,
                       levels = c("No choice","Exercise choice","Naloxone","Morphine")),
    Treatment = factor(Treatment, levels = c("Control","Exercise")),
    time_num  = as.numeric(as.character(timepoint))
  ) %>%
  droplevels()

write.csv(
  session_gd_summary,
  file.path(DIR_DATA,
            paste0("session_gd_summary_", TIMEPOINT_LABEL, "_", GROUP_LABEL, ".csv")),
  row.names = FALSE
)

message("Session summary: ", nrow(session_gd_summary), " rows × ",
        ncol(session_gd_summary), " columns")


# =============================================================================
# 7) NND FISH-LEVEL TABLE (filtered and factored)
# =============================================================================
nnd_model_data <- NULL

if (!is.null(nnd_fish_tp) && nrow(nnd_fish_tp) > 0L) {
  nnd_model_data <- nnd_fish_tp %>%
    filter(timepoint %in% ANALYSIS_TIMEPOINTS,
           is.finite(mean_nnd_cm)) %>%
    { if (!isTRUE(USE_ALL_GROUPS)) filter(., Group %in% GROUPS_RESTRICTED) else . } %>%
    mutate(
      timepoint = factor(timepoint, levels = ANALYSIS_TIMEPOINTS),
      Group     = factor(Group,
                         levels = c("No choice","Exercise choice","Naloxone","Morphine")),
      Treatment = factor(Treatment, levels = c("Control","Exercise")),
      time_num  = as.numeric(as.character(timepoint))
    ) %>%
    droplevels()

  write.csv(nnd_model_data,
            file.path(DIR_DATA,
                      paste0("NND_model_data_", TIMEPOINT_LABEL, "_", GROUP_LABEL, ".csv")),
            row.names = FALSE)
}


# =============================================================================
# 8) STATISTICAL MODELLING
# =============================================================================
# Same procedure as Steps 3 and 5:
#   1. Build a candidate set of LMMs combining different fixed-effect
#      structures (Group, Treatment, timepoint, and their interactions) with
#      two random-intercept options: (1|trial_id) and (1|trial_date).
#   2. Select the best model by AICc; prefer timepoint-as-factor over
#      time_num when both are within 2 AICc units.
#   3. Type III Wald chi-square ANOVA on the selected model.
#   4. Tukey-adjusted pairwise post-hoc tests where ANOVA is significant.

# RE forms for NND (per-fish outcome): use fish_id_matched_trial when idmatcher qualifies.
# Session-level outcomes (polarisation, zone-sync, centroid speed) must use RE_FORMS_NO_FISH.
RE_FORMS_FULL <- c(
  trial_id           = "(1|trial_id)",
  trial_date         = "(1|trial_date)",
  fish_matched       = "(1|trial_id) + (1|fish_id_matched_trial)",
  fish_matched_nested = "(1|trial_id/fish_id_matched_trial)"
)

RE_FORMS_NO_FISH <- c(
  trial_id   = "(1|trial_id)",
  trial_date = "(1|trial_date)"
)

RE_FORMS <- RE_FORMS_FULL   # overwritten by gate_re_forms_gd() below

# Gate fish RE candidates for NND based on idmatcher reliability tier (primary)
# or Hungarian BL-match tier (fallback). Session-level outcomes must NOT use
# RE_FORMS_FULL — call gate_re_forms_gd() only for NND.
gate_re_forms_gd <- function(bl_match_df = NULL, min_tier = NULL) {
  if (is.null(min_tier)) {
    min_tier <- if (exists("MIN_BL_RELIABILITY_TIER", envir = .GlobalEnv, inherits = FALSE))
      get("MIN_BL_RELIABILITY_TIER", envir = .GlobalEnv)
    else "reliable"
  }

  # ── Idmatcher path (primary) ─────────────────────────────────────────────
  idmatcher_min_tier <- if (exists("IDMATCHER_MIN_TIER", envir = .GlobalEnv, inherits = FALSE))
    get("IDMATCHER_MIN_TIER", envir = .GlobalEnv) else min_tier

  if (exists("idmatcher_tier_df", envir = .GlobalEnv, inherits = FALSE) &&
      is.data.frame(idmatcher_tier_df) &&
      "idmatcher_tier" %in% names(idmatcher_tier_df)) {
    qualifying_ids <- qualifying_trial_ids(idmatcher_tier_df, idmatcher_min_tier,
                                           tier_col = "idmatcher_tier")
    n_ok <- length(qualifying_ids)
    if (n_ok == 0) {
      message("STEP6: No trials meet idmatcher tier '", idmatcher_min_tier,
              "' — dropping fish RE candidates (NND will use trial_id RE only).")
      return(RE_FORMS_NO_FISH)
    }
    message("STEP6: ", n_ok, " trial(s) meet idmatcher tier '", idmatcher_min_tier,
            "' — using (1|fish_id_matched_trial) for NND.")
    return(RE_FORMS_FULL)
  }

  # ── Hungarian fallback ───────────────────────────────────────────────────
  if (is.null(bl_match_df) || !is.data.frame(bl_match_df)) {
    message("STEP6: No idmatcher or BL match table available — using full RE_FORMS.")
    return(RE_FORMS_FULL)
  }
  tier_col <- if ("hungarian_tier" %in% names(bl_match_df)) "hungarian_tier"
              else if ("reliability_tier" %in% names(bl_match_df)) "reliability_tier"
              else NULL
  if (is.null(tier_col)) {
    message("STEP6: BL match table has no tier column — using full RE_FORMS.")
    return(RE_FORMS_FULL)
  }
  tier_order <- c("unreliable", "borderline", "reliable", "conservative")
  min_idx    <- match(tolower(trimws(min_tier)), tier_order)
  if (is.na(min_idx)) min_idx <- 3L
  tier_vals  <- match(tolower(trimws(bl_match_df[[tier_col]])), tier_order)
  n_ok <- sum(!is.na(tier_vals) & tier_vals >= min_idx)
  if (n_ok == 0) {
    message("STEP6: No trials meet Hungarian tier '", min_tier,
            "' — dropping fish RE candidates from STEP6.")
    return(RE_FORMS_NO_FISH)
  }
  message("STEP6 (Hungarian fallback): ", n_ok, " trial(s) meet tier '", min_tier,
          "' — including fish_id_matched_trial RE in STEP6 NND.")
  RE_FORMS_FULL
}

# Pre-specified full 3-way fixed structure (AICc selects RE only).
GD_FIXED_FORMULA_FACTOR  <- "Group * Treatment * timepoint"
GD_FIXED_FORMULA_NUMERIC <- "Group * Treatment * time_num"

# Fit the pre-specified fixed structure with each RE candidate; return best by AICc.
# Differs from fit_candidate_set_gd() in that the fixed structure is NOT searched.
# use_gamma: when TRUE, fits Gamma LMM (glmmTMB) on the raw response instead of
#            Gaussian lmer on a (possibly transformed) column.
fit_candidate_set_gd_prespecified <- function(data, response, outcome_prefix,
                                               re_forms = RE_FORMS,
                                               use_gamma = FALSE) {
  candidates <- list()
  meta_rows  <- list()

  if (use_gamma && requireNamespace("glmmTMB", quietly = TRUE)) {
    # ---- Gamma LMM candidates -----------------------------------------------
    for (time_var in c("timepoint", "time_num")) {
      fixed_formula <- if (time_var == "timepoint") GD_FIXED_FORMULA_FACTOR else GD_FIXED_FORMULA_NUMERIC
      for (re_name in names(re_forms)) {
        re_form <- re_forms[[re_name]]
        if (re_name == "trial_date" && !"trial_date" %in% names(data)) next

        fml <- tryCatch(
          stats::as.formula(sprintf("%s ~ %s + %s", response, fixed_formula, re_form)),
          error = function(e) NULL
        )
        if (is.null(fml)) next

        mname <- paste(outcome_prefix, paste0("time=", time_var),
                       paste0("RE=", re_name), "Gamma", sep = "__")

        fit <- tryCatch(
          glmmTMB::glmmTMB(fml, data = data, family = Gamma(link = "log"),
                            control = glmmTMB::glmmTMBControl(
                              optCtrl = list(iter.max = 1000L))),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          aicc_g <- tryCatch(MuMIn::AICc(fit), error = function(e) NA_real_)
          if (is.finite(aicc_g)) {
            candidates[[mname]] <- fit
            meta_rows[[mname]] <- data.frame(
              model         = mname, time_var = time_var, re_name = re_name,
              fixed_formula = fixed_formula,
              singular      = FALSE,   # isSingular not applicable to glmmTMB here
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    if (length(candidates) == 0L) {
      warning("All Gamma LMM candidates failed for ", outcome_prefix,
              " — falling back to Gaussian.", call. = FALSE)
      use_gamma <- FALSE   # fall through to Gaussian below
    }
  }

  if (!use_gamma) {
    # ---- Gaussian LMM candidates (default) ----------------------------------
    for (time_var in c("timepoint", "time_num")) {
      fixed_formula <- if (time_var == "timepoint") GD_FIXED_FORMULA_FACTOR else GD_FIXED_FORMULA_NUMERIC
      for (re_name in names(re_forms)) {
        re_form <- re_forms[[re_name]]
        if (re_name == "trial_date" && !"trial_date" %in% names(data)) next

        fml <- tryCatch(
          stats::as.formula(sprintf("%s ~ %s + %s", response, fixed_formula, re_form)),
          error = function(e) NULL
        )
        if (is.null(fml)) next

        mname <- paste(outcome_prefix, paste0("time=", time_var),
                       paste0("RE=", re_name), sep = "__")

        fit <- tryCatch(lme4::lmer(fml, data = data, REML = FALSE), error = function(e) NULL)
        if (!is.null(fit)) {
          candidates[[mname]] <- fit
          meta_rows[[mname]] <- data.frame(
            model         = mname, time_var = time_var, re_name = re_name,
            fixed_formula = fixed_formula,
            singular      = lme4::isSingular(fit, tol = 1e-4),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(candidates) == 0L) {
    warning("All pre-specified candidate models failed for: ", outcome_prefix, call. = FALSE)
    return(NULL)
  }

  meta_df <- dplyr::bind_rows(meta_rows)
  aic_vals <- sapply(candidates, function(m)
    tryCatch(MuMIn::AICc(m), error = function(e) NA_real_))
  aic_tbl <- dplyr::tibble(model = names(candidates), AICc = aic_vals) %>%
    dplyr::left_join(meta_df, by = "model") %>%
    dplyr::arrange(AICc) %>%
    dplyr::mutate(delta = AICc - min(AICc, na.rm = TRUE))

  best_name <- aic_tbl$model[1L]

  if (!use_gamma) {
    # Prefer factor timepoint if within 2 AICc
    best_tp <- aic_tbl[aic_tbl$time_var == "timepoint", ]
    best_tn <- aic_tbl[aic_tbl$time_var == "time_num",  ]
    if (nrow(best_tp) > 0 && nrow(best_tn) > 0 &&
        aic_tbl$model[1L] == best_tn$model[1L] &&
        is.finite(best_tp$AICc[1L] - best_tn$AICc[1L]) &&
        (best_tp$AICc[1L] - best_tn$AICc[1L]) < 2) {
      best_name <- best_tp$model[1L]
    }
    # Prefer non-singular within delta 2
    best_aicc <- aic_tbl$AICc[aic_tbl$model == best_name]
    aic_tbl$singular <- vapply(aic_tbl$model, function(nm)
      lme4::isSingular(candidates[[nm]], tol = 1e-4), logical(1))
    near_ns <- aic_tbl[!aic_tbl$singular & aic_tbl$AICc <= best_aicc + 2, , drop = FALSE]
    if (nrow(near_ns) > 0 && isTRUE(aic_tbl$singular[aic_tbl$model == best_name]))
      best_name <- near_ns$model[which.min(near_ns$AICc)]
  }

  aic_tbl <- dplyr::mutate(aic_tbl, selected = (model == best_name))
  write.csv(aic_tbl, file.path(DIR_AICC, paste0("AICc_", outcome_prefix, ".csv")), row.names = FALSE)

  list(candidates = candidates, aic = aic_tbl,
       best = candidates[[best_name]], best_name = best_name)
}

# Beta GLMM for bounded [0, 1] outcomes (polarisation, zone_synchrony, ang_momentum).
# Models on proportion scale; estimates are on the logit scale internally.
fit_beta_gd <- function(data, response, outcome_prefix, re_forms = RE_FORMS) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    warning("glmmTMB not available — skipping Beta GLMM for ", outcome_prefix, call. = FALSE)
    return(NULL)
  }
  prop_col <- ".prop_beta"
  data[[prop_col]] <- pmax(pmin(data[[response]], 1 - 1e-4), 1e-4)
  data <- data[is.finite(data[[prop_col]]), , drop = FALSE]
  if (nrow(data) < 6L) return(NULL)

  candidates <- list()
  meta_rows  <- list()
  for (re_name in names(re_forms)) {
    if (re_name == "trial_date" && !"trial_date" %in% names(data)) next
    re_form <- re_forms[[re_name]]
    fml_str <- paste(prop_col, "~ Group * Treatment * timepoint +", re_form)
    fit <- tryCatch(
      glmmTMB::glmmTMB(stats::as.formula(fml_str), data = data,
                        family = glmmTMB::beta_family(link = "logit"),
                        control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1000))),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      mname <- paste(outcome_prefix, paste0("RE=", re_name), "beta", sep = "__")
      candidates[[mname]] <- fit
      meta_rows[[mname]]  <- data.frame(model = mname, re_name = re_name, stringsAsFactors = FALSE)
    }
  }
  if (length(candidates) == 0L) {
    warning("All Beta GLMM candidates failed for: ", outcome_prefix, call. = FALSE)
    return(NULL)
  }
  aics <- sapply(candidates, function(m) tryCatch(stats::AIC(m), error = function(e) NA_real_))
  best_name <- names(which.min(aics))
  list(candidates = candidates, aic = data.frame(model = names(aics), AIC = aics, stringsAsFactors = FALSE),
       best = candidates[[best_name]], best_name = best_name)
}

# Build the hierarchical fixed-effects candidate set (kept for legacy reference)
make_fixed_set <- function(time_var) {
  G <- "Group"; Tr <- "Treatment"; t <- time_var
  unique(c(
    sprintf("%s + %s + %s",            G, Tr, t),
    sprintf("%s*%s + %s",              G, Tr, t),
    sprintf("%s*%s + %s",              G, t,  Tr),
    sprintf("%s*%s + %s",              Tr, t, G),
    sprintf("%s*%s + %s*%s + %s",     G, Tr, G, t, Tr),
    sprintf("%s*%s + %s*%s + %s",     G, Tr, Tr, t, G),
    sprintf("%s*%s + %s*%s + %s",     G, t, Tr, t, G),
    sprintf("%s*%s + %s*%s + %s*%s",  G, Tr, G, t, Tr, t),
    sprintf("%s*%s*%s",               G, Tr, t)
  ))
}

fit_candidate_set_gd <- function(data, response, outcome_prefix) {
  # Try every combination of time parameterisation × random structure × fixed structure.
  # Returns NULL if no model converges.
  candidates <- list()
  meta_rows  <- list()

  for (time_var in c("timepoint","time_num")) {
    for (re_name in names(RE_FORMS)) {
      re_form <- RE_FORMS[[re_name]]
      # Skip trial_date RE if that column is absent
      if (re_name == "trial_date" && !"trial_date" %in% names(data)) next

      for (fx in make_fixed_set(time_var)) {
        fml <- tryCatch(
          stats::as.formula(sprintf("%s ~ %s + %s", response, fx, re_form)),
          error = function(e) NULL
        )
        if (is.null(fml)) next

        mname <- paste(outcome_prefix,
                       paste0("time=", time_var),
                       paste0("RE=", re_name),
                       gsub("[^A-Za-z0-9]+", "_", fx),
                       sep = "__")

        fit <- tryCatch(
          lme4::lmer(fml, data = data, REML = FALSE),
          error = function(e) NULL
        )

        if (!is.null(fit)) {
          candidates[[mname]] <- fit
          meta_rows[[mname]] <- data.frame(
            model         = mname,
            time_var      = time_var,
            re_name       = re_name,
            fixed_formula = fx,
            singular      = lme4::isSingular(fit, tol = 1e-4),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  if (length(candidates) == 0L) {
    warning("All candidate models failed for: ", outcome_prefix, call. = FALSE)
    return(NULL)
  }

  meta_df <- dplyr::bind_rows(meta_rows)

  aic_tbl <- tibble(model = names(candidates),
                    AICc  = sapply(candidates, MuMIn::AICc)) %>%
    left_join(meta_df, by = "model") %>%
    arrange(AICc) %>%
    mutate(delta = AICc - min(AICc))

  best_overall <- aic_tbl$model[1L]

  # Prefer timepoint-as-factor if it is within 2 AICc units of the best time_num model
  best_tp <- aic_tbl %>% filter(time_var == "timepoint") %>% slice(1L)
  best_tn <- aic_tbl %>% filter(time_var == "time_num")  %>% slice(1L)
  best_name <- best_overall

  if (nrow(best_tp) == 1L && nrow(best_tn) == 1L &&
      best_tn$model == best_overall &&
      is.finite(best_tp$AICc - best_tn$AICc) &&
      (best_tp$AICc - best_tn$AICc) < 2) {
    best_name <- best_tp$model
  }

  # Among models within 2 AICc units of best_name, prefer non-singular (delta_tol = 2)
  best_aicc <- aic_tbl$AICc[aic_tbl$model == best_name]
  aic_tbl$singular <- vapply(aic_tbl$model, function(nm) {
    m <- candidates[[nm]]
    if (inherits(m, "merMod")) lme4::isSingular(m, tol = 1e-4) else NA
  }, logical(1))

  near_nonsingular <- aic_tbl[
    !is.na(aic_tbl$singular) & !aic_tbl$singular &
    aic_tbl$AICc <= best_aicc + 2,
  ]
  if (nrow(near_nonsingular) > 0 && isTRUE(aic_tbl$singular[aic_tbl$model == best_name])) {
    best_name <- near_nonsingular$model[which.min(near_nonsingular$AICc)]
  }

  aic_tbl <- aic_tbl %>% mutate(selected = (model == best_name))
  write.csv(aic_tbl,
            file.path(DIR_AICC, paste0("AICc_", outcome_prefix, ".csv")),
            row.names = FALSE)

  list(candidates = candidates, aic = aic_tbl,
       best = candidates[[best_name]], best_name = best_name)
}

# Poisson / negative-binomial candidate fitter for count outcomes (e.g.
# n_col_mot_events).  Both families are tried via glmmTMB; AICc selects.
# An overdispersion diagnostic (Pearson chi-sq / residual df) is attached to
# the Poisson rows of the combined AICc table.
fit_candidate_set_gd_count <- function(data, response, outcome_prefix) {
  fit_one <- function(fml, family_obj, mname) {
    tryCatch(
      glmmTMB::glmmTMB(fml, data = data, family = family_obj),
      error = function(e) NULL
    )
  }

  candidates <- list()
  meta_rows  <- list()

  for (family_str in c("poisson", "nb")) {
    family_obj <- if (family_str == "poisson") stats::poisson() else glmmTMB::nbinom2()
    prefix_fam <- paste0(outcome_prefix, "_", family_str)

    for (time_var in c("timepoint", "time_num")) {
      for (re_name in names(RE_FORMS)) {
        re_form <- RE_FORMS[[re_name]]
        if (re_name == "trial_date" && !"trial_date" %in% names(data)) next

        for (fx in make_fixed_set(time_var)) {
          fml <- tryCatch(
            stats::as.formula(sprintf("%s ~ %s + %s", response, fx, re_form)),
            error = function(e) NULL
          )
          if (is.null(fml)) next

          mname <- paste(prefix_fam,
                         paste0("time=", time_var),
                         paste0("RE=", re_name),
                         gsub("[^A-Za-z0-9]+", "_", fx),
                         sep = "__")

          fit <- fit_one(fml, family_obj, mname)
          if (!is.null(fit)) {
            candidates[[mname]] <- fit
            meta_rows[[mname]] <- data.frame(
              model         = mname,
              family        = family_str,
              time_var      = time_var,
              re_name       = re_name,
              fixed_formula = fx,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  if (length(candidates) == 0L) {
    warning("All candidate count models failed for: ", outcome_prefix, call. = FALSE)
    return(NULL)
  }

  meta_df <- dplyr::bind_rows(meta_rows)

  aic_tbl <- tibble::tibble(model = names(candidates),
                     AICc  = sapply(candidates, MuMIn::AICc)) %>%
    dplyr::left_join(meta_df, by = "model") %>%
    dplyr::arrange(AICc) %>%
    dplyr::mutate(delta = AICc - min(AICc))

  # Overdispersion of best Poisson model
  best_pois_name <- aic_tbl %>% dplyr::filter(family == "poisson") %>% dplyr::slice(1L) %>%
    dplyr::pull(model)
  overdisp_ratio <- if (length(best_pois_name) == 1L) {
    tryCatch({
      m    <- candidates[[best_pois_name]]
      pr   <- sum(residuals(m, type = "pearson")^2, na.rm = TRUE)
      df_r <- df.residual(m)
      if (is.finite(df_r) && df_r > 0) round(pr / df_r, 3) else NA_real_
    }, error = function(e) NA_real_)
  } else NA_real_

  aic_tbl$overdisp_pearson_ratio <- ifelse(aic_tbl$family == "poisson", overdisp_ratio, NA_real_)

  best_name <- aic_tbl$model[1L]

  # Prefer timepoint-as-factor if within 2 AICc units of best time_num model
  best_tp <- aic_tbl %>% dplyr::filter(time_var == "timepoint") %>% dplyr::slice(1L)
  best_tn <- aic_tbl %>% dplyr::filter(time_var == "time_num")  %>% dplyr::slice(1L)
  if (nrow(best_tp) == 1L && nrow(best_tn) == 1L &&
      best_tn$model == best_name &&
      is.finite(best_tp$AICc - best_tn$AICc) &&
      (best_tp$AICc - best_tn$AICc) < 2) {
    best_name <- best_tp$model
  }

  aic_tbl <- aic_tbl %>% dplyr::mutate(selected = (model == best_name))
  write.csv(aic_tbl,
            file.path(DIR_AICC, paste0("AICc_", outcome_prefix, "_count.csv")),
            row.names = FALSE)

  list(candidates = candidates, aic = aic_tbl,
       best = candidates[[best_name]], best_name = best_name,
       overdisp_ratio = overdisp_ratio)
}

run_anova_gd <- function(model, prefix) {
  # Use Kenward-Roger F-tests when available; fall back to Wald chi-square.
  use_kr <- inherits(model, "lmerMod") &&
    requireNamespace("lmerTest", quietly = TRUE) &&
    requireNamespace("pbkrtest", quietly = TRUE)

  a <- if (use_kr) {
    tryCatch(
      car::Anova(lmerTest::as.lmerModLmerTest(model), type = 3, test.statistic = "F"),
      error = function(e) car::Anova(model, type = 3)
    )
  } else {
    car::Anova(model, type = 3)
  }

  adf <- as.data.frame(a)
  adf$term <- rownames(adf); rownames(adf) <- NULL

  pcol <- intersect(c("Pr(>Chisq)","Pr(>F)","p.value","p"), names(adf))[1L]
  if (!is.na(pcol)) {
    pv <- suppressWarnings(as.numeric(adf[[pcol]]))
    pv_bh <- p.adjust(pv, method = "BH")
    adf <- adf %>% mutate(
      p_value_numeric    = pv,
      p_value_BH         = pv_bh,
      significant        = is.finite(pv) & pv < ALPHA,
      significant_BH     = is.finite(pv_bh) & pv_bh < ALPHA,
      significance_label = dplyr::case_when(
        !is.finite(pv) ~ NA_character_,
        pv < 0.001 ~ "***", pv < 0.01 ~ "**", pv < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    )
  }

  # Attach attributes required by build_anova_caption_multiline() in pretty_labels.R
  stat_col_attr <- grep("^(F$|F value$|Chisq$|Chi\\.sq$|LR Chisq$)",
                        names(adf), value = TRUE)[1L]
  p_col_attr    <- grep("^Pr\\(", names(adf), value = TRUE)[1L]
  if (is.na(p_col_attr))
    p_col_attr <- intersect(c("p_value_numeric", "p.value", "p"), names(adf))[1L]
  df_cols_attr  <- grep("Df|df|NumDF|DenDF", names(adf), value = TRUE)
  df_cols_attr  <- setdiff(df_cols_attr,
                           c(p_col_attr, "p_value_numeric", "p_value_BH"))
  attr(adf, "stat_col") <- stat_col_attr
  attr(adf, "p_col")    <- p_col_attr
  attr(adf, "df_cols")  <- df_cols_attr

  write.csv(adf,
            file.path(DIR_ANOVA, paste0("ANOVA_typeIII_", prefix, ".csv")),
            row.names = FALSE)
  adf
}

run_posthoc_gd <- function(model, prefix, anova_df) {
  if (is.null(anova_df)) return(NULL)

  pcol <- intersect(c("Pr(>Chisq)","Pr(>F)","p.value","p","p_value_numeric"),
                    names(anova_df))[1L]
  if (is.na(pcol)) return(NULL)

  pv       <- suppressWarnings(as.numeric(anova_df[[pcol]]))
  ftxt     <- paste(deparse(stats::formula(model)), collapse = " ")
  uses_tp  <- grepl("\\btimepoint\\b", ftxt)
  allowed  <- c("Group","Treatment", if (uses_tp) "timepoint" else "time_num")

  # Helper: is every component of a term in the allowed set?
  is_eligible <- function(z) {
    parts <- strsplit(z, ":")[[1L]]
    length(parts) > 0L && all(parts %in% allowed)
  }

  # Full interaction term: always run CLD for this regardless of significance.
  full_int <- if (uses_tp) "timepoint:Treatment:Group" else "Treatment:Group"

  # All significant non-intercept terms that are eligible for post-hoc
  sig_terms <- anova_df$term[is.finite(pv) & pv < ALPHA & !is.na(anova_df$term)]
  sig_terms <- sig_terms[
    !grepl("Intercept", sig_terms, ignore.case = TRUE) &
    vapply(sig_terms, is_eligible, logical(1L))
  ]

  # Union: always include full interaction + all significant eligible terms
  terms_to_run <- unique(c(full_int, sig_terms))
  # Preserve ANOVA table order; fall back to the highest-order eligible term
  all_eligible <- anova_df$term[
    !is.na(anova_df$term) &
    !grepl("Intercept", anova_df$term, ignore.case = TRUE) &
    vapply(anova_df$term, is_eligible, logical(1L))
  ]
  terms_to_run <- intersect(c(full_int, all_eligible), terms_to_run)
  if (length(terms_to_run) == 0L) {
    terms_to_run <- all_eligible[which.max(
      vapply(all_eligible, function(z) stringr::str_count(z, ":"), integer(1L))
    )]
  }
  if (length(terms_to_run) == 0L) {
    write.csv(data.frame(outcome = prefix, posthoc_run = FALSE,
                         reason = "No eligible terms in ANOVA table",
                         stringsAsFactors = FALSE),
              file.path(DIR_POSTHOC, paste0("POSTHOC_status_", prefix, ".csv")),
              row.names = FALSE)
    return(NULL)
  }

  # Helper: get emmeans for one term
  get_emm <- function(term_str) {
    spec <- paste(strsplit(term_str, ":")[[1L]], collapse = " * ")
    tryCatch(
      if (uses_tp)
        emmeans::emmeans(model, stats::as.formula(paste0("~ ", spec)))
      else
        emmeans::emmeans(model, stats::as.formula(paste0("~ ", spec)),
                         at = list(time_num = ANALYSIS_TIMEPOINTS)),
      error = function(e) NULL
    )
  }

  # --- Primary pairwise contrasts: use the full interaction term ---
  emm_primary <- get_emm(full_int)
  if (is.null(emm_primary)) {
    # Fall back to highest-order eligible term
    highest <- terms_to_run[which.max(
      vapply(terms_to_run, function(z) stringr::str_count(z, ":"), integer(1L))
    )]
    emm_primary <- get_emm(highest)
    if (is.null(emm_primary)) return(NULL)
    full_int <- highest
  }

  pw_df  <- tryCatch(as.data.frame(emmeans::contrast(emm_primary, method = "pairwise",
                                                      adjust = "tukey")),
                     error = function(e) NULL)
  emm_df <- as.data.frame(emm_primary)

  # --- CLD tables: one per term in terms_to_run ---
  cld_tables <- list()
  for (trm in terms_to_run) {
    emm_trm <- if (identical(trm, full_int)) emm_primary else get_emm(trm)
    if (is.null(emm_trm)) next
    cld_trm <- tryCatch(
      as.data.frame(multcomp::cld(emm_trm, Letters = letters, adjust = "tukey")),
      error = function(e) NULL
    )
    if (!is.null(cld_trm)) {
      cld_tables[[trm]] <- cld_trm
      safe_trm <- gsub("[^A-Za-z0-9]+", "_", trm)
      write.csv(cld_trm,
                file.path(DIR_POSTHOC, paste0("CLD_", prefix, "_", safe_trm, ".csv")),
                row.names = FALSE)
    }
  }

  # Primary CLD (full interaction) for backward compatibility
  cld_df <- if (full_int %in% names(cld_tables)) cld_tables[[full_int]]
             else if (length(cld_tables) > 0L)   cld_tables[[1L]]
             else NULL

  if (!is.null(pw_df))
    write.csv(pw_df,  file.path(DIR_POSTHOC, paste0("POSTHOC_Tukey_", prefix, ".csv")),
              row.names = FALSE)
  if (!is.null(cld_df))
    write.csv(cld_df, file.path(DIR_POSTHOC, paste0("CLD_", prefix, ".csv")),
              row.names = FALSE)
  write.csv(emm_df,   file.path(DIR_POSTHOC, paste0("EMMEANS_", prefix, ".csv")),
            row.names = FALSE)

  list(chosen_term = full_int, emm = emm_df, pairs = pw_df,
       cld = cld_df, cld_tables = cld_tables)
}

# ---- run stats for each session-level indicator ----
SESSION_OUTCOMES <- list(
  list(col = "mean_polarization",       prefix = "GD_polarization",
       label = "Group polarization (0–1)",          count = FALSE, bounded = TRUE),
  list(col = "mean_ang_momentum",       prefix = "GD_ang_momentum",
       label = "Angular momentum (0–1)",            count = FALSE, bounded = TRUE),
  list(col = "mean_iid_cm",             prefix = "GD_IID_cm",
       label = "Mean inter-individual distance (cm)",   count = FALSE, bounded = FALSE),
  list(col = "zone_synchrony",          prefix = "GD_zone_synchrony",
       label = "Zone synchrony (fraction of frames)",   count = FALSE, bounded = TRUE),
  list(col = "mean_centroid_speed_cm_s",prefix = "GD_centroid_speed",
       label = "Group centroid speed (cm/s)",           count = FALSE, bounded = FALSE),
  list(col = "mean_hull_area_cm2",      prefix = "GD_hull_area",
       label = "Convex hull area (cm²)",           count = FALSE, bounded = FALSE)
  ## DISABLED: collective motion and leadership indicators removed.
  # list(col = "n_col_mot_events",        prefix = "GD_n_col_mot_events",
  #      label = "Collective motion events (count)",      count = TRUE,  bounded = FALSE),
)

# Gate RE_FORMS based on BL-match reliability tier (once at analysis start)
.bl_match_df_gd <- if (exists("bl_match_result", envir = .GlobalEnv, inherits = FALSE))
  get("bl_match_result", envir = .GlobalEnv) else NULL
RE_FORMS <- gate_re_forms_gd(.bl_match_df_gd)
message("STEP6 active RE_FORMS: ", paste(names(RE_FORMS), collapse = ", "))

gd_model_results <- list()

for (out in SESSION_OUTCOMES) {
  if (!(out$col %in% names(session_gd_summary))) {
    message("Skipping stats [", out$label, "] — column absent.")
    next
  }
  dat <- session_gd_summary %>% dplyr::filter(is.finite(.data[[out$col]])) %>% droplevels()
  if (nrow(dat) < 6L) {
    message("Skipping stats [", out$label, "] — fewer than 6 observations.")
    next
  }
  message("Fitting models: ", out$label)

  is_count   <- isTRUE(out$count)
  is_bounded <- isTRUE(out$bounded)

  # ---- Pre-modelling assumption diagnostics (continuous indicators only) ----
  .diag_gd <- NULL
  .response_col_gd <- out$col   # may be overwritten by a transformation

  if (!is_count && !is_bounded &&
      exists("run_diag_for_indicator", mode = "function")) {
    .gd_grouping <- tryCatch(
      interaction(dat$Group, dat$Treatment, drop = TRUE),
      error = function(e) rep("all", nrow(dat))
    )
    .diag_gd <- tryCatch(
      run_diag_for_indicator(
        y             = dat[[out$col]],
        data          = dat,
        grouping      = .gd_grouping,
        re_form_str   = "(1 | trial_id)",
        ind_nm        = out$col,
        analysis_unit = "session",
        try_gamma     = TRUE
      ),
      error = function(e) {
        warning("[DIAG] Diagnostics failed for ", out$col, ": ", e$message, call. = FALSE)
        NULL
      }
    )
    # Apply transformation to dat if needed
    if (!is.null(.diag_gd) && .diag_gd$transform != "identity" &&
        .diag_gd$transform != "gamma") {
      .response_col_gd <- .diag_gd$response_col
      dat[[.response_col_gd]] <- tryCatch(
        switch(.diag_gd$transform,
               log   = log(dat[[out$col]]),
               log1p = log1p(dat[[out$col]]),
               sqrt  = sqrt(dat[[out$col]]),
               dat[[out$col]]),
        error = function(e) dat[[out$col]]
      )
    } else if (!is.null(.diag_gd)) {
      .response_col_gd <- .diag_gd$response_col
    }
  }

  if (is_count) {
    dat[[out$col]] <- as.integer(round(dat[[out$col]]))
    fit_obj <- tryCatch(
      fit_candidate_set_gd_count(dat, out$col, out$prefix),
      error = function(e) NULL
    )
  } else if (is_bounded) {
    # Beta GLMM for polarisation / zone_synchrony / ang_momentum
    fit_obj <- tryCatch(
      fit_beta_gd(dat, out$col, out$prefix, re_forms = RE_FORMS),
      error = function(e) NULL
    )
  } else {
    # Pre-specified full 3-way + RE AICc selection
    .use_gamma_gd <- isTRUE(!is.null(.diag_gd) && .diag_gd$use_gamma)
    fit_obj <- tryCatch(
      fit_candidate_set_gd_prespecified(dat, .response_col_gd, out$prefix,
                                        re_forms = RE_FORMS,
                                        use_gamma = .use_gamma_gd),
      error = function(e) NULL
    )
  }

  if (is.null(fit_obj)) next
  anova_df <- tryCatch(run_anova_gd(fit_obj$best, out$prefix), error = function(e) NULL)
  post_obj <- tryCatch(run_posthoc_gd(fit_obj$best, out$prefix, anova_df), error = function(e) NULL)

  # DHARMa diagnostics
  dharma_summ_gd <- NULL
  if (requireNamespace("DHARMa", quietly = TRUE)) {
    sim_gd <- tryCatch(DHARMa::simulateResiduals(fit_obj$best, n = 500L, plot = FALSE),
                       error = function(e) NULL)
    if (!is.null(sim_gd)) {
      png_gd <- file.path(DIR_DIAG, paste0("dharma_", out$prefix, ".png"))
      grDevices::png(png_gd, width = 1600, height = 800, res = 150)
      try(DHARMa::plot(sim_gd, title = out$label), silent = TRUE)
      grDevices::dev.off()
      dharma_summ_gd <- tryCatch(
        data.frame(outcome = out$col,
                   KS_p         = DHARMa::testUniformity(sim_gd, plot = FALSE)$p.value,
                   dispersion_p = DHARMa::testDispersion(sim_gd, plot = FALSE)$p.value,
                   stringsAsFactors = FALSE),
        error = function(e) NULL
      )
    }
  }

  # R²m / R²c
  r2_gd <- tryCatch(MuMIn::r.squaredGLMM(fit_obj$best), error = function(e) NULL)

  # Build the diagnostics row (fallback for bounded/count that bypassed run_diag_for_indicator)
  .diag_row_gd <- if (!is.null(.diag_gd)) {
    .diag_gd$row
  } else if (is_bounded) {
    if (exists("diag_row_beta", mode = "function")) {
      clipped_any <- any(dat[[out$col]] %in% c(0, 1), na.rm = TRUE)
      diag_row_beta(nrow(dat), "session", out$col, clipped = clipped_any)
    } else NULL
  } else if (is_count) {
    if (exists("diag_row_count", mode = "function")) {
      diag_row_count(nrow(dat), "session", out$col,
                     overdisp = fit_obj$overdisp_ratio)
    } else NULL
  } else NULL

  gd_model_results[[out$prefix]] <- list(
    col = out$col, label = out$label, bounded = is_bounded,
    fit = fit_obj, anova = anova_df, post = post_obj,
    overdisp_ratio = if (is_count) fit_obj$overdisp_ratio else NULL,
    dharma = dharma_summ_gd,
    r2m = if (!is.null(r2_gd) && "R2m" %in% colnames(r2_gd)) r2_gd[1L, "R2m"] else NA_real_,
    r2c = if (!is.null(r2_gd) && "R2c" %in% colnames(r2_gd)) r2_gd[1L, "R2c"] else NA_real_,
    diagnostics = .diag_row_gd
  )
}

# ---- NND: fish-level model (with fish_id RE, gated on BL tier) ----
nnd_model_results <- NULL
if (!is.null(nnd_model_data) && nrow(nnd_model_data) >= 6L) {
  message("Fitting models: Nearest-neighbour distance (NND)")
  nnd_fit   <- tryCatch(
    fit_candidate_set_gd_prespecified(nnd_model_data, "mean_nnd_cm", "GD_NND_cm",
                                       re_forms = RE_FORMS),
    error = function(e) NULL
  )
  nnd_anova <- if (!is.null(nnd_fit))
    tryCatch(run_anova_gd(nnd_fit$best, "GD_NND_cm"), error = function(e) NULL)
  nnd_post  <- if (!is.null(nnd_anova))
    tryCatch(run_posthoc_gd(nnd_fit$best, "GD_NND_cm", nnd_anova), error = function(e) NULL)
  nnd_model_results <- list(
    col = "mean_nnd_cm",
    label = "Nearest-neighbour distance (cm, fish-level)",
    fit = nnd_fit, anova = nnd_anova, post = nnd_post
  )
}


# =============================================================================
# 8B) INDIVIDUAL-LEVEL ANALYSIS
# =============================================================================
# Builds a fish x session x timepoint dataset combining:
#   - Individual movement indicators from fish_activity_summary (STEP 4/5)
#   - NND per fish (mean distance to nearest neighbour, from section 5A)
#   - Per-fish IID (mean distance to all other fish, from section 5B)
# LMMs use the same pre-specified fixed structure (Group*Treatment*timepoint)
# with fish-level random effects where idmatcher qualifies.

ind_model_results <- list()

fish_act_sum <- get_global("fish_activity_summary")
if (is.null(fish_act_sum)) {
  # Try the most recent STEP4_output folder
  .step4_root <- file.path(dirname(normalizePath(PIPELINE_DATA_DIR, mustWork = FALSE)),
                            "STEP4_output")
  .step4_dirs <- if (dir.exists(.step4_root))
    list.dirs(.step4_root, recursive = FALSE, full.names = TRUE) else character(0)
  if (length(.step4_dirs) > 0L) {
    .fas_path <- file.path(sort(.step4_dirs, decreasing = TRUE)[1L],
                           "fish_activity_summary.csv")
    if (file.exists(.fas_path)) {
      message("[STEP 6] Loading fish_activity_summary from disk:\n  ", .fas_path)
      fish_act_sum <- data.table::fread(.fas_path, data.table = FALSE)
    }
  }
}

if (!is.null(fish_act_sum) && nrow(fish_act_sum) > 0L) {

  # Individual movement columns — subset to those present in fish_act_sum
  .ind_move_cols <- c(
    "distance_m", "distance_BL",
    "speed_cm_s_mean", "speed_cm_s_median", "speed_cm_s_max",
    "speed_BL_s_mean", "speed_BL_s_median", "speed_BL_s_max",
    "accel_cm_s2_mean", "accel_cm_s2_median", "accel_cm_s2_max",
    "n_zone_transitions",
    "time_moving_TS0.5_s", "time_moving_TS1_s", "time_spent_zone_flow_s"
  )

  # Resolve fish ID column (fish_activity_summary may use different names)
  .fid_col <- pick_col(fish_act_sum, c("fish_id", "fish_id_matched_trial", "fish_id_raw"))
  if (is.na(.fid_col)) {
    warning("[STEP 6] Individual analysis: fish_id column not found — skipped.", call. = FALSE)
    fish_act_sum <- NULL
  }
}

if (!is.null(fish_act_sum) && nrow(fish_act_sum) > 0L) {

  # Resolve timepoint column name (may be 'segment_id' or 'timepoint')
  .tp_col_fas <- pick_col(fish_act_sum, c("timepoint", "segment_id"))

  # Standardize column names using base-R to avoid dplyr no-op rename edge cases
  .fas_work <- fish_act_sum
  if (!is.na(.fid_col) && .fid_col != "fish_id")
    names(.fas_work)[names(.fas_work) == .fid_col]    <- "fish_id"
  if (!is.na(.tp_col_fas) && .tp_col_fas != "timepoint")
    names(.fas_work)[names(.fas_work) == .tp_col_fas] <- "timepoint"

  ind_base <- .fas_work %>%
    mutate(trial_id  = as.integer(as.character(trial_id)),
           timepoint = as.integer(as.character(timepoint))) %>%
    filter(timepoint %in% ANALYSIS_TIMEPOINTS) %>%
    { if (!isTRUE(USE_ALL_GROUPS)) filter(., Group %in% GROUPS_RESTRICTED) else . } %>%
    mutate(
      timepoint = factor(timepoint, levels = ANALYSIS_TIMEPOINTS),
      Group     = factor(Group,     levels = c("No choice","Exercise choice",
                                               "Naloxone","Morphine")),
      Treatment = factor(Treatment, levels = c("Control","Exercise")),
      time_num  = as.numeric(as.character(timepoint))
    ) %>%
    droplevels()

  # Join NND per fish (trial_id x timepoint x fish_id)
  if (!is.null(nnd_fish_tp) && nrow(nnd_fish_tp) > 0L &&
      "fish_id" %in% names(nnd_fish_tp)) {
    ind_base <- ind_base %>%
      mutate(.tp_int = as.integer(as.character(timepoint))) %>%
      left_join(
        nnd_fish_tp %>%
          select(trial_id, timepoint, fish_id, mean_nnd_cm) %>%
          mutate(trial_id = as.integer(trial_id),
                 fish_id  = as.character(fish_id),
                 .tp_int  = as.integer(timepoint)) %>%
          select(trial_id, fish_id, mean_nnd_cm, .tp_int),
        by = c("trial_id", ".tp_int", "fish_id")
      ) %>%
      select(-.tp_int)
  }

  # Join per-fish IID (trial_id x timepoint x fish_id)
  if (!is.null(iid_fish_tp) && nrow(iid_fish_tp) > 0L &&
      "fish_id" %in% names(iid_fish_tp)) {
    ind_base <- ind_base %>%
      mutate(.tp_int = as.integer(as.character(timepoint))) %>%
      left_join(
        iid_fish_tp %>%
          select(trial_id, timepoint, fish_id,
                 mean_iid_cm_fish = mean_iid_cm) %>%
          mutate(trial_id = as.integer(trial_id),
                 fish_id  = as.character(fish_id),
                 .tp_int  = as.integer(timepoint)) %>%
          select(trial_id, fish_id, mean_iid_cm_fish, .tp_int),
        by = c("trial_id", ".tp_int", "fish_id")
      ) %>%
      select(-.tp_int)
  }

  # Ensure fish_id_matched_trial exists for fish-level RE candidate
  if (!"fish_id_matched_trial" %in% names(ind_base))
    ind_base$fish_id_matched_trial <- paste(ind_base$trial_id, ind_base$fish_id, sep = "_")

  write.csv(ind_base,
            file.path(DIR_DATA, paste0("individual_analysis_dataset_",
                                       TIMEPOINT_LABEL, "_", GROUP_LABEL, ".csv")),
            row.names = FALSE)
  message("Individual-level dataset: ", nrow(ind_base), " rows x ",
          ncol(ind_base), " columns")

  # All individual indicators to model (movement + group-level per fish)
  .ind_all_cols <- unique(c(
    intersect(.ind_move_cols, names(ind_base)),
    if ("mean_nnd_cm"      %in% names(ind_base)) "mean_nnd_cm",
    if ("mean_iid_cm_fish" %in% names(ind_base)) "mean_iid_cm_fish"
  ))

  IND_OUTCOMES <- lapply(.ind_all_cols, function(col) {
    list(col     = col,
         prefix  = paste0("IND_", gsub("[^A-Za-z0-9]+", "_", col)),
         label   = pretty_indicator_name(col),   # manuscript-friendly label
         bounded = col %in% c("time_spent_zone_flow_s"),
         count   = col %in% c("n_zone_transitions"))
  })

  for (out in IND_OUTCOMES) {
    if (!(out$col %in% names(ind_base))) next
    dat <- ind_base %>% dplyr::filter(is.finite(.data[[out$col]])) %>% droplevels()
    if (nrow(dat) < 6L) {
      message("Individual analysis — skipping [", out$label, "] — fewer than 6 rows.")
      next
    }
    message("Individual analysis — fitting: ", out$label)

    # ---- Pre-modelling assumption diagnostics (continuous indicators only) ----
    .diag_ind <- NULL
    .response_col_ind <- out$col

    if (!isTRUE(out$bounded) && !isTRUE(out$count) &&
        exists("run_diag_for_indicator", mode = "function")) {
      .ind_grouping <- tryCatch(
        interaction(dat$Group, dat$Treatment, drop = TRUE),
        error = function(e) rep("all", nrow(dat))
      )
      .diag_ind <- tryCatch(
        run_diag_for_indicator(
          y             = dat[[out$col]],
          data          = dat,
          grouping      = .ind_grouping,
          re_form_str   = "(1 | trial_id)",
          ind_nm        = out$col,
          analysis_unit = "fish",
          try_gamma     = TRUE
        ),
        error = function(e) {
          warning("[DIAG] Diagnostics failed for ", out$col, ": ", e$message, call. = FALSE)
          NULL
        }
      )
      if (!is.null(.diag_ind) && .diag_ind$transform != "identity" &&
          .diag_ind$transform != "gamma") {
        .response_col_ind <- .diag_ind$response_col
        dat[[.response_col_ind]] <- tryCatch(
          switch(.diag_ind$transform,
                 log   = log(dat[[out$col]]),
                 log1p = log1p(dat[[out$col]]),
                 sqrt  = sqrt(dat[[out$col]]),
                 dat[[out$col]]),
          error = function(e) dat[[out$col]]
        )
      } else if (!is.null(.diag_ind)) {
        .response_col_ind <- .diag_ind$response_col
      }
    }

    .use_gamma_ind <- isTRUE(!is.null(.diag_ind) && .diag_ind$use_gamma)

    fit_obj <- tryCatch({
      if (isTRUE(out$bounded)) {
        fit_beta_gd(dat, out$col, out$prefix, re_forms = RE_FORMS)
      } else if (isTRUE(out$count)) {
        dat[[out$col]] <- as.integer(round(dat[[out$col]]))
        fit_candidate_set_gd_count(dat, out$col, out$prefix)
      } else {
        fit_candidate_set_gd_prespecified(dat, .response_col_ind, out$prefix,
                                          re_forms = RE_FORMS,
                                          use_gamma = .use_gamma_ind)
      }
    }, error = function(e) NULL)

    if (is.null(fit_obj)) next
    anova_ind  <- tryCatch(run_anova_gd(fit_obj$best, out$prefix), error = function(e) NULL)
    post_ind   <- tryCatch(run_posthoc_gd(fit_obj$best, out$prefix, anova_ind),
                           error = function(e) NULL)

    r2_ind <- tryCatch(MuMIn::r.squaredGLMM(fit_obj$best), error = function(e) NULL)

    # Build the diagnostics row
    .diag_row_ind <- if (!is.null(.diag_ind)) {
      .diag_ind$row
    } else if (isTRUE(out$bounded)) {
      if (exists("diag_row_beta", mode = "function")) {
        clipped_any_ind <- any(dat[[out$col]] %in% c(0, 1), na.rm = TRUE)
        diag_row_beta(nrow(dat), "fish", out$col, clipped = clipped_any_ind)
      } else NULL
    } else if (isTRUE(out$count)) {
      if (exists("diag_row_count", mode = "function")) {
        diag_row_count(nrow(dat), "fish", out$col,
                       overdisp = fit_obj$overdisp_ratio)
      } else NULL
    } else NULL

    ind_model_results[[out$prefix]] <- list(
      col = out$col, label = out$label,
      bounded = isTRUE(out$bounded),
      fit = fit_obj, anova = anova_ind, post = post_ind,
      r2m = if (!is.null(r2_ind) && "R2m" %in% colnames(r2_ind)) r2_ind[1L, "R2m"] else NA_real_,
      r2c = if (!is.null(r2_ind) && "R2c" %in% colnames(r2_ind)) r2_ind[1L, "R2c"] else NA_real_,
      diagnostics = .diag_row_ind
    )
  }

  message("Individual-level analysis complete. ",
          length(ind_model_results), " indicators fitted.")

} else {
  message("Individual-level analysis skipped — fish_activity_summary not available.")
}


# =============================================================================
# 9) VISUALISATION
# =============================================================================
# Two plot types per indicator, matching the style of Step 4:
#   Scatter  — x = Treatment, points coloured by Group, mean±SEM crossbars,
#              faceted by timepoint.
#   Trend    — x = timepoint, one line per Group×Treatment, error bars = SEM.

GROUP_COLORS <- c(
  "No choice"       = "#56B4E9",
  "Exercise choice" = "#D55E00",
  "Morphine"        = "#8DBB1F",
  "Naloxone"        = "#C084F0"
)
GROUP_ORDER     <- c("No choice","Exercise choice","Naloxone","Morphine")
TREATMENT_ORDER <- c("Control","Exercise")

BASE_THEME <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.background  = element_blank(),
    axis.line.x  = element_line(color = "black", linewidth = 1.05),
    axis.line.y  = element_line(color = "black", linewidth = 1.05),
    axis.ticks   = element_line(color = "black", linewidth = 0.8),
    axis.text.x  = element_text(size = 16, face = "bold",  color = "black"),
    axis.text.y  = element_text(size = 12,                 color = "black"),
    axis.title.x = element_text(size = 16, face = "bold",  color = "black"),
    axis.title.y = element_text(size = 16, face = "bold",  color = "black"),
    strip.text   = element_text(size = 14, face = "bold",  color = "black"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text  = element_text(size = 12),
    panel.spacing = grid::unit(1.2, "lines"),
    plot.margin   = ggplot2::margin(t = 8, r = 10, b = 18, l = 10),
    plot.caption  = element_text(size = 10.5, hjust = 0,
                                  margin = margin(t = 18), color = "black",
                                  lineheight = 1.05),
    plot.caption.position = "plot"
  )

# Build a multi-line ANOVA caption string from a run_anova_gd() result data frame.
# Uses build_anova_caption_multiline() from pretty_labels.R (always sourced first).
# Returns "" gracefully when anova_df is NULL / lacks required attributes.
build_gd_anova_caption <- function(anova_df) {
  if (is.null(anova_df) || !is.data.frame(anova_df) || nrow(anova_df) == 0L) return("")
  # Ensure required attributes exist (set by run_anova_gd(); add fallback here)
  if (is.null(attr(anova_df, "stat_col"))) {
    attr(anova_df, "stat_col") <- grep("^(F$|F value$|Chisq$)", names(anova_df),
                                        value = TRUE)[1L]
  }
  if (is.null(attr(anova_df, "p_col"))) {
    attr(anova_df, "p_col") <- grep("^Pr\\(", names(anova_df), value = TRUE)[1L]
    if (is.na(attr(anova_df, "p_col")))
      attr(anova_df, "p_col") <- intersect(c("p_value_numeric","p.value","p"),
                                            names(anova_df))[1L]
  }
  if (is.null(attr(anova_df, "df_cols"))) {
    attr(anova_df, "df_cols") <- grep("Df|df|NumDF|DenDF", names(anova_df), value = TRUE)
  }
  tryCatch(build_anova_caption_multiline(anova_df, caption_format = "plain"),
           error = function(e) "")
}

# Build a sample-size caption string: "Interval 2: n = 12, Interval 3: n = 14"
build_n_caption_gd <- function(df) {
  ns <- df %>%
    mutate(.tp = as.integer(as.character(timepoint))) %>%
    group_by(.tp) %>%
    summarise(n = dplyr::n(), .groups = "drop") %>%
    arrange(.tp)
  paste(paste0("Interval ", ns$.tp, ": n = ", ns$n), collapse = ", ")
}

make_scatter_gd <- function(df, y_col, y_label, title = y_label,
                            anova_df = NULL, cld_df = NULL) {
  summ <- df %>%
    group_by(Treatment, Group, timepoint) %>%
    summarise(m = safe_mean(.data[[y_col]]),
              s = safe_sem(.data[[y_col]]),
              .groups = "drop") %>%
    filter(is.finite(m))

  pos_jd <- position_jitterdodge(jitter.width = 0.14, jitter.height = 0,
                                  dodge.width  = 0.75)
  pos_d  <- position_dodge(width = 0.75)

  df_p <- df %>%
    filter(is.finite(.data[[y_col]])) %>%
    mutate(Group     = factor(Group,     levels = GROUP_ORDER),
           Treatment = factor(Treatment, levels = TREATMENT_ORDER))

  # Build caption: sample sizes + ANOVA results
  n_cap    <- build_n_caption_gd(df_p)
  anova_cap <- build_gd_anova_caption(anova_df)
  caption_txt <- if (nzchar(anova_cap))
    paste0(n_cap, "\n", anova_cap)
  else
    n_cap

  p <- ggplot(df_p, aes(x = Treatment, y = .data[[y_col]],
                         color = Group, group = Group)) +
    geom_point(position = pos_jd, size = 3.6, alpha = 0.75) +
    geom_errorbar(data = summ,
                  aes(x = Treatment, y = m, ymin = m - s, ymax = m + s,
                      color = Group, group = Group),
                  position = pos_d, width = 0.18, linewidth = 0.7,
                  inherit.aes = FALSE) +
    geom_crossbar(data = summ,
                  aes(x = Treatment, y = m, ymin = m, ymax = m,
                      color = Group, group = Group),
                  position = pos_d, width = 0.55, linewidth = 0.9,
                  inherit.aes = FALSE) +
    scale_color_manual(values = GROUP_COLORS,
                       limits = names(GROUP_COLORS), drop = FALSE) +
    facet_wrap(~ timepoint, nrow = 1L,
               labeller = as_labeller(function(x) paste0("Interval ", x))) +
    labs(x = "Treatment", y = y_label, title = title, caption = caption_txt) +
    BASE_THEME

  # --- CLD letter annotation (suppressed when all groups share one letter) ---
  if (!is.null(cld_df) && exists("add_cld_to_scatter", mode = "function")) {
    tp_vals <- as.character(unique(df_p$timepoint))
    p <- tryCatch(
      add_cld_to_scatter(p,
                         y_data      = df_p[[y_col]],
                         cld_df      = cld_df,
                         summ_df     = summ,
                         x_col       = "Treatment",
                         grp_col     = "Group",
                         tp_col      = "timepoint",
                         tp_values   = tp_vals,
                         dodge_width = 0.75,
                         font_size   = 4.5),
      error = function(e) p    # silent fallback: return plot without CLD
    )
  }

  p
}

make_trend_gd <- function(df, y_col, y_label,
                           title = paste0(y_label, " — trend"), anova_df = NULL) {
  summ <- df %>%
    mutate(tp_num = as.numeric(as.character(timepoint))) %>%
    group_by(Treatment, Group, tp_num) %>%
    summarise(m = safe_mean(.data[[y_col]]),
              s = safe_sem(.data[[y_col]]),
              .groups = "drop") %>%
    filter(is.finite(m)) %>%
    mutate(Group     = factor(Group,     levels = GROUP_ORDER),
           Treatment = factor(Treatment, levels = TREATMENT_ORDER))

  pos_d <- position_dodge(width = 0.10)

  # Build caption: sample sizes + ANOVA results
  n_cap     <- build_n_caption_gd(df %>% filter(is.finite(.data[[y_col]])))
  anova_cap <- build_gd_anova_caption(anova_df)
  caption_txt <- if (nzchar(anova_cap))
    paste0(n_cap, "\n", anova_cap)
  else
    n_cap

  ggplot(summ, aes(x = tp_num, y = m,
                   color = Group, shape = Treatment, linetype = Treatment,
                   group = interaction(Group, Treatment))) +
    geom_line(linewidth = 1.2, position = pos_d) +
    geom_point(size = 3.8, position = pos_d) +
    geom_errorbar(aes(ymin = m - s, ymax = m + s),
                  width = 0.10, linewidth = 0.7, position = pos_d) +
    scale_color_manual(values = GROUP_COLORS,
                       limits = names(GROUP_COLORS), drop = FALSE) +
    scale_linetype_manual(values = c("Control" = "dashed", "Exercise" = "solid")) +
    scale_shape_manual(values   = c("Control" = 17L,       "Exercise" = 16L)) +
    scale_x_continuous(
      breaks = ANALYSIS_TIMEPOINTS,
      labels = pretty_interval_name(as.character(ANALYSIS_TIMEPOINTS))
    ) +
    labs(x = "Interval", y = y_label, title = title, caption = caption_txt) +
    BASE_THEME
}

save_to_pptx <- function(plots_named, pptx_path,
                          slide_w = 13.33, slide_h = 7.5) {
  ppt <- officer::read_pptx()
  for (nm in names(plots_named)) {
    ppt <- officer::add_slide(ppt, layout = "Blank", master = "Office Theme")
    ppt <- officer::ph_with(
      ppt,
      value    = rvg::dml(ggobj = plots_named[[nm]]),
      location = officer::ph_location(
        left   = 0.2,  top    = 0.2,
        width  = slide_w - 0.4,
        height = slide_h - 0.4
      )
    )
  }
  print(ppt, target = pptx_path)
  invisible(pptx_path)
}

# ---- all indicators to plot ----
# Labels match the canonical GOT_INDICATOR_LABELS (defined in pretty_labels.R).
PLOT_INDICATORS_SESSION <- list(
  list(col = "mean_polarization",        label = "Group polarization (0–1)"),
  list(col = "mean_ang_momentum",        label = "Angular momentum (0–1)"),
  list(col = "mean_iid_cm",              label = "Mean inter-individual distance (cm)"),
  list(col = "zone_synchrony",           label = "Zone synchrony (fraction of frames)"),
  list(col = "mean_centroid_speed_cm_s", label = "Group centroid speed (cm/s)"),
  list(col = "mean_hull_area_cm2",       label = "Convex hull area (cm²)")
  ## DISABLED indicators removed from list — they are absent from session_gd_summary:
  # list(col = "n_col_mot_events",  label = "Collective motion events (count)"),
  # list(col = "mean_event_dur_s",  label = "Collective motion event duration (s)"),
  # list(col = "leadership_index",  label = "Leadership index")
)

PLOT_NND <- list(col = "mean_nnd_cm", label = "Mean nearest-neighbour distance (cm)")

# OPT-0: build branch list dynamically to match upstream filter switches.
branch_specs <- if (isTRUE(EXCLUDE_PHARMA_GROUPS)) {
  message("OPT-0 (STEP6): EXCLUDE_PHARMA_GROUPS=TRUE -> Exe_choice branch only")
  list(list(branch_name = "Exe_choice", group_filter = c("Exercise choice", "No choice")))
} else {
  list(
    list(branch_name = "All_groups", group_filter = NULL),
    list(branch_name = "Exe_choice", group_filter = c("Exercise choice", "No choice"))
  )
}
subbranch_specs <- if (isTRUE(USE_TP1)) {
  list(
    list(sb_name = "All_timepoints", tp_filter = NULL),
    list(sb_name = "Timepoints_23",  tp_filter = c(2L, 3L)),
    list(sb_name = "Timepoint1",     tp_filter = c(1L))
  )
} else {
  message("OPT-0 (STEP6): USE_TP1=FALSE -> dropping Timepoint1 subbranch")
  list(
    list(sb_name = "All_timepoints", tp_filter = NULL),
    list(sb_name = "Timepoints_23",  tp_filter = c(2L, 3L))
  )
}

# Pre-build ANOVA lookup: col name → anova_df (with attributes set by run_anova_gd)
gd_anova_lookup <- setNames(
  lapply(gd_model_results, function(r) r$anova),
  vapply(gd_model_results, function(r) if (!is.null(r$col)) r$col else "", character(1L))
)
# Post-hoc lookup keyed by col (mirrors gd_anova_lookup)
gd_post_lookup <- setNames(
  lapply(gd_model_results, function(r) r$post),
  vapply(gd_model_results, function(r) if (!is.null(r$col)) r$col else "", character(1L))
)
# Add NND
if (!is.null(nnd_model_results) && !is.null(nnd_model_results$anova))
  gd_anova_lookup[["mean_nnd_cm"]] <- nnd_model_results$anova
if (!is.null(nnd_model_results) && !is.null(nnd_model_results$post))
  gd_post_lookup[["mean_nnd_cm"]] <- nnd_model_results$post

for (br in branch_specs) {
  for (sb in subbranch_specs) {

    sub_dir <- file.path(DIR_GRAPHS, br$branch_name, sb$sb_name)
    dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)

    png_dir <- file.path(sub_dir, "interval_level_png")
    dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

    scatter_plots <- list()
    trend_plots   <- list()

    # --- session-level indicators ---
    for (ind in PLOT_INDICATORS_SESSION) {
      if (!(ind$col %in% names(session_gd_summary))) next

      df_sub <- session_gd_summary %>% filter(is.finite(.data[[ind$col]]))
      if (!is.null(br$group_filter))
        df_sub <- df_sub %>% filter(Group %in% br$group_filter)
      if (!is.null(sb$tp_filter))
        df_sub <- df_sub %>%
          filter(as.integer(as.character(timepoint)) %in% sb$tp_filter)
      if (nrow(df_sub) < 2L) next
      df_sub <- df_sub %>% droplevels()

      anova_ind <- gd_anova_lookup[[ind$col]]

      # Build CLD annotation from model post-hoc result
      .cld_sess <- tryCatch({
        if (exists("extract_cld_for_plot", mode = "function") &&
            exists("gd_post_lookup", mode = "list") &&
            !is.null(gd_post_lookup[[ind$col]])) {
          tp_vals_sess <- as.character(unique(df_sub$timepoint))
          extract_cld_for_plot(gd_post_lookup[[ind$col]],
                               tp_values = tp_vals_sess)
        } else NULL
      }, error = function(e) NULL)

      tryCatch({
        scatter_plots[[ind$col]] <- make_scatter_gd(df_sub, ind$col, ind$label,
                                                     anova_df = anova_ind,
                                                     cld_df   = .cld_sess)
      }, error = function(e) NULL)

      has_multiple_tp <- length(unique(as.integer(as.character(df_sub$timepoint)))) >= 2L
      if (has_multiple_tp) {
        tryCatch({
          trend_plots[[ind$col]] <- make_trend_gd(df_sub, ind$col, ind$label,
                                                   anova_df = anova_ind)
        }, error = function(e) NULL)
      }
    }

    # --- NND (fish-level) ---
    if (!is.null(nnd_model_data) && nrow(nnd_model_data) > 0L) {
      df_nnd <- nnd_model_data
      if (!is.null(br$group_filter))
        df_nnd <- df_nnd %>% filter(Group %in% br$group_filter)
      if (!is.null(sb$tp_filter))
        df_nnd <- df_nnd %>%
          filter(as.integer(as.character(timepoint)) %in% sb$tp_filter)

      if (nrow(df_nnd) >= 2L) {
        df_nnd <- df_nnd %>% droplevels()
        nnd_anova <- gd_anova_lookup[["mean_nnd_cm"]]
        .cld_nnd <- tryCatch({
          if (exists("extract_cld_for_plot", mode = "function") &&
              exists("gd_post_lookup", mode = "list") &&
              !is.null(gd_post_lookup[["mean_nnd_cm"]])) {
            tp_vals_nnd <- as.character(unique(df_nnd$timepoint))
            extract_cld_for_plot(gd_post_lookup[["mean_nnd_cm"]],
                                 tp_values = tp_vals_nnd)
          } else NULL
        }, error = function(e) NULL)
        tryCatch({
          scatter_plots[["mean_nnd_cm"]] <- make_scatter_gd(
            df_nnd, "mean_nnd_cm", PLOT_NND$label, anova_df = nnd_anova,
            cld_df = .cld_nnd)
        }, error = function(e) NULL)

        if (length(unique(as.integer(as.character(df_nnd$timepoint)))) >= 2L) {
          tryCatch({
            trend_plots[["mean_nnd_cm"]] <- make_trend_gd(
              df_nnd, "mean_nnd_cm", PLOT_NND$label, anova_df = nnd_anova)
          }, error = function(e) NULL)
        }
      }
    }

    if (length(scatter_plots) > 0L) {
      sc_path <- file.path(
        sub_dir,
        paste0("GROUP_DYNAMICS_scatter_", br$branch_name, "_", sb$sb_name, ".pptx")
      )
      tryCatch(save_to_pptx(scatter_plots, sc_path),
               error = function(e)
                 warning("Scatter PPT failed [", br$branch_name, "/",
                         sb$sb_name, "]: ", e$message, call. = FALSE))

      # Save scatter PNGs
      for (nm in names(scatter_plots)) {
        safe_indic <- gsub("[^A-Za-z0-9_]", "_", nm)
        png_file   <- file.path(png_dir, paste0(safe_indic, "_scatter.png"))
        tryCatch(
          ggplot2::ggsave(png_file, plot = scatter_plots[[nm]], width = 10, height = 7,
                          dpi = 150, bg = "white"),
          error = function(e) warning("PNG save failed: ", conditionMessage(e), call. = FALSE)
        )
      }
    }

    if (length(trend_plots) > 0L) {
      tr_path <- file.path(
        sub_dir,
        paste0("GROUP_DYNAMICS_trend_", br$branch_name, "_", sb$sb_name, ".pptx")
      )
      tryCatch(save_to_pptx(trend_plots, tr_path),
               error = function(e)
                 warning("Trend PPT failed [", br$branch_name, "/",
                         sb$sb_name, "]: ", e$message, call. = FALSE))

      # Save trend PNGs
      for (nm in names(trend_plots)) {
        safe_indic <- gsub("[^A-Za-z0-9_]", "_", nm)
        png_file   <- file.path(png_dir, paste0(safe_indic, "_trend.png"))
        tryCatch(
          ggplot2::ggsave(png_file, plot = trend_plots[[nm]], width = 10, height = 7,
                          dpi = 150, bg = "white"),
          error = function(e) warning("PNG save failed: ", conditionMessage(e), call. = FALSE)
        )
      }
    }

    message("Graphs (session-level): ", br$branch_name, " / ", sb$sb_name,
            "  [", length(scatter_plots), " scatter + ",
            length(trend_plots), " trend]")

    # ---- Individual-level graphs (one per indicator in ind_model_results) ----
    if (length(ind_model_results) > 0L &&
        exists("ind_base", inherits = FALSE) &&
        !is.null(ind_base) && nrow(ind_base) > 0L) {

      ind_png_dir <- file.path(sub_dir, "individual_level_png")
      dir.create(ind_png_dir, showWarnings = FALSE, recursive = TRUE)

      ind_anova_lookup <- setNames(
        lapply(ind_model_results, function(r) r$anova),
        vapply(ind_model_results, function(r) if (!is.null(r$col)) r$col else "",
               character(1L))
      )

      n_ind_sc <- 0L; n_ind_tr <- 0L
      for (ind_nm in names(ind_model_results)) {
        res_ind <- ind_model_results[[ind_nm]]
        col_nm  <- res_ind$col
        lbl     <- res_ind$label

        if (!col_nm %in% names(ind_base)) next

        df_ind <- ind_base %>% filter(is.finite(.data[[col_nm]]))
        if (!is.null(br$group_filter))
          df_ind <- df_ind %>% filter(Group %in% br$group_filter)
        if (!is.null(sb$tp_filter))
          df_ind <- df_ind %>%
            filter(as.integer(as.character(timepoint)) %in% sb$tp_filter)
        if (nrow(df_ind) < 2L) next
        df_ind <- df_ind %>% droplevels()

        anova_ind2 <- ind_anova_lookup[[col_nm]]
        safe_nm    <- gsub("[^A-Za-z0-9_]", "_", col_nm)

        # Build CLD annotation for individual-level indicator
        .cld_ind <- tryCatch({
          if (exists("extract_cld_for_plot", mode = "function") &&
              !is.null(ind_model_results[[ind_nm]]$post)) {
            tp_vals_ind <- as.character(unique(df_ind$timepoint))
            extract_cld_for_plot(ind_model_results[[ind_nm]]$post,
                                 tp_values = tp_vals_ind)
          } else NULL
        }, error = function(e) NULL)

        tryCatch({
          p_sc <- make_scatter_gd(df_ind, col_nm, lbl, anova_df = anova_ind2,
                                   cld_df = .cld_ind)
          ggplot2::ggsave(
            file.path(ind_png_dir, paste0(safe_nm, "_scatter.png")),
            plot = p_sc, width = 10, height = 7, dpi = 150, bg = "white")
          n_ind_sc <- n_ind_sc + 1L
        }, error = function(e) NULL)

        has_multi_tp <- length(unique(as.integer(as.character(df_ind$timepoint)))) >= 2L
        if (has_multi_tp) {
          tryCatch({
            p_tr <- make_trend_gd(df_ind, col_nm, lbl, anova_df = anova_ind2)
            ggplot2::ggsave(
              file.path(ind_png_dir, paste0(safe_nm, "_trend.png")),
              plot = p_tr, width = 10, height = 7, dpi = 150, bg = "white")
            n_ind_tr <- n_ind_tr + 1L
          }, error = function(e) NULL)
        }
      }
      message("Graphs (individual-level): ", br$branch_name, " / ", sb$sb_name,
              "  [", n_ind_sc, " scatter + ", n_ind_tr, " trend]")
    }
  }
}


# =============================================================================
# 9b) MARKDOWN REPORT WRITER
# =============================================================================

# Convert a data.frame to a Markdown table string.
# Numeric columns are formatted to 3 decimal places.
# Row names are included as a "Term" column if they are non-trivial.
df_to_md_table <- function(df) {
  if (!is.data.frame(df) || nrow(df) == 0L) return("*No data.*")
  df2 <- as.data.frame(lapply(df, function(col) {
    if (is.numeric(col)) formatC(col, digits = 3L, format = "f") else as.character(col)
  }), stringsAsFactors = FALSE, check.names = FALSE)
  rn <- rownames(df)
  if (!is.null(rn) && any(rn != as.character(seq_len(nrow(df))))) {
    df2 <- cbind(Term = rn, df2)
  }
  header <- paste0("| ", paste(names(df2), collapse = " | "), " |")
  sep    <- paste0("| ", paste(rep("---", ncol(df2)), collapse = " | "), " |")
  rows   <- apply(df2, 1L, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

# Write a Markdown report for STEP6 group-dynamics and individual-level results.
# all_results_gd : named list — c(gd_model_results, nnd entry); each item has
#   $label, $anova (data.frame), $post (list with $chosen_term, $pairs, $cld),
#   $fit (list with $best_name, $aic)
# all_results_ind: named list — ind_model_results (same structure as above)
# n_sessions     : integer, number of sessions analysed
# out_dir        : base output directory (OUT_DIR)
write_report_markdown_step6 <- function(all_results_gd, all_results_ind, n_sessions, out_dir) {
  md_dir <- file.path(out_dir, "reports_md")
  dir.create(md_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(md_dir,
    paste0("report_group_dynamics_",
           if (exists("TIMEPOINT_LABEL")) TIMEPOINT_LABEL else "TP",
           "_",
           if (exists("GROUP_LABEL"))     GROUP_LABEL     else "groups",
           ".md"))

  lines <- c(
    "# GOT Group Dynamics -- Statistical Report", "",
    paste0("**N sessions:** ", n_sessions, "  "),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "", "---", ""
  )

  # Helper: render one result item (gd or ind) to Markdown lines
  render_result_item <- function(res, key) {
    out_lines <- c(paste0("## ", res$label), "")

    if (!is.null(res$fit) && !is.null(res$fit$best_name)) {
      out_lines <- c(out_lines,
        paste0("**Selected model:** `", res$fit$best_name, "`"), "")
    }

    # ANOVA table
    at <- res$anova
    if (!is.null(at) && is.data.frame(at) && nrow(at) > 0L) {
      # Prettify: apply pretty term names, drop dot-columns, rename headers
      at_pr <- at
      if ("term" %in% names(at_pr) &&
          exists("pretty_term_name", mode = "function")) {
        at_pr$term <- vapply(as.character(at_pr$term),
                             pretty_term_name, character(1L))
      }
      # Drop internal helper columns (names starting with ".")
      dot_c <- grep("^\\.", names(at_pr), value = TRUE)
      if (length(dot_c) > 0L)
        at_pr <- at_pr[, !names(at_pr) %in% dot_c, drop = FALSE]
      rownames(at_pr) <- NULL
      out_lines <- c(out_lines, "### Type III ANOVA", "", df_to_md_table(at_pr), "")
    } else {
      out_lines <- c(out_lines, "*ANOVA table not available.*", "")
    }

    # Post-hoc
    po <- res$post
    if (!is.null(po)) {
      pretty_chosen <- if (exists("pretty_term_name", mode = "function"))
        pretty_term_name(po$chosen_term) else po$chosen_term
      out_lines <- c(out_lines,
        "### Post-hoc comparisons", "",
        paste0("**Post-hoc term:** ", pretty_chosen), "")
      if (!is.null(po$pairs) && is.data.frame(po$pairs) && nrow(po$pairs) > 0L) {
        out_lines <- c(out_lines,
          "#### Pairwise contrasts (Tukey)", "",
          df_to_md_table(po$pairs), "")
      }
      # Iterate all CLD tables (one per term); fall back to single $cld if needed
      if (!is.null(po$cld_tables) && length(po$cld_tables) > 0L) {
        for (cld_term in names(po$cld_tables)) {
          cld_tbl <- po$cld_tables[[cld_term]]
          if (!is.data.frame(cld_tbl) || nrow(cld_tbl) == 0L) next
          p_cld_term <- if (exists("pretty_term_name", mode = "function"))
            pretty_term_name(cld_term) else gsub(":", " × ", cld_term)
          out_lines <- c(out_lines,
            paste0("#### CLD: ", p_cld_term), "",
            df_to_md_table(cld_tbl), "")
        }
      } else if (!is.null(po$cld) && is.data.frame(po$cld) && nrow(po$cld) > 0L) {
        out_lines <- c(out_lines,
          "#### Compact letter display", "",
          df_to_md_table(po$cld), "")
      }
    } else {
      out_lines <- c(out_lines,
        "*Post-hoc not run (no significant ANOVA terms at alpha = 0.05).*", "")
    }

    out_lines
  }

  # Group-dynamics section
  lines <- c(lines, "# Group-Dynamics Indicators", "")
  valid_gd <- names(all_results_gd)[!vapply(all_results_gd, is.null, logical(1L))]
  for (key in valid_gd) {
    res <- all_results_gd[[key]]
    if (is.null(res$fit)) next
    lines <- c(lines, render_result_item(res, key), "---", "")
  }

  # Individual-level section
  if (length(all_results_ind) > 0L) {
    lines <- c(lines, "# Individual-Level Indicators", "")
    valid_ind <- names(all_results_ind)[!vapply(all_results_ind, is.null, logical(1L))]
    for (key in valid_ind) {
      res <- all_results_ind[[key]]
      if (is.null(res$fit)) next
      lines <- c(lines, render_result_item(res, key), "---", "")
    }
  }

  tryCatch(
    writeLines(lines, out_file, useBytes = FALSE),
    error = function(e) warning("STEP6 Markdown report write failed: ",
                                conditionMessage(e), call. = FALSE)
  )
  message("STEP6 Markdown report written: ", out_file)
  invisible(out_file)
}

# =============================================================================
# 10) WORD REPORT
# =============================================================================
ps_a4 <- officer::prop_section(
  page_size    = officer::page_size(orient = "portrait", width = 8.27, height = 11.69),
  page_margins = officer::page_mar(top = 0.9843, bottom = 0.9843,
                                   left = 0.9843, right = 0.9843)
)

PAGE_W    <- 8.27
MARGIN    <- 0.9843
USABLE_W  <- PAGE_W - 2 * MARGIN

doc <- officer::read_docx()
doc <- officer::body_set_default_section(doc, ps_a4)

add_p   <- function(doc, text, style = "Normal")
  officer::body_add_par(doc, text, style = style)

fmt_sig3 <- function(x) {
  # Format to exactly 3 decimal places (max 3 digits after the decimal point).
  x_num <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x_num), "NA",
         formatC(round(x_num, 3L), format = "f", digits = 3L))
}

fmt_df_sig3 <- function(df) {
  out <- df
  for (nm in names(out)) {
    if (is.numeric(out[[nm]])) out[[nm]] <- fmt_sig3(out[[nm]])
  }
  out
}

fmt_anova_journal <- function(df) {
  d    <- as.data.frame(df)
  # Format every p-value column: 3 decimal places; use "< 0.0001" only below that threshold.
  p_fmt_col <- function(pv) {
    ifelse(is.na(pv), "NA",
           ifelse(pv < 0.0001, "< 0.0001",
                  formatC(pv, format = "f", digits = 3L)))
  }
  for (pcol in intersect(c("Pr(>Chisq)","Pr(>F)","p_value_numeric","p.value","p",
                            "p_value_BH"), names(d))) {
    pv        <- suppressWarnings(as.numeric(d[[pcol]]))
    d[[pcol]] <- p_fmt_col(pv)
  }
  # Format test statistics to 3 decimal places
  for (cc in intersect(c("Chisq","F value","F"), names(d))) {
    x       <- suppressWarnings(as.numeric(d[[cc]]))
    d[[cc]] <- ifelse(is.na(x), "NA", fmt_sig3(x))
  }
  d
}

# Rename ANOVA data frame columns to manuscript-friendly headings and apply
# pretty_term_name() to the "term" column values.  Also strips internal
# dot-prefixed columns (e.g. ".stat_col", ".p_col") before rendering.
# Handles duplicate column names by numbering duplicates.
prettify_anova_cols_gd <- function(df) {
  if (!is.data.frame(df)) return(df)

  # Apply pretty term names
  if ("term" %in% names(df) && exists("pretty_term_name", mode = "function")) {
    df$term <- vapply(as.character(df$term), pretty_term_name, character(1L))
  }

  # Drop internal dot-prefixed helper columns
  dot_cols <- grep("^\\.", names(df), value = TRUE)
  if (length(dot_cols) > 0L)
    df <- df[, !names(df) %in% dot_cols, drop = FALSE]

  # Format p-value columns to exactly 3 decimal places (as character) BEFORE renaming,
  # so flextable cannot display them with floating-point rounding artefacts.
  p_val_cols <- intersect(c("p_value_numeric", "p_value_BH",
                             "Pr(>F)", "Pr(>Chisq)", "p.value", "p"),
                          names(df))
  for (pc in p_val_cols) {
    raw <- suppressWarnings(as.numeric(df[[pc]]))
    df[[pc]] <- ifelse(is.na(raw), NA_character_,
                       ifelse(raw <= 0.0001, "< 0.0001",
                              formatC(raw, format = "f", digits = 3L)))
  }
  # Format numeric test-statistic columns to 3 decimal places
  stat_cols <- intersect(c("F value", "Chisq", "F"), names(df))
  for (sc in stat_cols) {
    raw <- suppressWarnings(as.numeric(df[[sc]]))
    df[[sc]] <- ifelse(is.na(raw), NA_character_,
                       formatC(round(raw, 3L), format = "f", digits = 3L))
  }

  # Canonical column-name map (exact match first)
  col_map <- c(
    "term"            = "Term",
    "Pr(>Chisq)"      = "p-value",
    "Pr(>F)"          = "p-value",
    "Chisq"           = "χ²",  # χ²
    "F value"         = "F",
    "F"               = "F",
    "NumDF"           = "df (num)",
    "DenDF"           = "df (den)",
    "Df"              = "df",
    "df"              = "df",
    "p_value_numeric" = "p-value",
    "p_value_BH"      = "p (BH-adj)",
    "significant"     = "p < 0.05",
    "significant_BH"  = "p (BH) < 0.05",
    "AICc"            = "AICc",
    "AIC"             = "AIC"
  )

  new_names <- names(df)
  for (i in seq_along(new_names)) {
    nm <- new_names[i]
    if (nm %in% names(col_map)) {
      new_names[i] <- col_map[[nm]]
    } else {
      # Fallback: replace underscores with spaces + title-case
      new_names[i] <- gsub("_", " ", nm)
      new_names[i] <- paste0(toupper(substr(new_names[i], 1, 1)),
                             substr(new_names[i], 2, nchar(new_names[i])))
    }
  }

  # Handle duplicates by appending a numeric suffix
  dup_idx <- which(duplicated(new_names))
  if (length(dup_idx) > 0L) {
    for (idx in dup_idx) {
      base_name <- new_names[idx]
      # Count how many times this name appears before this index
      prev_count <- sum(new_names[seq_len(idx - 1L)] == base_name)
      new_names[idx] <- paste0(base_name, " (", prev_count + 1L, ")")
    }
  }

  names(df) <- new_names
  df
}

make_ft <- function(df, font_size = 9L) {
  df2 <- fmt_df_sig3(df)
  ft  <- flextable::flextable(df2) %>%
    flextable::theme_booktabs() %>%
    flextable::font(fontname = "Calibri", part = "all") %>%
    flextable::padding(padding = 2L, part = "all") %>%
    flextable::fontsize(size = font_size, part = "body") %>%
    flextable::fontsize(size = max(7L, font_size - 1L), part = "header") %>%
    flextable::bold(part = "header") %>%
    flextable::align(align = "left", part = "body") %>%
    flextable::align(align = "center", part = "header")
  n   <- ncol(df2)
  flextable::width(ft, width = rep(USABLE_W / n, n)) %>%
    flextable::set_table_properties(layout = "fixed", width = 1L)
}

add_tbl <- function(doc, title, df, font_size = 9L, max_rows = NULL) {
  doc <- add_p(doc, title, style = "heading 3")
  if (!is.null(max_rows) && nrow(df) > max_rows) {
    doc <- add_p(doc, paste0("First ", max_rows, " rows; full table in CSV."))
    df  <- df[seq_len(max_rows), , drop = FALSE]
  }
  flextable::body_add_flextable(doc, make_ft(df, font_size))
}

# =============================================================================
# Trial-aggregated 2-way branch for STEP6 (pseudoreplication-free)
# Must run BEFORE section 10b so trial_agg_gd_results is available for Section 4.
# =============================================================================
tryCatch({
  if (exists("aggregate_trial_level", mode = "function") &&
      exists("session_gd_summary") && is.data.frame(session_gd_summary)) {

    trial_agg_gd_dir <- file.path(OUT_DIR, "trial_aggregated_2way")
    dir.create(trial_agg_gd_dir, recursive = TRUE, showWarnings = FALSE)

    gd_outcome_cols <- vapply(SESSION_OUTCOMES, `[[`, character(1), "col")
    gd_outcome_cols <- gd_outcome_cols[gd_outcome_cols %in% names(session_gd_summary)]

    dat_gd_agg_full <- session_gd_summary
    if (!isTRUE(USE_ALL_GROUPS)) {
      dat_gd_agg_full <- dat_gd_agg_full[
        dat_gd_agg_full$Group %in% GROUPS_RESTRICTED, , drop = FALSE]
    }

    # Convert bounded outcomes to proportion scale for aggregation
    for (out in SESSION_OUTCOMES) {
      if (isTRUE(out$bounded) && out$col %in% names(dat_gd_agg_full)) {
        dat_gd_agg_full[[paste0(out$col, "_prop")]] <-
          pmax(pmin(dat_gd_agg_full[[out$col]], 1 - 1e-4), 1e-4)
      }
    }

    dat_trial_gd <- tryCatch(
      aggregate_trial_level(dat_gd_agg_full, gd_outcome_cols, tps_to_avg = c(2L, 3L)),
      error = function(e) NULL
    )
    if (!is.null(dat_trial_gd) && nrow(dat_trial_gd) >= 5) {
      trial_agg_gd_results <- lapply(gd_outcome_cols, function(ind) {
        if (!ind %in% names(dat_trial_gd)) return(NULL)
        d <- dat_trial_gd[is.finite(dat_trial_gd[[ind]]), , drop = FALSE]
        if (nrow(d) < 5) return(NULL)
        fit <- tryCatch(stats::lm(stats::as.formula(paste(ind, "~ Treatment * Group")), data = d),
                        error = function(e) NULL)
        if (is.null(fit)) return(NULL)
        anova_tab <- tryCatch(as.data.frame(car::Anova(fit, type = 3)), error = function(e) data.frame())
        # Post-hoc assumption diagnostics for aggregated lm
        .agg_diag_row <- if (exists("run_agg_diag", mode = "function")) {
          tryCatch({
            .agg_grp <- tryCatch(
              interaction(d$Group, d$Treatment, drop = TRUE),
              error = function(e) rep("all", nrow(d))
            )
            run_agg_diag(fit, .agg_grp, ind, n_obs = nrow(d))
          }, error = function(e) NULL)
        } else NULL
        list(indicator = ind, model = fit, anova_table = anova_tab,
             n_trials = nrow(d), model_type = class(fit)[1L],
             diagnostics = .agg_diag_row)
      })
      names(trial_agg_gd_results) <- gd_outcome_cols

      anova_sum_gd <- do.call(rbind, lapply(names(trial_agg_gd_results), function(ind) {
        r <- trial_agg_gd_results[[ind]]
        if (is.null(r) || !is.data.frame(r$anova_table) || nrow(r$anova_table) == 0) return(NULL)
        at <- r$anova_table
        p_col <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), names(at))
        cbind(data.frame(indicator = ind, term = rownames(at), n_trials = r$n_trials,
                         stringsAsFactors = FALSE),
              at[, union(p_col, names(at)), drop = FALSE])
      }))
      if (!is.null(anova_sum_gd) && nrow(anova_sum_gd) > 0) {
        p_col_s <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), names(anova_sum_gd))
        if (length(p_col_s) > 0)
          anova_sum_gd$p_BH <- p.adjust(anova_sum_gd[[p_col_s[1]]], method = "BH")
        write.csv(anova_sum_gd,
                  file.path(trial_agg_gd_dir, "trial_aggregated_anova_summary_STEP6.csv"),
                  row.names = FALSE)
      }
      write.csv(dat_trial_gd, file.path(trial_agg_gd_dir, "trial_aggregated_data_STEP6.csv"),
                row.names = FALSE)
      message("STEP6 trial-aggregated branch written to: ", trial_agg_gd_dir)
    }
  }
}, error = function(e) warning("STEP6 trial-aggregated branch failed: ", conditionMessage(e), call. = FALSE))

# =============================================================================
# 10a) REPORT HELPER FUNCTIONS
# =============================================================================

# Returns the display label for a result list item (works for both GD and IND
# result formats — tries $label, then $col, then the key name).
.res_label <- function(res, key = NULL) {
  if (!is.null(res$label) && nzchar(res$label)) return(res$label)
  col <- res$col %||% key
  if (!is.null(col) && exists("pretty_indicator_name", mode = "function"))
    return(pretty_indicator_name(col))
  gsub("_", " ", as.character(col))
}
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Extract R²m (marginal, fixed effects only) and R²c (conditional, fixed +
# random) from any fitted model object.  For plain lm, R²m = R²c = R².
compute_r2_from_model <- function(fit) {
  if (is.null(fit)) return(list(r2m = NA_real_, r2c = NA_real_))
  if (inherits(fit, c("lmerMod", "lmerModLmerTest", "glmerMod", "glmmTMB"))) {
    r2 <- tryCatch(MuMIn::r.squaredGLMM(fit), error = function(e) NULL)
    if (!is.null(r2) && all(c("R2m", "R2c") %in% colnames(r2)))
      return(list(r2m = r2[1L, "R2m"], r2c = r2[1L, "R2c"]))
  }
  if (inherits(fit, "lm")) {
    rs <- tryCatch(summary(fit)$r.squared, error = function(e) NA_real_)
    return(list(r2m = rs, r2c = rs))
  }
  list(r2m = NA_real_, r2c = NA_real_)
}

# Format R²m / R²c scalars into a 1-row display data frame for add_tbl().
format_r2_df <- function(r2m, r2c, model_name = NULL) {
  fmt <- function(x) if (is.finite(x)) formatC(x, format = "f", digits = 3L) else "NA"
  data.frame(
    "Model"                        = if (!is.null(model_name)) model_name else "",
    "R²m (marginal)"         = fmt(r2m),
    "R²c (conditional)"      = fmt(r2c),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Build the standardised AICc preview data frame from a $fit$aic table.
.aic_preview_df <- function(aic_df) {
  data.frame(
    Candidate       = aic_df$model,
    "Time variable" = if ("time_var"      %in% names(aic_df)) aic_df$time_var      else NA_character_,
    "Random effects"= if ("re_name"       %in% names(aic_df)) aic_df$re_name       else NA_character_,
    "Fixed formula" = if ("fixed_formula" %in% names(aic_df)) aic_df$fixed_formula else NA_character_,
    Singular        = if ("singular" %in% names(aic_df))
                        ifelse(is.na(aic_df$singular), "NA",
                               ifelse(aic_df$singular, "YES", "NO"))
                      else NA_character_,
    AICc  = if ("AICc" %in% names(aic_df)) aic_df$AICc
            else if ("AIC" %in% names(aic_df)) aic_df$AIC else NA_real_,
    Delta = if ("delta" %in% names(aic_df)) aic_df$delta else NA_real_,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# Render one interval-level result item (GD or individual) into the doc.
# Works for any result produced by run_posthoc_gd() / select_best_gd_model().
.render_result_block <- function(doc, res) {
  # ---- Assumption-diagnostics table (placed first, before AICc) -----------
  if (!is.null(res$diagnostics) && is.data.frame(res$diagnostics) &&
      nrow(res$diagnostics) > 0L)
    doc <- add_tbl(doc, "Assumption diagnostics",
                   res$diagnostics, font_size = 8L)

  if (!is.null(res$fit$aic))
    doc <- add_tbl(doc, "AICc candidate model ranking",
                   .aic_preview_df(res$fit$aic), font_size = 8L, max_rows = 15L)

  doc <- add_p(doc, paste0("Selected model: ", res$fit$best_name))

  # R² table — uses stored r2m/r2c when available, falls back to model object
  {
    r2m_val <- if (!is.null(res$r2m) && is.finite(res$r2m)) res$r2m else NA_real_
    r2c_val <- if (!is.null(res$r2c) && is.finite(res$r2c)) res$r2c else NA_real_
    if (!is.finite(r2m_val) || !is.finite(r2c_val)) {
      .r2_fallback <- compute_r2_from_model(res$fit$best)
      if (!is.finite(r2m_val)) r2m_val <- .r2_fallback$r2m
      if (!is.finite(r2c_val)) r2c_val <- .r2_fallback$r2c
    }
    if (is.finite(r2m_val) || is.finite(r2c_val))
      doc <- add_tbl(doc, "R² (Nakagawa & Schielzeth 2013)",
                     format_r2_df(r2m_val, r2c_val, res$fit$best_name),
                     font_size = 9L)
  }

  if (!is.null(res$anova))
    doc <- add_tbl(doc, "Type III Wald ANOVA",
                   prettify_anova_cols_gd(fmt_anova_journal(res$anova)),
                   font_size = 9L)

  if (is.null(res$post)) {
    doc <- add_p(doc,
      "Post-hoc comparisons: not run (no significant ANOVA terms at α = 0.05).")
  } else {
    doc <- add_p(doc, paste0("Post-hoc comparison term: ",
                             pretty_term_name(res$post$chosen_term)))
    if (!is.null(res$post$pairs) && nrow(res$post$pairs) > 0L)
      doc <- add_tbl(doc, "Pairwise contrasts (Tukey-adjusted)",
                     res$post$pairs, font_size = 8L, max_rows = 30L)

    if (!is.null(res$post$cld_tables) && length(res$post$cld_tables) > 0L) {
      for (cld_term in names(res$post$cld_tables)) {
        cld_tbl <- res$post$cld_tables[[cld_term]]
        if (!is.data.frame(cld_tbl) || nrow(cld_tbl) == 0L) next
        doc <- add_tbl(doc,
          paste0("Compact letter display — ", pretty_term_name(cld_term)),
          cld_tbl, font_size = 9L)
        # Remapped CLD table for graph-annotation transparency
        if (exists("build_remapped_report_tbl", mode = "function")) {
          .remap_tbl <- tryCatch(build_remapped_report_tbl(cld_tbl), error = function(e) NULL)
          if (!is.null(.remap_tbl) && nrow(.remap_tbl) > 0L)
            doc <- add_tbl(doc, "Remapped CLD letters used in graphs",
                           .remap_tbl, font_size = 9L)
        }
      }
    } else if (!is.null(res$post$cld) && nrow(res$post$cld) > 0L) {
      doc <- add_tbl(doc, "Compact letter display", res$post$cld, font_size = 9L)
      if (exists("build_remapped_report_tbl", mode = "function")) {
        .remap_tbl <- tryCatch(build_remapped_report_tbl(res$post$cld), error = function(e) NULL)
        if (!is.null(.remap_tbl) && nrow(.remap_tbl) > 0L)
          doc <- add_tbl(doc, "Remapped CLD letters used in graphs",
                         .remap_tbl, font_size = 9L)
      }
    }
  }
  doc
}

# Render one trial-aggregated result item (STEP5 or STEP6 aggregated branch).
# Renders: assumption diagnostics, model summary, R², ANOVA, post-hoc contrasts, CLD.
.render_agg_block <- function(doc, res, ind_key) {
  fit    <- res$model
  n_tri  <- if (!is.null(res$n_trials))   as.integer(res$n_trials)   else NA_integer_
  m_type <- if (!is.null(res$model_type)) res$model_type else
            if (!is.null(fit)) class(fit)[1L] else "lm"
  fmla   <- if (!is.null(fit))
    tryCatch(paste(deparse(stats::formula(fit)), collapse = " "), error = function(e) "unknown")
  else "unknown"
  doc <- add_p(doc, paste0("N trials: ", n_tri,
                            "   Model: ", m_type,
                            "   Formula: ", fmla))

  # ---- Assumption-diagnostics table ----------------------------------------
  if (!is.null(res$diagnostics) && is.data.frame(res$diagnostics) &&
      nrow(res$diagnostics) > 0L)
    doc <- add_tbl(doc, "Assumption diagnostics",
                   res$diagnostics, font_size = 8L)

  # R² table
  r2_vals <- compute_r2_from_model(fit)
  if (is.finite(r2_vals$r2m) || is.finite(r2_vals$r2c))
    doc <- add_tbl(doc, "R² (Nakagawa & Schielzeth 2013)",
                   format_r2_df(r2_vals$r2m, r2_vals$r2c, m_type),
                   font_size = 9L)

  # ANOVA table
  at <- res$anova_table
  if (is.data.frame(at) && nrow(at) > 0L) {
    at_pr <- cbind(Term = rownames(at), as.data.frame(at))
    rownames(at_pr) <- NULL
    at_pr$Term <- vapply(as.character(at_pr$Term), pretty_term_name, character(1L))
    # Add p_value_numeric so prettify_anova_cols_gd can format it
    pcol_at <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), names(at_pr))
    if (length(pcol_at) > 0L && !"p_value_numeric" %in% names(at_pr))
      at_pr$p_value_numeric <- suppressWarnings(as.numeric(at_pr[[pcol_at[1L]]]))
    doc <- add_tbl(doc, "Type III ANOVA", prettify_anova_cols_gd(at_pr), font_size = 9L)
  } else {
    doc <- add_p(doc, "ANOVA table not available for this indicator.")
    return(doc)
  }

  # Post-hoc contrasts and CLD — build an anova_df with $term and $p_value_numeric
  # columns that run_posthoc_gd() expects, then call it directly.
  if (!is.null(fit)) {
    anova_ph <- cbind(term = rownames(at), as.data.frame(at))
    rownames(anova_ph) <- NULL
    pcol_ph <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), names(anova_ph))
    if (length(pcol_ph) > 0L)
      anova_ph$p_value_numeric <- suppressWarnings(as.numeric(anova_ph[[pcol_ph[1L]]]))
    safe_key <- gsub("[^A-Za-z0-9]+", "_", ind_key)
    post <- tryCatch(
      run_posthoc_gd(fit, paste0("AGG_", safe_key), anova_ph),
      error = function(e) {
        warning("Agg post-hoc failed (", ind_key, "): ", e$message, call. = FALSE)
        NULL
      }
    )
    if (is.null(post)) {
      doc <- add_p(doc,
        "Post-hoc: not run (no significant terms at α = 0.05 or model error).")
    } else {
      doc <- add_p(doc, paste0("Post-hoc term: ", pretty_term_name(post$chosen_term)))
      if (!is.null(post$pairs) && nrow(post$pairs) > 0L)
        doc <- add_tbl(doc, "Pairwise contrasts (Tukey-adjusted)",
                       post$pairs, font_size = 8L, max_rows = 30L)
      if (!is.null(post$cld_tables) && length(post$cld_tables) > 0L) {
        for (.cld_trm in names(post$cld_tables)) {
          .cld_tbl <- post$cld_tables[[.cld_trm]]
          if (!is.data.frame(.cld_tbl) || nrow(.cld_tbl) == 0L) next
          doc <- add_tbl(doc,
            paste0("Compact letter display — ", pretty_term_name(.cld_trm)),
            .cld_tbl, font_size = 9L)
          if (exists("build_remapped_report_tbl", mode = "function")) {
            .remap_tbl <- tryCatch(build_remapped_report_tbl(.cld_tbl), error = function(e) NULL)
            if (!is.null(.remap_tbl) && nrow(.remap_tbl) > 0L)
              doc <- add_tbl(doc, "Remapped CLD letters used in graphs",
                             .remap_tbl, font_size = 9L)
          }
        }
      } else if (!is.null(post$cld) && nrow(post$cld) > 0L) {
        doc <- add_tbl(doc, "Compact letter display", post$cld, font_size = 9L)
        if (exists("build_remapped_report_tbl", mode = "function")) {
          .remap_tbl <- tryCatch(build_remapped_report_tbl(post$cld), error = function(e) NULL)
          if (!is.null(.remap_tbl) && nrow(.remap_tbl) > 0L)
            doc <- add_tbl(doc, "Remapped CLD letters used in graphs",
                           .remap_tbl, font_size = 9L)
        }
      }
    }
  }
  doc
}

# Add a numbered section heading followed by bullet-point summary.
.add_section_header <- function(doc, num, title, bullets) {
  doc <- officer::body_add_break(doc)
  doc <- add_p(doc, paste0(num, ". ", title), style = "heading 1")
  for (b in bullets) doc <- add_p(doc, paste0("•  ", b))
  doc
}

# Load STEP5 trial-aggregated individual results for Section 3.
# Checks GlobalEnv first (set when STEP5 runs in the same session),
# then falls back to the stats_cache saved by STEP5.
.load_trial_agg_ind <- function() {
  if (exists("trial_agg_results_step5", envir = .GlobalEnv, inherits = FALSE)) {
    res <- get("trial_agg_results_step5", envir = .GlobalEnv)
    message("Combined report: STEP5 trial-agg results loaded from GlobalEnv (",
            length(res), " indicators)")
    return(res)
  }
  if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
    rds_path <- file.path(PIPELINE_DATA_DIR, "dataset", "stats_cache",
                          "trial_agg_results_step5.rds")
    if (file.exists(rds_path)) {
      res <- tryCatch(readRDS(rds_path), error = function(e) NULL)
      if (!is.null(res)) {
        message("Combined report: STEP5 trial-agg results loaded from stats_cache (",
                length(res), " indicators)")
        return(res)
      }
    }
  }
  message("Combined report: STEP5 trial-agg results not found — Section 3 will be empty.")
  NULL
}

# =============================================================================
# 10b) BUILD COMBINED FOUR-SECTION WORD REPORT
# =============================================================================

# =============================================================================
# Report cover page
# =============================================================================
doc <- add_p(doc,
             paste0("GOT Behavioural Statistics — Combined Report"),
             style = "heading 1")
doc <- add_p(doc, paste0(TIMEPOINT_LABEL, " | ", GROUP_LABEL,
                          "  —  Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M")))
doc <- add_p(doc, paste0(
  "Analysed intervals: ",
  paste(pretty_interval_name(as.character(ANALYSIS_TIMEPOINTS)), collapse = ", "), ".  ",
  "All models selected by AICc.  ",
  "Inference: Type III Wald ANOVA (car::Anova, type = 3).  ",
  "Post-hoc: Tukey-adjusted pairwise comparisons (emmeans).  ",
  "Compact letter display (CLD) always computed for the full interaction term ",
  "and for every additionally significant non-intercept term."
))

# =============================================================================
# SECTION 1: Individual-Level Analysis
# Model structure: Group × Treatment × Interval | Experimental unit: fish
# =============================================================================
.sec1_labels <- if (length(ind_model_results) > 0L) {
  vapply(names(ind_model_results),
         function(k) .res_label(ind_model_results[[k]], k), character(1L))
} else character(0L)

doc <- .add_section_header(doc, 1L, "Individual-Level Analysis",
  c("Model structure: Group × Treatment × Interval (three-way mixed model).",
    "Experimental unit: fish (one observation per fish per interval per trial).",
    paste0("Analysed intervals: ",
           paste(pretty_interval_name(as.character(ANALYSIS_TIMEPOINTS)),
                 collapse = ", "), "."),
    "Random effects: (1|trial_id), (1|trial_date); (1|fish_id_matched_trial) where idmatcher reliability qualifies.",
    "Interval modelled as factor and as numeric trend; AICc selects the better parameterisation.",
    if (length(.sec1_labels) > 0L)
      paste0("Indicators: ", paste(unique(.sec1_labels), collapse = "; "), ".")
    else "No individual-level interval results available in this run."
  ))

for (.key1 in names(ind_model_results)) {
  .res1 <- ind_model_results[[.key1]]
  if (is.null(.res1$fit)) next
  doc <- officer::body_add_break(doc)
  doc <- add_p(doc, .res_label(.res1, .key1), style = "heading 2")
  doc <- .render_result_block(doc, .res1)
}

# =============================================================================
# SECTION 2: Group Dynamics
# Model structure: Interval × Group × Treatment | Experimental unit: trial
# =============================================================================
.all_gd2 <- c(gd_model_results,
               if (!is.null(nnd_model_results)) list(GD_NND_cm = nnd_model_results))
.sec2_labels <- vapply(names(.all_gd2),
                       function(k) .res_label(.all_gd2[[k]], k), character(1L))

doc <- .add_section_header(doc, 2L, "Group Dynamics",
  c("Model structure: Interval × Group × Treatment (three-way mixed model).",
    "Experimental unit: trial (one session-level observation per interval per trial).",
    paste0("Analysed intervals: ",
           paste(pretty_interval_name(as.character(ANALYSIS_TIMEPOINTS)),
                 collapse = ", "), "."),
    "Random effects: (1|trial_id), (1|trial_date).",
    "Interval modelled as factor and as numeric trend; AICc selects the better parameterisation.",
    if (length(.sec2_labels) > 0L)
      paste0("Indicators: ", paste(unique(.sec2_labels), collapse = "; "), ".")
    else "No group-dynamics interval results available in this run."
  ))

for (.key2 in names(.all_gd2)) {
  .res2 <- .all_gd2[[.key2]]
  if (is.null(.res2$fit)) next
  doc <- officer::body_add_break(doc)
  doc <- add_p(doc, .res_label(.res2, .key2), style = "heading 2")
  doc <- .render_result_block(doc, .res2)
}

# =============================================================================
# SECTION 3: Interval-Aggregated Individual-Level Analysis
# Model structure: Group × Treatment | Experimental unit: fish
# =============================================================================
trial_agg_ind <- .load_trial_agg_ind()

.sec3_labels <- if (!is.null(trial_agg_ind) && length(trial_agg_ind) > 0L) {
  vapply(names(trial_agg_ind),
         function(nm) {
           r <- trial_agg_ind[[nm]]
           if (!is.null(r$label) && nzchar(r$label)) r$label
           else pretty_indicator_name(nm)
         }, character(1L))
} else character(0L)

doc <- .add_section_header(doc, 3L,
  "Interval-Aggregated Individual-Level Analysis",
  c("Model structure: Group × Treatment (two-way; interval factor removed).",
    "Experimental unit: fish (values averaged across Intervals 2 and 3 per fish per trial).",
    "Model: lm(y ~ Treatment × Group); (1|trial_date) added where not collinear with trial.",
    "Pseudoreplication-free sensitivity analysis complementing Section 1.",
    if (length(.sec3_labels) > 0L)
      paste0("Indicators: ", paste(unique(.sec3_labels), collapse = "; "), ".")
    else "STEP5 trial-aggregated results not available — run STEP5 before STEP6."
  ))

if (!is.null(trial_agg_ind) && length(trial_agg_ind) > 0L) {
  for (.ind3 in names(trial_agg_ind)) {
    .r3 <- trial_agg_ind[[.ind3]]
    if (is.null(.r3) || !is.data.frame(.r3$anova_table) ||
        nrow(.r3$anova_table) == 0L) next
    .lbl3 <- if (!is.null(.r3$label) && nzchar(.r3$label)) .r3$label
             else pretty_indicator_name(.ind3)
    if (is.null(.r3$model_type)) .r3$model_type <- "lm"
    doc <- officer::body_add_break(doc)
    doc <- add_p(doc, .lbl3, style = "heading 2")
    doc <- .render_agg_block(doc, .r3, .ind3)
  }
} else {
  doc <- add_p(doc,
    "Section 3 data not available: ensure STEP5 runs before STEP6 in the same pipeline session.")
}

# =============================================================================
# SECTION 4: Interval-Aggregated Group Dynamics
# Model structure: Group × Treatment | Experimental unit: trial
# =============================================================================
.sec4_ok <- exists("trial_agg_gd_results") &&
  is.list(trial_agg_gd_results) && length(trial_agg_gd_results) > 0L
.sec4_labels <- if (.sec4_ok) {
  vapply(names(trial_agg_gd_results), pretty_indicator_name, character(1L))
} else character(0L)

doc <- .add_section_header(doc, 4L,
  "Interval-Aggregated Group Dynamics",
  c("Model structure: Group × Treatment (two-way; interval factor removed).",
    "Experimental unit: trial (session-level group metric averaged across Intervals 2 and 3).",
    "Model: lm(y ~ Treatment × Group).",
    "Pseudoreplication-free sensitivity analysis complementing Section 2.",
    if (length(.sec4_labels) > 0L)
      paste0("Indicators: ", paste(unique(.sec4_labels), collapse = "; "), ".")
    else "Trial-aggregated group dynamics not available in this run."
  ))

if (.sec4_ok) {
  for (.ind4 in names(trial_agg_gd_results)) {
    .r4 <- trial_agg_gd_results[[.ind4]]
    if (is.null(.r4) || !is.data.frame(.r4$anova_table) ||
        nrow(.r4$anova_table) == 0L) next
    if (is.null(.r4$model_type)) .r4$model_type <- "lm"
    doc <- officer::body_add_break(doc)
    doc <- add_p(doc, pretty_indicator_name(.ind4), style = "heading 2")
    doc <- .render_agg_block(doc, .r4, .ind4)
  }
} else {
  doc <- add_p(doc,
    "Section 4 data not available: trial-aggregated group dynamics branch did not run.")
}

# =============================================================================
# SECTION 5: Collated Assumption Diagnostics Summary
# One row per indicator × block, collated from Sections 1–4.
# =============================================================================

# Helper: extract $diagnostics (or $diag_row) data.frames from a named result list,
# and prepend Block + Indicator columns for the combined view.
.collect_diag_section <- function(res_list, block_label) {
  if (is.null(res_list) || length(res_list) == 0L) return(NULL)
  rows <- lapply(names(res_list), function(nm) {
    r <- res_list[[nm]]
    if (is.null(r)) return(NULL)
    # Accept both field names: STEP6 uses $diagnostics; STEP5 per-indicator uses $diag_row
    diag_df <- if (!is.null(r$diagnostics) && is.data.frame(r$diagnostics) &&
                   nrow(r$diagnostics) > 0L) {
      r$diagnostics
    } else if (!is.null(r$diag_row) && is.data.frame(r$diag_row) &&
               nrow(r$diag_row) > 0L) {
      r$diag_row
    } else {
      return(NULL)
    }
    diag_df[["Block"]]     <- block_label
    diag_df[["Indicator"]] <- tryCatch(pretty_indicator_name(nm),
                                        error = function(e) nm)
    diag_df
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(NULL)
  tryCatch(do.call(rbind, rows), error = function(e) NULL)
}

# Collect from each of the four sections
.diag_s1 <- .collect_diag_section(
  if (exists("ind_model_results")) ind_model_results else list(),
  "1 — Individual (interval-level)")

.diag_s2 <- .collect_diag_section(
  c(if (exists("gd_model_results")) gd_model_results else list(),
    if (exists("nnd_model_results") && !is.null(nnd_model_results))
      list(GD_NND_cm = nnd_model_results) else list()),
  "2 — Group dynamics (interval-level)")

# Section 3: STEP5 trial-aggregated individual results
# trial_agg_ind was loaded at line 3115; reuse it to avoid a second disk read.
.diag_s3 <- if (exists("trial_agg_ind") && !is.null(trial_agg_ind))
  .collect_diag_section(trial_agg_ind, "3 — Individual (trial-aggregated)") else NULL

# Section 4: STEP6 trial-aggregated group dynamics
.diag_s4 <- if (exists("trial_agg_gd_results") && is.list(trial_agg_gd_results)) {
  .collect_diag_section(trial_agg_gd_results, "4 — Group dynamics (trial-aggregated)")
} else NULL

.diag_all <- tryCatch(
  do.call(rbind, Filter(Negate(is.null), list(.diag_s1, .diag_s2, .diag_s3, .diag_s4))),
  error = function(e) NULL
)

if (!is.null(.diag_all) && nrow(.diag_all) > 0L) {
  # Reorder: Block and Indicator columns first, then remaining diagnostic columns
  .lead_cols  <- intersect(c("Block", "Indicator"), names(.diag_all))
  .other_cols <- setdiff(names(.diag_all), .lead_cols)
  .diag_all   <- .diag_all[, c(.lead_cols, .other_cols), drop = FALSE]
  rownames(.diag_all) <- NULL

  doc <- .add_section_header(doc, 5L, "Assumption Diagnostics Summary",
    c("One row per indicator per analysis block.",
      "SW p = Shapiro-Wilk p on model residuals (post-fitting for aggregated lm).",
      "Levene p = Levene test for homoscedasticity.",
      "'not applicable' = non-Gaussian family (Beta/Gamma/Poisson-NB); no transformation test needed.",
      "Aggregated-branch diagnostics are post-hoc only — no automatic refitting applied.",
      paste0(nrow(.diag_all), " rows across ",
             length(unique(.diag_all[["Block"]])), " blocks.")
    ))
  doc <- add_tbl(doc, "Assumption diagnostics — all blocks", .diag_all, font_size = 7L)
  message("[DIAG] Collated diagnostics table: ", nrow(.diag_all), " rows across ",
          length(unique(.diag_all[["Block"]])), " blocks.")
} else {
  message("[DIAG] No diagnostics rows collected for Section 5 — skipped.")
}

# =============================================================================
# Save combined report
# =============================================================================
REPORT_PATH <- file.path(
  DIR_REPORTS,
  paste0("REPORT_combined_", TIMEPOINT_LABEL, "_", GROUP_LABEL, ".docx")
)
print(doc, target = REPORT_PATH)
message("Combined Word report: ", normalizePath(REPORT_PATH))

# -- Markdown report (group dynamics + individual-level) ----------------------
tryCatch({
  .gd_all <- c(
    gd_model_results,
    if (!is.null(nnd_model_results)) list(GD_NND_cm = nnd_model_results)
  )
  .ind_all <- if (exists("ind_model_results")) ind_model_results else list()
  .n_sess  <- if (exists("session_gd_summary")) nrow(session_gd_summary) else NA_integer_
  write_report_markdown_step6(.gd_all, .ind_all, .n_sess, OUT_DIR)
}, error = function(e) warning("STEP6 MD report failed: ", conditionMessage(e), call. = FALSE))


# =============================================================================
# 11) EXPORT TO GLOBAL ENVIRONMENT + MANIFEST
# =============================================================================
# Expose summary tables so downstream scripts can access them

assign("session_group_dynamics_summary", session_gd_summary, envir = .GlobalEnv)
if (!is.null(nnd_model_data))
  assign("nnd_group_dynamics_fish_tp",   nnd_model_data,      envir = .GlobalEnv)
if (!is.null(iid_fish_tp) && nrow(iid_fish_tp) > 0L)
  assign("iid_group_dynamics_fish_tp",   iid_fish_tp,         envir = .GlobalEnv)
if (!is.null(hull_session_tp) && nrow(hull_session_tp) > 0L)
  assign("hull_area_session_tp",         hull_session_tp,     envir = .GlobalEnv)

# Output manifest — lists every file created by this script
manifest <- data.frame(
  file = list.files(OUT_DIR, recursive = TRUE, full.names = FALSE),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(subfolder = dirname(file), filename = basename(file))

write.csv(manifest, file.path(OUT_DIR, "output_manifest.csv"), row.names = FALSE)

# Settings log for reproducibility
write.csv(
  data.frame(
    RUN_STAMP           = RUN_STAMP,
    USE_TP1             = USE_TP1,
    USE_ALL_GROUPS      = USE_ALL_GROUPS,
    ANALYSIS_TIMEPOINTS = paste(ANALYSIS_TIMEPOINTS, collapse = ","),
    GROUP_MODE          = ifelse(USE_ALL_GROUPS, "ALL_GROUPS",
                                 paste(GROUPS_RESTRICTED, collapse = ",")),
    stringsAsFactors = FALSE
  ),
  file.path(OUT_DIR, "run_settings.csv"),
  row.names = FALSE
)

.got_write_session_info(OUT_DIR)

# =============================================================================
# BH-adjusted p-values across STEP6 outcome family (informational)
# =============================================================================
tryCatch({
  p_rows_gd <- lapply(names(gd_model_results), function(nm) {
    r  <- gd_model_results[[nm]]
    at <- r$anova
    if (!is.data.frame(at) || nrow(at) == 0) return(NULL)
    p_col <- intersect(c("Pr(>F)", "Pr(>Chisq)", "p.value"), names(at))
    if (length(p_col) == 0) return(NULL)
    primary_terms <- c("Group:Treatment:timepoint", "Treatment:timepoint",
                       "Group:Treatment", "Group*Treatment*timepoint")
    hit <- intersect(primary_terms, rownames(at))
    if (length(hit) == 0) hit <- rownames(at)[nrow(at)]
    data.frame(outcome = r$col, term = hit[1],
               p_value = as.numeric(at[hit[1], p_col[1]]),
               stringsAsFactors = FALSE)
  })
  p_df_gd <- do.call(rbind, Filter(Negate(is.null), p_rows_gd))
  if (!is.null(p_df_gd) && nrow(p_df_gd) > 0) {
    p_df_gd$p_BH <- p.adjust(p_df_gd$p_value, method = "BH")
    write.csv(p_df_gd, file.path(DIR_ANOVA, "BH_adjusted_pvalues_STEP6.csv"), row.names = FALSE)
  }
}, error = function(e) warning("BH computation failed: ", conditionMessage(e), call. = FALSE))

message("GROUP DYNAMICS PIPELINE COMPLETE.  Outputs: ", normalizePath(OUT_DIR))
