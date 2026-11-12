# =============================================================================
# MASTER PIPELINE — GOT ACTIVITY ANALYSIS
# =============================================================================
# This script runs the full analysis pipeline from start to finish.
# Instead of opening and running each script one by one, you just run this file.
#
# THE STEPS:
#   1  Parse raw idtracker tracking data into a clean fish-by-frame table
#   2  Add activity labels (moving/immobile) and body-length speed to the table
#   3  [REMOVED — superseded by Step 5] slot kept so step numbering is stable
#   4  Make scatter plots and timepoint-trend graphs → saved as PowerPoint files
#   5  Run the primary statistics pipeline → tables + Word report
#   6  Group dynamics pipeline → polarisation, NND, IID, zone sync, centroid speed
#
#   (1.5 and 1.6, run between steps 1 and 2: identity matching across intervals
#    with idmatcher.ai results, and body-length matching to the sampled fish.)
#   Steps 4 and 5 are exploratory; the statistics reported in the article come
#   from 05_statistics/behaviour_models.R, which reads the step 6 outputs.
#
# HOW TO USE:
#   1. Set the working directory to the repository root and put the input data
#      in data/ (see README.md, "Input data").
#   2. Set START_FROM_STEP if you have already run some steps and don't want
#      to repeat them (e.g. set to 2 if step 1 is already done).
#   3. Run: source("02_trajectory_pipeline/00_run_trajectory_pipeline.R")
# =============================================================================


# =============================================================================
# ==== CONFIG — THE ONLY SECTION YOU NEED TO EDIT ==============================
# =============================================================================

# Folder with one idtracker.ai session folder per trial x interval
# (relative to the repository root). Use forward slashes / even on Windows.
PIPELINE_DATA_DIR <- "data/tracking_sessions"

# ---- Metadata paths ----
# Paths to external metadata workbooks and the Word template.
# Edit here once; child scripts read from GlobalEnv with an inline fallback.
TRIAL_SUMMARY_XLSX  <- "data/trials summary GOT.xlsx"
SAMPLING_DATA_XLSX  <- "data/sampling data GOT.xlsx"
REFERENCE_DOCX_PATH <- "data/reference_a4_portrait.docx"

# Which step to start from.
# Set to 1 to run everything from scratch (default).
# Set to 2 if step 1 (idtracker parsing) is already done — master_outputs exists on disk.
# Set to 3 skips to step 4 (step 3 is deprecated; its slot is kept so step numbering is stable).
# Set to 4 to re-run only the graphs and stats (master_fish_by_frame must be in R).
# Set to 5 to re-run only the primary stats script (fish_activity_summary must be in R).
# Set to 6 to re-run only the group dynamics script (master_fish_by_frame + group tables must be in R).
if (!exists("START_FROM_STEP")) START_FROM_STEP <- 1

# Minimum body-length match reliability tier required to include the fish-level
# random effect (1|fish_id) / (1|Sample_ID) in model candidate sets.
# Tiers (in ascending reliability): "unreliable" < "borderline" < "reliable" < "conservative"
# Default "reliable": fish from borderline/unreliable BL matches only get (1|trial_id).
MIN_BL_RELIABILITY_TIER <- "reliable"

# ---- STEP 1.5: idmatcher.ai cross-segment identity annotation ----
# Set RUN_STEP_1_5 = FALSE to skip (e.g. if no Colab run exists yet).
RUN_STEP_1_5         <- TRUE
# Root folder that contains the YYYYMMDD_HHMMSS Colab output subfolders.
IDMATCHER_COLAB_ROOT <- file.path("data", "idmatcher",
                                   "idmatcherai outputs from colab")
# Minimum idmatcher tier for a trial to use (1|fish_id_matched_trial) in STEP5/6.
IDMATCHER_MIN_TIER   <- "reliable"

# ---- STEP 1.6: Hungarian body-length matching (R5: promoted from physiology) ----
# Set RUN_STEP_1_6 = FALSE to skip (e.g. when sampling data isn't available).
# When enabled, writes the BL reliability-tier CSV used by STEP 5/6 + plasma/
# monoamine/aggression side workflows for random-effect candidate gating.
RUN_STEP_1_6         <- TRUE

# Set to TRUE to stop the pipeline immediately if any step throws an error.
# Set to FALSE to attempt to continue to the next step even after a failure
# (not recommended — later steps may run on incomplete data).
STOP_ON_ERROR <- TRUE

# ---- Parallelism (OPT-1 / OPT-2 / OPT-4) ----
# Number of parallel workers for each step.
# 0L = serial (safe default, identical output to original pipeline).
# Set > 0 to enable parallelism; recommended starting value = detectCores() - 1.
STEP1_PARALLEL_WORKERS <- 0L   # STEP 1: per-session idtracker loop
STEP4_PARALLEL_WORKERS <- 6L   # STEP 4: per-plot ggplot+DML generation
STEP5_PARALLEL_WORKERS <- 6L   # STEP 5: per-indicator LMM fits

# ---- Scope filters (OPT-0) ----
# These switches propagate to every downstream step via GlobalEnv.
# STEP 1 uses them to skip sessions before the per-session loop starts
# (the OPT-0 filter message will confirm how many sessions are dropped).
# STEP 2, 4, 5, and 6 use them to filter rows from the data tables.

# Set to FALSE to drop Segment 1 / Timepoint 1 sessions throughout the pipeline.
# Set to TRUE to include them (useful for acclimation-period analysis).
KEEP_TP1 <- FALSE

# Set to TRUE to exclude Morphine / Naloxone trial groups throughout the pipeline.
# When TRUE, STEP 1 skips those sessions entirely (fastest); downstream steps
# also filter their data tables as a safety net.
DROP_PHARMA <- TRUE

# The group names that identify pharmacological intervention trials.
# Matched case-insensitively against the 'Group' column in TRIAL_SUMMARY_XLSX.
PHARMA_GROUPS <- c("Morphine", "Naloxone")


# =============================================================================
# ==== INTERNAL SETUP — do not edit below this line ============================
# =============================================================================

# Resolve this script's own directory NOW, at the moment the file is sourced.
# This is done at the top level (not inside a function) so that the call-stack
# frame that holds the file path is still accessible.
.pipeline_dir <- tryCatch({
  # Walk the call-stack frames and find the one that has an 'ofile' attribute,
  # which R sets when a script is being executed via source().
  frames     <- sys.frames()
  ofile_envs <- Filter(function(f) exists("ofile", envir = f, inherits = FALSE), frames)
  if (length(ofile_envs) > 0) {
    # Use the innermost frame that has ofile (the most recently sourced script)
    normalizePath(
      dirname(get("ofile", envir = ofile_envs[[length(ofile_envs)]])),
      winslash = "/", mustWork = FALSE
    )
  } else {
    getwd()   # fallback when running line-by-line instead of via source()
  }
}, error = function(e) getwd())

# Source shared output-management helpers so they are in GlobalEnv before any
# step script runs. Defensive: silently skip if the file cannot be found.
.helpers_path <- file.path(.pipeline_dir, "helpers/output_helpers.R")
if (file.exists(.helpers_path)) {
  source(.helpers_path, local = FALSE)
} else {
  message("Note: 00_output_helpers.R not found at ", .helpers_path,
          " — helper functions will not be available.")
}

.stats_helpers_path <- file.path(.pipeline_dir, "helpers/stats_helpers.R")
if (file.exists(.stats_helpers_path)) {
  source(.stats_helpers_path, local = FALSE)
} else {
  message("Note: 00_stats_helpers.R not found at ", .stats_helpers_path,
          " — statistical helper functions will not be available.")
}

.idmatcher_helpers_path <- file.path(.pipeline_dir, "helpers/idmatcher_helpers.R")
if (file.exists(.idmatcher_helpers_path)) {
  source(.idmatcher_helpers_path, local = FALSE)
} else {
  message("Note: 00_idmatcher_helpers.R not found at ", .idmatcher_helpers_path,
          " — STEP 1.5 idmatcher functions will not be available.")
}

# Centralised pretty-label dictionary (R3) — sourced once here so that every
# step + side-workflow script can rely on `pretty_indicator_name()`,
# `GOT_INDICATOR_LABELS`, etc. being available in `.GlobalEnv` without each
# script re-implementing its own defensive fallback block.
.pretty_labels_path <- file.path(.pipeline_dir, "helpers/pretty_labels.R")
if (file.exists(.pretty_labels_path)) {
  source(.pretty_labels_path, local = FALSE)
} else {
  message("Note: pretty_labels.R not found at ", .pretty_labels_path,
          " — pretty-label functions will not be available.")
}

# Transformation-diagnostics helpers — provides run_diag_for_indicator(),
# build_diag_table_row(), run_agg_diag() used by STEP5 and STEP6.
# Sourced after 00_stats_helpers so its SW/Levene wrappers are available.
.transdiag_path <- file.path(.pipeline_dir, "helpers/transformation_diagnostics.R")
if (file.exists(.transdiag_path)) {
  source(.transdiag_path, local = FALSE)
} else {
  message("Note: transformation_diagnostics.R not found at ",
          .transdiag_path, " — transformation diagnostics will not be available.")
}

# CLD helpers — provides remap_cld_letters(), all_cld_same(),
# extract_cld_for_plot(), add_cld_to_scatter(), build_remapped_report_tbl(),
# compute_agg_cld_df() used by STEP5 and STEP6 for graph annotation and reports.
.cld_helpers_path <- file.path(.pipeline_dir, "helpers/cld_helpers.R")
if (file.exists(.cld_helpers_path)) {
  source(.cld_helpers_path, local = FALSE)
} else {
  message("Note: cld_helpers.R not found at ",
          .cld_helpers_path, " — CLD annotation helpers will not be available.")
}


# Helper: run one step script, with timing, skip logic, and error handling.
.pipeline_source <- function(step_number, script_name, start_from, stop_on_error) {

  # --- Skip if below the chosen start step ---
  if (step_number < start_from) {
    cat(strrep("-", 70), "\n")
    cat("SKIPPED (START_FROM_STEP=", start_from, "): STEP ", step_number,
        " — ", script_name, "\n", sep = "")
    cat(strrep("-", 70), "\n\n")
    return(invisible(NULL))
  }

  # --- Locate the script file ---
  # Primary location: same folder as this master script (resolved at load time above)
  script_path_primary   <- file.path(.pipeline_dir, script_name)
  # Secondary location: 02_trajectory_pipeline/ subfolder of the working directory
  script_path_secondary <- file.path(getwd(), "02_trajectory_pipeline", script_name)

  if (file.exists(script_path_primary)) {
    script_path <- script_path_primary
  } else if (file.exists(script_path_secondary)) {
    script_path <- script_path_secondary
  } else {
    # Neither location worked — show both paths searched so the user can diagnose
    stop(
      "Script not found: ", script_name,
      "\n  Searched (1): ", script_path_primary,
      "\n  Searched (2): ", script_path_secondary,
      call. = FALSE
    )
  }

  # --- Print a visible header before running ---
  cat(strrep("=", 70), "\n")
  cat("STEP ", step_number, " — ", script_name, "\n", sep = "")
  cat(strrep("=", 70), "\n\n")

  # --- Record start time so we can report how long the step took ---
  t_start <- proc.time()[["elapsed"]]

  # --- Run the step script ---
  # local = FALSE (the default) means all objects the script creates go into
  # the global R environment so later steps can access them.
  result <- tryCatch(
    source(script_path, local = FALSE, echo = FALSE, verbose = FALSE),
    error = function(e) e   # capture any error instead of crashing immediately
  )

  elapsed <- round(proc.time()[["elapsed"]] - t_start, 1)

  # --- Report outcome ---
  if (inherits(result, "error")) {
    cat("\n", strrep("!", 70), "\n", sep = "")
    cat("ERROR in STEP ", step_number, " — ", script_name, "\n", sep = "")
    cat("Message: ", conditionMessage(result), "\n", sep = "")
    cat("Elapsed: ", elapsed, " s\n", sep = "")
    cat(strrep("!", 70), "\n\n", sep = "")

    if (isTRUE(stop_on_error)) {
      stop(
        "Pipeline halted at step ", step_number, " due to the error above.\n",
        "Fix the issue and re-run with START_FROM_STEP = ", step_number, ".",
        call. = FALSE
      )
    } else {
      warning("Step ", step_number, " failed but STOP_ON_ERROR = FALSE; continuing.",
              call. = FALSE)
    }
  } else {
    cat("\n", strrep("-", 70), "\n", sep = "")
    cat("DONE: STEP ", step_number, " — ", script_name,
        "  [", elapsed, " s]\n\n", sep = "")

    if (exists(".got_append_pipeline_log", mode = "function")) {
      tryCatch(
        .got_append_pipeline_log(step_number, script_name, NA_character_, elapsed),
        error = function(e) warning("Pipeline log write failed: ", conditionMessage(e),
                                    call. = FALSE)
      )
    }
  }

  invisible(NULL)
}


# =============================================================================
# ==== VALIDATE CONFIG BEFORE RUNNING ANYTHING =================================
# =============================================================================

# Stop immediately if the user forgot to set a real path
if (!nzchar(trimws(PIPELINE_DATA_DIR)) ||
    PIPELINE_DATA_DIR == "C:/path/to/your/session/folder") {
  stop(
    "PIPELINE_DATA_DIR has not been set.\n",
    "Open 00_master_pipeline_GOT.R and replace the placeholder with your real folder path.",
    call. = FALSE
  )
}

# Normalize slashes and check the folder actually exists on disk
PIPELINE_DATA_DIR <- normalizePath(PIPELINE_DATA_DIR, winslash = "/", mustWork = FALSE)

if (!dir.exists(PIPELINE_DATA_DIR)) {
  stop(
    "PIPELINE_DATA_DIR does not exist on disk:\n  ", PIPELINE_DATA_DIR, "\n",
    "Check the path for typos and make sure the drive is connected.",
    call. = FALSE
  )
}

# Stop if START_FROM_STEP is not a valid number (must be 1 through 6)
if (!is.numeric(START_FROM_STEP) || !START_FROM_STEP %in% 1:6) {
  stop("START_FROM_STEP must be a whole number between 1 and 6.", call. = FALSE)
}

# Warn the user if they are skipping steps but the required objects are missing.
# This catches the common mistake of starting from step 3 after a fresh R session.
if (START_FROM_STEP >= 3 &&
    !exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
  warning(
    "START_FROM_STEP = ", START_FROM_STEP,
    " but 'master_fish_by_frame' is not in the R environment.\n",
    "Steps 3 and 4 both require this object. ",
    "Either run from step 2, or load the object from disk first.",
    call. = FALSE
  )
}

if (START_FROM_STEP >= 5 &&
    !exists("fish_activity_summary", envir = .GlobalEnv, inherits = FALSE)) {
  warning(
    "START_FROM_STEP = ", START_FROM_STEP,
    " but 'fish_activity_summary' is not in the R environment.\n",
    "Step 5 requires this object. Either run from step 4, or load it from disk first.",
    call. = FALSE
  )
}

if (START_FROM_STEP >= 6 &&
    !exists("master_fish_by_frame", envir = .GlobalEnv, inherits = FALSE)) {
  warning(
    "START_FROM_STEP = 6 but 'master_fish_by_frame' is not in the R environment.\n",
    "Step 6 requires master_fish_by_frame plus the group dynamics tables from Step 1.\n",
    "Either run from step 1, or load the required objects from disk first.",
    call. = FALSE
  )
}

# Resolve metadata paths with a local fallback: if a path does not exist on
# disk, resolve_dataset_path() searches data/ for a file with the same basename.
TRIAL_SUMMARY_XLSX  <- resolve_dataset_path(TRIAL_SUMMARY_XLSX)
SAMPLING_DATA_XLSX  <- resolve_dataset_path(SAMPLING_DATA_XLSX)
REFERENCE_DOCX_PATH <- resolve_dataset_path(REFERENCE_DOCX_PATH)

# Make PIPELINE_DATA_DIR available globally so each step script can find the data
# (this replaces the clipboard paste that each step script normally requires)
assign("PIPELINE_DATA_DIR",       PIPELINE_DATA_DIR,       envir = .GlobalEnv)
assign("TRIAL_SUMMARY_XLSX",      TRIAL_SUMMARY_XLSX,      envir = .GlobalEnv)
assign("SAMPLING_DATA_XLSX",      SAMPLING_DATA_XLSX,      envir = .GlobalEnv)
assign("REFERENCE_DOCX_PATH",     REFERENCE_DOCX_PATH,     envir = .GlobalEnv)
assign("MIN_BL_RELIABILITY_TIER", MIN_BL_RELIABILITY_TIER, envir = .GlobalEnv)
# OPT-0 scope-filter flags — exposed so that every sourced step script and any
# side workflow run in the same R session can read them from GlobalEnv.
assign("KEEP_TP1",      KEEP_TP1,      envir = .GlobalEnv)
assign("DROP_PHARMA",   DROP_PHARMA,   envir = .GlobalEnv)
assign("PHARMA_GROUPS", PHARMA_GROUPS, envir = .GlobalEnv)

# Print a summary of what is about to run
cat(strrep("=", 70), "\n")
cat("MASTER PIPELINE — GOT\n")
cat("Data dir:        ", PIPELINE_DATA_DIR, "\n")
cat("Starting at step:", START_FROM_STEP, "\n")
cat("Stop on error:   ", STOP_ON_ERROR, "\n")
cat("KEEP_TP1:        ", KEEP_TP1, "  (FALSE = skip Segment1 sessions)\n")
cat("DROP_PHARMA:     ", DROP_PHARMA, "  (TRUE  = skip Morphine/Naloxone)\n")
cat(strrep("=", 70), "\n\n")

# Record overall pipeline start time
.pipeline_start_time <- proc.time()[["elapsed"]]


# =============================================================================
# ==== RUN THE STEPS IN ORDER ==================================================
# =============================================================================

# Step 1: Read raw idtracker output and build the master fish-by-frame table
.pipeline_source(1, "01_preprocess_trajectories.R",  START_FROM_STEP, STOP_ON_ERROR)

# Step 1.5: Annotate idmatcher.ai cross-segment identity matches; assign stable
#           fish_id_matched_trial; compute mean body length across seg2 + seg3;
#           merge with Hungarian tier into a single reliability_tier CSV.
if (isTRUE(RUN_STEP_1_5) && START_FROM_STEP <= 1) {
  .pipeline_source(1.5, "01_5_annotate_identity_matches.R", 1, STOP_ON_ERROR)
}

# Step 1.6 (R5): Hungarian body-length matching of sampled fish IDs to tracked
#                fish IDs. Promoted from 02_side_workflows/physiology/. Writes
#                bl_reliability_tier_<TS>.csv AND a stable canonical copy at
#                master_outputs/reliability/reliability_tier_LATEST.csv (R4),
#                and assigns BL_MATCH_PATH in .GlobalEnv. Set RUN_STEP_1_6 =
#                FALSE to skip (e.g. when sampling data isn't available yet).
if (!exists("RUN_STEP_1_6")) RUN_STEP_1_6 <- TRUE
if (isTRUE(RUN_STEP_1_6) && START_FROM_STEP <= 1) {
  .pipeline_source(1.6, "01_6_match_body_length.R", 1, STOP_ON_ERROR)
}

# Step 2: Attach activity status (moving/immobile) and body-length speed to
#         the table; export enriched master_fish_by_frame + diagnostic CSVs
.pipeline_source(2, "02_individual_indicators.R",    START_FROM_STEP, STOP_ON_ERROR)

# Step 3: REMOVED (R9). The legacy STEP 3 script
# (03_activity_analysis_DEPRECATED.R) was superseded by Step 5: every outcome
# it modelled plus nine more is covered by Step 5 with a pre-specified full
# Treatment × Group × timepoint fixed structure, mandatory (1|trial_id)
# random effect, and a parallel trial-aggregated branch. Audit confirmed no
# unique business logic in the deprecated script that is not also produced
# (or superseded) by Step 5; the script file has been removed from
# efficient_scripts/. For backward compatibility, START_FROM_STEP == 3 is
# silently treated as "jump to Step 4".

# Free large STEP2 objects before STEP4 to reclaim RAM.
# STEP4 only needs master_fish_by_frame (loaded from .GlobalEnv or disk).
# The activity_TS*, speed_by_frame, position_by_frame, zone_by_frame objects
# from STEP2 are no longer needed and together can occupy 10-20 GB of RAM.
if (START_FROM_STEP <= 2) {
  .step2_large <- c(
    "activity_TS0.5", "activity_TS1", "activity_TS1.5", "activity_TS2",
    "speed_by_frame", "position_by_frame", "zone_by_frame", "group_by_frame",
    "activity_base", "mf", "df_act2", "ppt_dist", "ppt",
    "fish_summaries", "distance_by_fish"
  )
  .to_rm <- intersect(.step2_large, ls(envir = .GlobalEnv))
  if (length(.to_rm) > 0) {
    rm(list = .to_rm, envir = .GlobalEnv)
    cat("Freed STEP2 objects before STEP4:", paste(.to_rm, collapse = ", "), "\n")
  }
  invisible(gc(verbose = FALSE))
}

# Step 4: Build fish-level summary table; make scatter + timepoint-trend plots
#         for every activity indicator; save as PowerPoint files
.pipeline_source(4, "04_exploratory_graphs.R",    START_FROM_STEP, STOP_ON_ERROR)

# Step 5: Run the second statistics pipeline on fish_activity_summary;
#         model selection, ANOVA, post-hoc tests, Word report + diagnostic plots
.pipeline_source(5, "05_exploratory_individual_stats.R",     START_FROM_STEP, STOP_ON_ERROR)

# Step 6: Compute group-level behavioural indicators (polarization, NND, IID,
#         zone synchrony, centroid speed, collective motion events, leadership);
#         run LMM stats + Tukey post-hoc; save scatter + trend PPTs and a Word report
.pipeline_source(6, "06_school_indicators.R", START_FROM_STEP, STOP_ON_ERROR)


# =============================================================================
# ==== FINISHED ================================================================
# =============================================================================

.total_elapsed <- round(proc.time()[["elapsed"]] - .pipeline_start_time, 1)

cat(strrep("=", 70), "\n")
cat("MASTER PIPELINE COMPLETE\n")
cat("Total elapsed:   ", .total_elapsed, " s\n", sep = "")
cat("Output folders created inside: ", getwd(), "\n")
cat(strrep("=", 70), "\n")
