# =============================================================================
# STEP 1.5 — idmatcher.ai cross-segment identity annotation
#
# Runs automatically after STEP1 (which produces master_bodylength_by_fish.csv)
# and before STEP2 (which joins body-length speed and activity status).
#
# What this step does:
#   1. Auto-detects the most recent idmatcher.ai Colab output folder.
#   2. Stages its files into STEP1_5_output/.
#   3. Calls annotate_idmatcher_scores_run() to label match quality.
#   4. Assigns stable trial-scoped fish_id_matched_trial identifiers
#      (format T<trial>_MF<n>; MF = matched fish identity).
#   5. Computes mean body length across seg2 + seg3 for each matched fish.
#   6. Derives a per-trial idmatcher_tier (conservative/reliable/borderline/unreliable).
#   7. Merges with the Hungarian matcher's hungarian_tier to produce a single
#      reliability_tier_<timestamp>.csv (used by STEP5/STEP6 and side-workflows).
#   8. Writes outputs and assigns key tables to GlobalEnv for downstream steps.
#
# GlobalEnv objects created:
#   idmatcher_id_map      – lookup (trial_id, segment, fish_in_session) →
#                           fish_id_matched_trial + body_length_mean_cm
#   idmatcher_bl_lookup   – per matched fish BL diagnostics
#   idmatcher_tier_df     – merged tier table (hungarian_tier + idmatcher_tier)
#   BL_MATCH_PATH         – updated to point to the new merged reliability_tier CSV
#                           (keeps reliability_tier alias = hungarian_tier for backward
#                           compatibility with side-workflows that use the default tier_col)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

t_start_1_5 <- proc.time()[["elapsed"]]

# ── Resolve required globals ───────────────────────────────────────────────────
if (!exists("PIPELINE_DATA_DIR", inherits = FALSE))
  stop("STEP 1.5: PIPELINE_DATA_DIR not found — run from 00_master_pipeline_GOT.R.")
if (!exists("IDMATCHER_COLAB_ROOT", inherits = FALSE))
  stop("STEP 1.5: IDMATCHER_COLAB_ROOT not found — add it to 00_master_pipeline_GOT.R config.")

bl_csv_path <- file.path(PIPELINE_DATA_DIR, "master_outputs", "morph",
                          "master_bodylength_by_fish.csv")
if (!file.exists(bl_csv_path))
  stop("STEP 1.5: master_bodylength_by_fish.csv not found at ", bl_csv_path,
       "\n  Run STEP1 first or check PIPELINE_DATA_DIR.")

# ── Source the annotation engine if not already loaded ───────────────────────
if (!exists("annotate_idmatcher_scores_run", inherits = FALSE)) {
  .annotate_path <- file.path(.pipeline_dir, "helpers", "annotate_idmatcher_scores.R")
  if (!file.exists(.annotate_path))
    stop("STEP 1.5: Cannot find annotate_idmatcher_scores.R at ", .annotate_path)
  ANNOTATE_IDMATCHER_CALLED_FROM_PIPELINE <- TRUE
  source(.annotate_path, local = FALSE)
}

# ── Create output directory ───────────────────────────────────────────────────
.step1_5_dir  <- .got_make_output_dir("STEP1_5")
.step1_5_ts   <- format(Sys.time(), "%Y%m%d_%H%M%S")

cat(strrep("-", 70), "\n")
cat("STEP 1.5 — idmatcher.ai identity annotation\n")
cat("Output dir: ", .step1_5_dir, "\n")
cat(strrep("-", 70), "\n\n")

# ── 1. Find latest Colab run ───────────────────────────────────────────────────
.colab_run_dir <- find_latest_idmatcher_run(root = IDMATCHER_COLAB_ROOT)
cat(sprintf("Latest Colab run: %s\n", basename(.colab_run_dir)))
cat(sprintf("  Has matrices/: %s\n", isTRUE(attr(.colab_run_dir, "has_matrices"))))

# ── 2. Stage files ────────────────────────────────────────────────────────────
.staged_inputs_dir   <- file.path(.step1_5_dir, "inputs")
.staged_matrices_dir <- file.path(.step1_5_dir, "matrices")
dir.create(.staged_inputs_dir, showWarnings = FALSE)

.corr_src  <- file.path(.colab_run_dir, "idmatcher_identity_correspondences.csv")
.corr_dest <- file.path(.staged_inputs_dir, "idmatcher_identity_correspondences.csv")
file.copy(.corr_src, .corr_dest, overwrite = TRUE)

.matrices_dir <- NULL
if (isTRUE(attr(.colab_run_dir, "has_matrices"))) {
  dir.create(.staged_matrices_dir, showWarnings = FALSE)
  .mat_files <- list.files(file.path(.colab_run_dir, "matrices"),
                            full.names = TRUE, pattern = "\\.csv$")
  if (length(.mat_files) > 0) {
    file.copy(.mat_files, .staged_matrices_dir, overwrite = TRUE)
    .matrices_dir <- .staged_matrices_dir
    cat(sprintf("  Staged %d matrix file(s).\n", length(.mat_files)))
  }
}

# ── 3. Run annotation engine ──────────────────────────────────────────────────
cat("\nRunning annotate_idmatcher_scores_run()...\n")
.annotated_raw <- annotate_idmatcher_scores_run(
  input_path   = .corr_src,
  bl_csv_path  = bl_csv_path,
  matrices_dir = if (!is.null(.matrices_dir)) file.path(.colab_run_dir, "matrices") else NULL,
  output_path  = NULL   # we write our own outputs below
)

# ── 4. Assign fish_id_matched_trial ───────────────────────────────────────────
cat("\nAssigning fish_id_matched_trial identifiers...\n")
.annotated <- assign_fish_id_matched_trial(.annotated_raw)

n_matched <- sum(!is.na(.annotated$fish_id_matched_trial))
n_total   <- nrow(.annotated)
cat(sprintf("  Matched fish: %d / %d rows received fish_id_matched_trial\n",
            n_matched, n_total))
cat(sprintf("  Unmatched   : %d rows have NA (rejected / no correspondence)\n",
            n_total - n_matched))

# ── 5. Recalculate body length ────────────────────────────────────────────────
cat("\nRecalculating mean body length across seg2 + seg3...\n")
.bl_lookup <- recalculate_body_length(.annotated, bl_csv_path)

cat(sprintf("  Mean BL computed for %d matched fish identities.\n",
            sum(!is.na(.bl_lookup$body_length_mean_cm))))
cat(sprintf("  BL source breakdown:\n"))
print(table(.bl_lookup$body_length_source, useNA = "ifany"))

# ── 6. Compute idmatcher tier ─────────────────────────────────────────────────
cat("\nComputing per-trial idmatcher reliability tier...\n")
.idmatcher_tier <- compute_idmatcher_reliability_tier(.annotated)
cat("  Tier distribution:\n")
print(table(.idmatcher_tier$idmatcher_tier))

# ── 7. Build idmatcher_id_map ─────────────────────────────────────────────────
cat("\nBuilding idmatcher_id_map (canonical seg × fish lookup)...\n")
.id_map <- build_idmatcher_id_map(.annotated, .bl_lookup)
cat(sprintf("  idmatcher_id_map: %d rows, %d trials\n",
            nrow(.id_map), n_distinct(.id_map$trial_id)))

# ── 8. Merge with Hungarian tier (if available) ───────────────────────────────
.hungarian_tier_df <- NULL

# Try to read the most recent Hungarian tier CSV
.bl_match_path_src <- if (exists("BL_MATCH_PATH", inherits = FALSE) &&
                           !is.null(BL_MATCH_PATH) && file.exists(BL_MATCH_PATH)) {
  BL_MATCH_PATH
} else {
  # Search in standard match BL output folder
  .match_bl_outdir <- file.path(getwd(), "match_BL_phyS_track_TP3_GOT_output")
  if (dir.exists(.match_bl_outdir)) {
    .bl_candidates <- list.files(.match_bl_outdir, pattern = "bl_reliability_tier.*\\.csv$",
                                  recursive = TRUE, full.names = TRUE)
    if (length(.bl_candidates) > 0)
      .bl_candidates[order(basename(.bl_candidates), decreasing = TRUE)][1]
    else NULL
  } else NULL
}

if (!is.null(.bl_match_path_src) && file.exists(.bl_match_path_src)) {
  .hungarian_raw <- tryCatch(
    read_csv(.bl_match_path_src, show_col_types = FALSE),
    error = function(e) NULL
  )
  if (!is.null(.hungarian_raw) && "trial_id" %in% names(.hungarian_raw) &&
      "reliability_tier" %in% names(.hungarian_raw)) {
    .hungarian_tier_df <- .hungarian_raw |>
      distinct(trial_id, .keep_all = TRUE) |>
      select(trial_id, hungarian_tier = reliability_tier)
    cat(sprintf("\nMerging Hungarian tier from: %s\n",
                basename(.bl_match_path_src)))
  }
}

if (is.null(.hungarian_tier_df)) {
  message("STEP 1.5: No Hungarian tier CSV found — merged CSV will contain only idmatcher_tier.")
  message("  Run 'match BL phyS_track TP3 GOT 11.03.R' first and set BL_MATCH_PATH.")
}

# Merge
.merged_tier <- .idmatcher_tier |>
  left_join(
    if (!is.null(.hungarian_tier_df)) .hungarian_tier_df
    else tibble(trial_id = integer(0), hungarian_tier = character(0)),
    by = "trial_id"
  ) |>
  mutate(
    # backward-compatibility alias: side-workflows using default tier_col="reliability_tier"
    # will get the hungarian_tier (unchanged from pre-integration behaviour)
    reliability_tier = coalesce(hungarian_tier, NA_character_)
  )

# ── 9. Write outputs ──────────────────────────────────────────────────────────

# 9a. Canonical id_map
.id_map_path <- file.path(.step1_5_dir,
                           paste0("idmatcher_id_map_", .step1_5_ts, ".csv"))
write_csv(.id_map, .id_map_path)

# 9b. BL diagnostics
.bl_diag_path <- file.path(.step1_5_dir,
                            paste0("body_length_diagnostics_", .step1_5_ts, ".csv"))
write_csv(.bl_lookup, .bl_diag_path)

# 9c. Merged tier CSV (the new BL_MATCH_PATH)
.tier_path <- file.path(.step1_5_dir,
                         paste0("reliability_tier_", .step1_5_ts, ".csv"))
write_csv(.merged_tier, .tier_path)

# 9d. Annotated table (full diagnostics)
.annotated_path <- file.path(.step1_5_dir,
                              paste0("idmatcher_annotated_", .step1_5_ts, ".csv"))
write_csv(.annotated, .annotated_path)

# 9e. Run settings + session info
.got_write_run_settings(
  .step1_5_dir,
  step              = "STEP1_5_annotate_idmatcher",
  colab_run_folder  = basename(.colab_run_dir),
  idmatcher_min_tier = if (exists("IDMATCHER_MIN_TIER")) IDMATCHER_MIN_TIER else "reliable",
  n_trials_processed = n_distinct(.id_map$trial_id),
  n_matched_fish    = sum(!is.na(.id_map$fish_id_matched_trial)),
  idmatcher_id_map_path  = .id_map_path,
  tier_csv_path     = .tier_path
)
.got_write_session_info(.step1_5_dir)

cat(sprintf("\nOutputs written to: %s\n", .step1_5_dir))
cat(sprintf("  idmatcher_id_map         : %s\n", basename(.id_map_path)))
cat(sprintf("  body_length_diagnostics  : %s\n", basename(.bl_diag_path)))
cat(sprintf("  reliability_tier (merged): %s\n", basename(.tier_path)))

# ── 10. Assign to GlobalEnv ───────────────────────────────────────────────────
idmatcher_id_map    <<- .id_map
idmatcher_bl_lookup <<- .bl_lookup
idmatcher_tier_df   <<- .merged_tier

# Update BL_MATCH_PATH so all downstream consumers (STEP5, side-workflows) get the merged CSV
BL_MATCH_PATH <<- .tier_path
cat(sprintf("\nBL_MATCH_PATH updated to: %s\n", basename(.tier_path)))

# ---- R4: stable canonical copy at master_outputs/reliability/ ---------------
# Side workflows (and STEP 5/6) can fall back to this fixed path when
# BL_MATCH_PATH is not explicitly set in the R session.
# Note: if STEP 1.6 also runs, it will overwrite this file with its own merged
# tier. Both writers point BL_MATCH_PATH to the same canonical path, so the
# fallback always finds the most recent producer's output.
.bl_pdd <- if (exists("PIPELINE_DATA_DIR", envir = .GlobalEnv, inherits = FALSE)) {
  get("PIPELINE_DATA_DIR", envir = .GlobalEnv)
} else {
  dirname(dirname(.step1_5_dir))   # heuristic fallback
}
.bl_canonical_dir  <- file.path(.bl_pdd, "master_outputs", "reliability")
dir.create(.bl_canonical_dir, recursive = TRUE, showWarnings = FALSE)
.bl_canonical_path <- file.path(.bl_canonical_dir, "reliability_tier_LATEST.csv")
file.copy(.tier_path, .bl_canonical_path, overwrite = TRUE)
cat(sprintf("R4: wrote canonical copy: %s\n", .bl_canonical_path))

# Point BL_MATCH_PATH to the canonical (stable) path so re-runs in fresh R
# sessions can find it without STEP 1.5 having been re-run.
BL_MATCH_PATH <<- .bl_canonical_path
cat(sprintf("R4: BL_MATCH_PATH set to canonical: %s\n", .bl_canonical_path))
rm(.bl_pdd, .bl_canonical_dir, .bl_canonical_path)

# Log
.elapsed_1_5 <- round(proc.time()[["elapsed"]] - t_start_1_5, 1)
.got_append_pipeline_log(1.5, "STEP1_5_annotate_idmatcher", .step1_5_dir, .elapsed_1_5)

cat(sprintf("\nSTEP 1.5 complete in %.1f s\n", .elapsed_1_5))
cat(strrep("-", 70), "\n\n")
